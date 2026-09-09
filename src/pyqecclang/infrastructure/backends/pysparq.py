"可选 PySparQ 执行适配器。RIR 的根寄存器对应原生整数寄存器。"

from __future__ import annotations

import cmath
from threading import Lock

from pyqecclang.infrastructure.execution import (
    LocalEnter,
    LocalExit,
    RegisterState,
    check_memory,
    events,
    gate_matrix,
)
from pyqecclang.infrastructure.ir import Load, ValidationError
from pyqecclang.infrastructure.native import NativeContext, NativeSite
from pyqecclang.infrastructure.validation import locations, validate

_LOCK = Lock()


def run_pysparq(
    program,
    memory=None,
    *,
    max_steps=1_000_000,
    max_states=65536,
    native_registry=None,
    report=None,
):
    program = validate(program, require_closed=native_registry is None)
    native_modules = frozenset() if native_registry is None else native_registry.matching(program)
    if native_registry is not None and native_registry.missing(program):
        raise ValidationError("PySparQ 缺少实现：" + ", ".join(native_registry.missing(program)))
    report = {} if report is None else report
    report.update(native_calls=0, gate_events=0, native_labels=[], correctness="pending")
    memories = check_memory(program, memory)
    if any(r.type.address_width > 20 for r in program.main.resources):
        raise ValidationError("当前适配器拒绝物化超过 2^20 项的 QRAM")
    # 在接管全局注册表前完成展开预算检查。
    from pyqecclang.infrastructure.execution import expanded_steps

    if expanded_steps(program, max_steps, native_modules) > max_steps:
        raise ValidationError("PySparQ 执行超过展开预算")
    try:
        import pysparq as ps
    except ImportError as exc:
        raise ValidationError("PySparQ 执行需要已安装 pysparq 的环境") from exc

    with _LOCK:
        if ps.System.get_activated_register_size():
            raise ValidationError("PySparQ 的全局寄存器表非空；请先结束已有模拟")
        try:
            state = ps.SparseState()
            names = {reg.name: f"pyqec_{i}" for i, reg in enumerate(program.main.registers)}
            kinds = {
                "bits": ps.StateStorageType.General,
                "uint": ps.StateStorageType.UnsignedInteger,
                "sint": ps.StateStorageType.SignedInteger,
                "rational": ps.StateStorageType.Rational,
            }
            for reg in program.main.registers:
                if reg.type.width:
                    ps.AddRegister(names[reg.name], kinds[reg.type.kind], reg.type.width)(state)
            qrams = {}
            for res in program.main.resources:
                data = [0] * (1 << res.type.address_width)
                for address, value in memories[res.name].items():
                    data[address] = value
                qrams[res.name] = ps.QRAMCircuit_qutrit(
                    res.type.address_width, res.type.data_width, data
                )

            def apply(operator, controls=()):
                if controls:
                    operator.conditioned_by_bit(list(controls))
                operator(state)
                if len(state.basis_states) > max_states:
                    raise ValidationError("PySparQ 执行超过稀疏态数量预算")

            def actual(ref):
                return [(names[name], bit) for name, bit in locations(ref)]

            for node, conditions, inverse in events(
                program, max_steps=max_steps, native_modules=native_modules
            ):
                if isinstance(node, LocalEnter):
                    names[node.name] = "pyqec_" + node.name
                    ps.AddRegister(names[node.name], kinds[node.type.kind], node.type.width)(state)
                    continue
                if isinstance(node, LocalExit):
                    native_name = names[node.name]
                    rid = ps.System.get_id(native_name)
                    if any(
                        int(basis.get(rid).value) & ((1 << node.width) - 1)
                        for basis in state.basis_states
                    ):
                        raise ValidationError(f"PySparQ 局部工作区未复净：{node.name}")
                    ps.RemoveRegister(native_name)(state)
                    del names[node.name]
                    continue
                controls, zeros = [], []
                for ref, expected in conditions:
                    for bit_index, pair in enumerate(actual(ref)):
                        controls.append(pair)
                        if not ((expected >> bit_index) & 1):
                            zeros.append(pair)
                for pair in zeros:
                    ps.Xgate_Bool(*pair)(state)
                if isinstance(node, NativeSite):
                    entry = native_registry.entries[node.module.name]
                    context = NativeContext(ps, node, names, qrams, memories)
                    untouched = (
                        ps.split_systems(
                            state,
                            [],
                            [],
                            [(ps.System.get_id(name), bit) for name, bit in controls],
                            [],
                        )
                        if controls
                        else None
                    )
                    try:
                        if state.size():
                            operator = entry.factory(context)
                            (operator.dag if inverse else operator)(state)
                            if len(state.basis_states) > max_states:
                                raise ValidationError("原生算子超过稀疏态预算")
                    finally:
                        if untouched is not None:
                            ps.combine_systems(state, untouched)
                    report["native_calls"] += 1
                    if entry.label not in report["native_labels"]:
                        report["native_labels"].append(entry.label)
                elif isinstance(node, Load):
                    # 视图通过可逆 XOR 复制到临时整数寄存器；载入后完整反算。
                    a, d = "pyqec_tmp_address", "pyqec_tmp_data"
                    ps.AddRegister(a, ps.StateStorageType.General, node.address.width)(state)
                    ps.AddRegister(d, ps.StateStorageType.General, node.data.width)(state)
                    for i, pair in enumerate(actual(node.address)):
                        apply(ps.Xgate_Bool(a, i), (pair,))
                    query = ps.QRAMLoad(qrams[node.resource], a, d)
                    query(state)
                    for i, pair in enumerate(actual(node.data)):
                        apply(ps.Xgate_Bool(*pair), tuple(controls) + ((d, i),))
                    query(state)
                    for i, pair in reversed(list(enumerate(actual(node.address)))):
                        apply(ps.Xgate_Bool(a, i), (pair,))
                    ps.RemoveRegister(d)(state)
                    ps.RemoveRegister(a)(state)
                elif node.op == "gphase":
                    theta = node.angle * (-1 if inverse else 1)
                    if controls:
                        apply(ps.Phase_Bool(*controls[-1], theta), controls[:-1])
                    else:
                        root = next(r for r in program.main.registers if r.type.width)
                        phase = cmath.exp(1j * theta)
                        apply(ps.Rot_Bool(names[root.name], 0, [phase, 0j, 0j, phase]))
                elif node.op in {"xor", "swap"}:
                    for a, b in zip(
                        actual(node.operands[0]), actual(node.operands[1]), strict=True
                    ):
                        apply(ps.Xgate_Bool(*b), tuple(controls) + (a,))
                        if node.op == "swap":
                            apply(ps.Xgate_Bool(*a), tuple(controls) + (b,))
                            apply(ps.Xgate_Bool(*b), tuple(controls) + (a,))
                elif node.op == "add_const":
                    bits = actual(node.operands[0])
                    value = (-node.value if inverse else node.value) % (1 << len(bits))
                    for offset in range(len(bits)):
                        if (value >> offset) & 1:
                            for i in reversed(range(offset + 1, len(bits))):
                                apply(
                                    ps.Xgate_Bool(*bits[i]), tuple(controls) + tuple(bits[offset:i])
                                )
                            apply(ps.Xgate_Bool(*bits[offset]), controls)
                else:
                    matrix = gate_matrix(node.op, node.angle, inverse)
                    for pair in actual(node.operands[0]):
                        apply(ps.Rot_Bool(*pair, [v for row in matrix for v in row]), controls)
                if not isinstance(node, NativeSite):
                    report["gate_events"] += 1
                for pair in reversed(zeros):
                    ps.Xgate_Bool(*pair)(state)
            result = {}
            for basis in state.basis_states:
                key = tuple(
                    (
                        int(basis.get(ps.System.get_id(names[r.name])).value)
                        & ((1 << r.type.width) - 1)
                    )
                    if r.type.width
                    else 0
                    for r in program.main.registers
                )
                result[key] = result.get(key, 0j) + complex(basis.amplitude)
            return RegisterState(
                program.main.registers, {k: v for k, v in result.items() if abs(v) > 1e-15}
            )
        finally:
            ps.System.clear()
