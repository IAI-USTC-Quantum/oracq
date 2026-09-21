"可选 PySparQ 执行适配器。RIR 的根寄存器对应原生整数寄存器。"

from __future__ import annotations

import cmath
from collections.abc import Iterator, Mapping, Sequence
from threading import Lock
from typing import cast

from pyqecclang.infrastructure.execution import (
    LocalEnter,
    LocalExit,
    RegisterState,
    check_memory,
    events,
    gate_matrix,
)
from pyqecclang.infrastructure.ir import Load, Primitive, Program, Ref, ValidationError
from pyqecclang.infrastructure.native import NativeContext, NativeRegistry, NativeSite
from pyqecclang.infrastructure.validation import locations, validate

_LOCK = Lock()


def run_pysparq(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
    native_registry: NativeRegistry | None = None,
    report: dict[str, int | list[str] | str] | None = None,
) -> RegisterState:
    """经 PySparQ 稀疏寄存器模拟器按事件执行程序。

    RIR 根寄存器映射为 PySparQ 的原生命名寄存器，门与 QRAM 查询在稀疏态上
    逐个应用。需要已安装 pysparq 的解释器，依赖缺失时抛 ValidationError。
    整个执行持有模块级互斥锁，接管前要求 PySparQ 全局寄存器表为空，并在
    结束或异常退出时清理本适配器创建的寄存器。

    Args:
        program: 待执行的程序；未提供 ``native_registry`` 时必须闭合。
        memory: 按资源名提供的 QRAM 数据，值为完整字序列或地址到字的稀疏字典；省略时全部为零，提供时必须恰好覆盖入口的全部资源。
        max_steps: 展开预算，估计步骤数超过该值即拒绝执行。
        max_states: 稀疏基态数量预算，门或原生算子应用后超出即中止。
        native_registry: 模块级原生实现注册表；提供后，开放 oracle 可由原生实现闭合。
        report: 传入的字典会被就地填充执行统计，键为 ``native_calls``、
            ``gate_events``、``native_labels`` 与 ``correctness``（后者初始化为 ``pending``）。

    Returns:
        RegisterState: 以入口各寄存器整数值元组为键的稀疏复幅度，幅度低于 1e-15 的分量被丢弃。

    Raises:
        ValidationError: 程序非法、含运行期 Store、原生实现缺失或与模块描述不匹配、
            QRAM 物化超过 2^20 项、超出展开或稀疏态预算、局部工作区未复净、
            PySparQ 全局寄存器表非空，或 pysparq 未安装。
    """
    from pyqecclang.infrastructure.linking import uses_store

    program = validate(program, require_closed=native_registry is None)
    if uses_store(program):
        raise ValidationError("PySparQ 适配器暂不支持运行期 QRAM 写（Store）")
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
            qrams: dict[str, object] = {}
            for res in program.main.resources:
                data = [0] * (1 << res.type.address_width)
                for address, value in memories[res.name].items():
                    data[address] = value
                qrams[res.name] = ps.QRAMCircuit_qutrit(
                    res.type.address_width, res.type.data_width, data
                )

            def apply(operator: object, controls: Sequence[tuple[str, int]] = ()) -> None:
                """在当前稀疏态上施加算子并检查稀疏态数量预算。"""
                if controls:
                    # pysparq 原生算子为动态后端对象，无类型存根可用。
                    operator.conditioned_by_bit(list(controls))  # type: ignore[attr-defined]
                operator(state)  # type: ignore[operator]
                if len(state.basis_states) > max_states:
                    raise ValidationError("PySparQ 执行超过稀疏态数量预算")

            def actual(ref: Ref) -> list[tuple[str, int]]:
                """把视图解析为 PySparQ 的 ``(寄存器名, 位)`` 坐标列表。"""
                return [(names[name], bit) for name, bit in locations(ref)]

            # 程序已在上方拒绝 Store，事件流中实际只会出现其余五类节点。
            for node, conditions, inverse in cast(
                "Iterator[tuple[Primitive | Load | NativeSite | LocalEnter | LocalExit, tuple[tuple[Ref, int], ...], bool]]",
                events(program, max_steps=max_steps, native_modules=native_modules),
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
                    # NativeSite 事件仅在提供了原生注册表（native_modules 非空）时产生。
                    entry = cast(NativeRegistry, native_registry).entries[node.module.name]
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
                    cast("dict[str, int]", report)["native_calls"] += 1
                    if entry.label not in cast("list[str]", report["native_labels"]):
                        cast("list[str]", report["native_labels"]).append(entry.label)
                elif isinstance(node, Load):
                    # 视图通过可逆 XOR 复制到临时整数寄存器；载入后完整反算。
                    a: str | tuple[str, int]
                    d: str
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
                    theta = cast(float, node.angle) * (-1 if inverse else 1)
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
                    value = (-cast(int, node.value) if inverse else cast(int, node.value)) % (
                        1 << len(bits)
                    )
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
                    cast("dict[str, int]", report)["gate_events"] += 1
                for pair in reversed(zeros):
                    ps.Xgate_Bool(*pair)(state)
            result: dict[tuple[int, ...], complex] = {}
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


def run_pysparq_rir(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None = None,
    *,
    max_steps: int = 1_000_000,
    max_states: int = 65536,
) -> RegisterState:
    """经 PySparQ 原生 RIR 解释器执行；与 run_pysparq 互为独立实现，用于交叉验证。

    Args:
        program: 待执行的闭合 RIR 程序；不得含运行期 QRAM 写（Store）。
        memory: 资源名到 QRAM 内容的绑定；序列按下标、映射按地址给数据字。
        max_steps: 解释器展开执行的指令步数预算上限。
        max_states: 解释器维护的基矢数预算上限。

    Returns:
        RegisterState: 入口公开寄存器空间上的末态，以各寄存器整数值元组
        为键的稀疏复振幅。
    """
    from pyqecclang.infrastructure.linking import uses_store
    from pyqecclang.infrastructure.serialization import dumps

    program = validate(program, require_closed=True)
    if uses_store(program):
        raise ValidationError("PySparQ RIR 解释器暂不支持运行期 QRAM 写（Store）")
    memories = check_memory(program, memory)
    try:
        import pysparq as ps
    except ImportError as exc:
        raise ValidationError("PySparQ RIR 执行需要已安装 pysparq 的环境") from exc
    result = ps.run_rir(dumps(program), memories, max_steps=max_steps, max_states=max_states)
    return RegisterState(program.main.registers, dict(result.amplitudes))
