"""不变量断言库：跨算法测试共享的见证原语。

四个原语均以 ``unittest.TestCase`` 实例为首参数，失败时抛 ``AssertionError``。
本文件名不带 ``test_`` 前缀，不会被 unittest discover 收集；自证测试见
``test_witness.py``。
"""

import random

from oracq import ValidationError, bind, simulate, unresolved

# simulate 会丢弃 |a| < 1e-15 的幅度，复净检查用更宽松的阈值判零。
_ZERO_AMPLITUDE = 1e-9


def _register_index(program):
    return {r.name: i for i, r in enumerate(program.main.registers)}


def _basis_indices(program, samples):
    """全部基态（小寄存器空间）或固定种子的随机抽样（大空间），确定性优先。"""
    registers = program.main.registers
    width = sum(r.type.width for r in registers)
    size = 1 << width
    if samples is not None:
        return list(samples)
    if size <= 16:
        return list(range(size))
    return sorted(random.Random(0).sample(range(size), 16))


def _initial_for_index(program, index):
    """把扁平基态下标按寄存器声明序（小端）拆成 simulate 的 initial 字典。"""
    initial = {}
    offset = 0
    for reg in program.main.registers:
        initial[reg.name] = (index >> offset) & ((1 << reg.type.width) - 1)
        offset += reg.type.width
    return initial


def _inner_product(first, second):
    return sum(
        amplitude.conjugate() * second.get(key, 0j) for key, amplitude in first.items()
    )


def assert_unitary(case, program, *, places=9, samples=None):
    """W†W = I 的抽样见证：各基态列归一（Σ|a|²=1）且两两正交。

    samples 显式给出基态下标列表；缺省时小寄存器空间（≤16 个基态）取全部，
    大空间用 ``random.Random(0)`` 固定种子抽 16 个。
    """
    indices = _basis_indices(program, samples)
    columns = [
        simulate(program, initial=_initial_for_index(program, index)).amplitudes
        for index in indices
    ]
    for index, amplitudes in zip(indices, columns, strict=True):
        norm = sum(abs(a) ** 2 for a in amplitudes.values())
        case.assertAlmostEqual(
            norm, 1.0, places=places, msg=f"基态 {index} 的输出列不归一：Σ|a|²={norm}"
        )
    for i, first in enumerate(columns):
        for j in range(i + 1, len(columns)):
            inner = _inner_product(first, columns[j])
            case.assertAlmostEqual(
                inner.real,
                0.0,
                places=places,
                msg=f"列 {indices[i]} 与列 {indices[j]} 不正交：内积 {inner}",
            )
            case.assertAlmostEqual(
                inner.imag,
                0.0,
                places=places,
                msg=f"列 {indices[i]} 与列 {indices[j]} 不正交：内积 {inner}",
            )


def assert_uncomputation(case, program, *, initial=None, work_registers=None):
    """复净见证：局部寄存器必须复净，且指定根寄存器在输出中恒为 0。

    simulate 在 LocalExit 处强制复净（未复净抛 ValidationError），此处将其
    转为 AssertionError；work_registers 列出的根寄存器须在所有非零幅度
    基态中取值为 0。返回 simulate 的 RegisterState 便于调用方继续断言。
    """
    try:
        state = simulate(program, initial=initial)
    except ValidationError as exc:
        case.fail(f"复净失败：{exc}")
    if work_registers:
        index = _register_index(program)
        for name in work_registers:
            position = index[name]
            for key, amplitude in state.amplitudes.items():
                if abs(amplitude) > _ZERO_AMPLITUDE:
                    case.assertEqual(
                        key[position],
                        0,
                        msg=f"工作寄存器 {name} 未复净：基态 {key} 幅度 {amplitude}",
                    )
    return state


def assert_bind_invariant(case, abstract_program, bindings, *, tolerance=0.0):
    """绑定不变量：同一抽象槽位的各候选实现须给出一致的可观察分布。

    abstract_program 须恰有一个未绑定槽位（经 unresolved 解析）；bindings 的
    每个 (标签, Operation) 视为该槽位的一个候选实现，逐一 bind 后 simulate，
    逐基态比较 |a|² 分布。tolerance=0 时精确对拍；QRAM 等有量化误差的绑定
    传入误差界（如 0.02，沿用现有约定）。
    """
    names = [r.name for r in unresolved(abstract_program)]
    case.assertEqual(len(names), 1, msg=f"抽象程序须恰有一个未绑定槽位：{names}")
    slot = names[0]
    reference = None
    reference_label = None
    for label, operation in bindings.items():
        bound = bind(abstract_program, {slot: operation})
        distribution = {
            key: abs(amplitude) ** 2 for key, amplitude in simulate(bound).amplitudes.items()
        }
        if reference is None:
            reference, reference_label = distribution, label
            continue
        for key in set(reference) | set(distribution):
            expected = reference.get(key, 0.0)
            actual = distribution.get(key, 0.0)
            message = (
                f"候选 {label} 与 {reference_label} 在基态 {key} 的分布不一致："
                f"{actual} != {expected}"
            )
            if tolerance:
                case.assertAlmostEqual(actual, expected, delta=tolerance, msg=message)
            else:
                case.assertEqual(actual, expected, msg=message)


def block_column(be, column):
    """BE 的 (0,0) 块第 column 列（乘以 alpha 后）。"""
    state = simulate(be.operation.program(), initial={"target": column})
    size = 1 << be.width
    return [state.amplitudes.get((row, 0), 0) * be.alpha for row in range(size)]


def assert_block_equals(case, be, matrix, *, places=9):
    """块编码的 (0,0) 块逐列与稠密矩阵对拍（乘回 alpha）。"""
    for column in range(len(matrix)):
        actual = block_column(be, column)
        for row in range(len(matrix)):
            case.assertAlmostEqual(actual[row], matrix[row][column], places=places)
