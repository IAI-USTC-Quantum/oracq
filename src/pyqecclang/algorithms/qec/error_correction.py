"""三位重复码的相干编码与单错误恢复电路。"""

from __future__ import annotations

from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


def repetition_encode(*, error: str = "bit") -> Operation:
    """将一位逻辑态编码到三位重复码。

    Args:
        error: ``bit`` 对应单 X 错误，``phase`` 对应单 Z 错误。

    Returns:
        Operation: target 为一位逻辑输入，syndrome 为两位且输入必须为零。

    编码后的三个物理位按 target、syndrome[0]、syndrome[1] 排列。"""
    if error not in {"bit", "phase"}:
        raise ValidationError("重复码 error 只能是 bit 或 phase")
    b = Builder(
        "repetition_encode_" + error,
        {"target": Bits(1), "syndrome": Bits(2)},
        attributes={"algorithm": "repetition_encode", "error_kind": error},
    )
    b.xor(b["target"], b["syndrome"][0])
    b.xor(b["target"], b["syndrome"][1])
    if error == "phase":
        b.h(fuse(b["target"], b["syndrome"]))
    return b.finish()


def repetition_recover(*, error: str = "bit") -> Operation:
    """相干恢复三位重复码中的单个指定类型错误。

    Args:
        error: 与编码器一致的 ``bit`` 或 ``phase``。

    Returns:
        Operation: target 恢复逻辑态，错误信息保留在 syndrome。

    不包含测量或重置，不能把非零 syndrome 当作已经复净的工作区。"""
    if error not in {"bit", "phase"}:
        raise ValidationError("重复码 error 只能是 bit 或 phase")
    b = Builder(
        "repetition_recover_" + error,
        {"target": Bits(1), "syndrome": Bits(2)},
        attributes={
            "algorithm": "repetition_recover",
            "error_kind": error,
            "syndrome_policy": "retained; host reset required before reuse",
        },
    )
    if error == "phase":
        b.h(fuse(b["target"], b["syndrome"]))
    b.xor(b["target"], b["syndrome"][1])
    b.xor(b["target"], b["syndrome"][0])
    with b.control(b["syndrome"], 3):
        b.x(b["target"])
    return b.finish()
