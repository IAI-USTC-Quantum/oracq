"QRAM 内存定义文件（``*.qram.yaml``）的解析、加载与写出。"

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from oracq.infrastructure.execution import check_memory
from oracq.infrastructure.ir import Program, ValidationError

# 支持的数据类型；文件中写源数据类型，加载时统一编码为无符号字。
QRAM_DATA_TYPES = ("uint", "sint", "fixedpoint")


@dataclass(frozen=True)
class QramSegment:
    """qram YAML 中的一个段。

    Attributes:
        name: 存储名，执行时按入口资源名索引。
        address_length: 地址位宽，1..64，对应 RIR QRAM 的 address_width。
        word_length: 数据字位宽，1..64，对应 RIR QRAM 的 data_width。
        type: 数据类型，见 QRAM_DATA_TYPES；data 按该类型解释。
        data: 稠密源数据数组，下标即地址；短于 2^address_length 时高位单元按零补齐。
    """

    name: str
    address_length: int
    word_length: int
    type: str
    data: tuple[int, ...] | tuple[float, ...]


def _require_width(value: object, what: str) -> int:
    if type(value) is not int or not 1 <= value <= 64:
        raise ValidationError(f"{what} 必须是 1..64 的整数")
    return value


def _parse_words(name: str, word_length: int, data_type: str, words: object) -> tuple[int, ...] | tuple[float, ...]:
    limit = 1 << word_length
    if not isinstance(words, list):
        raise ValidationError(f"段 {name} 的 data 必须是数组")
    for value in words:
        if data_type == "uint":
            if type(value) is not int or not 0 <= value < limit:
                raise ValidationError(f"段 {name} 的字必须是 word_length 位宽内的无符号整数")
        elif data_type == "sint":
            if type(value) is not int or not -(limit >> 1) <= value < limit >> 1:
                raise ValidationError(f"段 {name} 的字必须是 word_length 位宽内的有符号整数")
        elif type(value) not in (int, float) or not 0 <= value < 1:
            raise ValidationError(f"段 {name} 的 fixedpoint 值必须在 [0, 1) 内")
    return tuple(words)  # type: ignore[return-value]  # 元素类型已按 data_type 逐一校验。


def _encode_word(segment: QramSegment, value: int | float) -> int:
    """把单个源数据编码为无符号字（与 FixedFormat.encode 同为向零截断）。"""
    if segment.type == "uint":
        return int(value)
    if segment.type == "sint":
        return int(value) & ((1 << segment.word_length) - 1)
    return int(value * (1 << segment.word_length))


def _encode_words(segment: QramSegment) -> list[int]:
    return [_encode_word(segment, value) for value in segment.data]


def parse_qram_yaml(document: str) -> tuple[QramSegment, ...]:
    """解析并校验一份 qram YAML 内存定义文本。

    Args:
        document: YAML 文本，顶层仅含 ``qram_segments`` 列表。

    Returns:
        tuple: 段元组，按文件内出现顺序；data 保留源数据（按 type 解释）。

    Raises:
        ValidationError: 文本不是合法 YAML、结构不符、段名缺失或重复、
            位宽越界、数据类型不受支持，或源数据越出 type 与 word_length
            规定的范围。
    """
    try:
        root = yaml.safe_load(document)
    except yaml.YAMLError as exc:
        raise ValidationError(f"QRAM 内存定义不是合法 YAML：{exc}") from exc
    if not isinstance(root, Mapping) or set(root) != {"qram_segments"}:
        raise ValidationError("QRAM 内存定义必须且只能包含 qram_segments")
    segments = root["qram_segments"]
    if not isinstance(segments, list):
        raise ValidationError("qram_segments 必须是段列表")
    fields = {"name", "address_length", "word_length", "type", "data"}
    result: list[QramSegment] = []
    seen: set[str] = set()
    for index, raw in enumerate(segments):
        where = f"第 {index} 个段"
        if not isinstance(raw, Mapping) or set(raw) != fields:
            raise ValidationError(
                f"{where} 必须且只能包含 name、address_length、word_length、type 与 data"
            )
        name = raw["name"]
        if not isinstance(name, str) or not name or name in seen:
            raise ValidationError("段名必须是非空字符串且在文件内唯一")
        seen.add(name)
        address_length = _require_width(raw["address_length"], f"段 {name} 的 address_length")
        word_length = _require_width(raw["word_length"], f"段 {name} 的 word_length")
        data_type = raw["type"]
        if data_type not in QRAM_DATA_TYPES:
            raise ValidationError(f"段 {name} 的数据类型必须是 {QRAM_DATA_TYPES} 之一")
        words = raw["data"]
        if isinstance(words, list) and len(words) > 1 << address_length:
            raise ValidationError(f"段 {name} 的 data 长度不能超过 2^address_length")
        data = _parse_words(name, word_length, data_type, words)
        result.append(QramSegment(name, address_length, word_length, data_type, data))
    return tuple(result)


def load_qram_yaml(path: str | Path) -> dict[str, list[int]]:
    """加载 qram YAML 内存定义文件为执行器接受的内存映射。

    源数据按段类型编码为无符号字（sint 取补码位模式，fixedpoint 乘以
    2^word_length 后向零截断），并补零到 2^address_length 长；资源名集合
    与数据范围在执行时由 check_memory 与程序声明交叉校验。

    Args:
        path: YAML 文件路径。

    Returns:
        dict: 资源名到稠密字数组的映射，与 ``simulate`` 等执行入口的
        ``memory`` 参数同形。
    """
    result = {}
    for segment in parse_qram_yaml(Path(path).read_text(encoding="utf-8")):
        words = _encode_words(segment)
        words.extend([0] * ((1 << segment.address_length) - len(words)))
        result[segment.name] = words
    return result


def dump_qram_yaml(
    program: Program,
    memory: Mapping[str, Sequence[int] | Mapping[int, int]] | None,
) -> str:
    """把入口的 QRAM 绑定数据写为 qram YAML 文本。

    先经 check_memory 规整与交叉校验（稀疏字典稠密化，零单元与缺省高位单元
    均按零），再按入口资源声明顺序生成段；执行器侧的字已是位宽内无符号整数，
    段类型恒为 uint；稠密数组避免稀疏字典的字符串键歧义。

    Args:
        program: 封闭 RIR 程序，提供资源名与位宽声明。
        memory: 按资源名提供的数据；每项为字序列或 ``地址 -> 字`` 的稀疏字典，
            ``None`` 视为全零。资源集合必须与入口声明完全一致。

    Returns:
        str: qram YAML 文本，由调用方写盘或哈希。
    """
    cells = check_memory(program, memory)
    segments = []
    for resource in program.main.resources:
        bank = cells[resource.name]
        words = [bank.get(address, 0) for address in range(1 << resource.type.address_width)]
        segments.append(
            {
                "name": resource.name,
                "address_length": resource.type.address_width,
                "word_length": resource.type.data_width,
                "type": "uint",
                "data": words,
            }
        )
    return yaml.safe_dump(
        {"qram_segments": segments}, sort_keys=False, allow_unicode=True, default_flow_style=False
    )
