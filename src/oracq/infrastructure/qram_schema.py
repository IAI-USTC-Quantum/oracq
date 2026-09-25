"Parsing, loading and writing of QRAM memory definition files (``*.qram.yaml``)."

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from oracq.infrastructure.execution import check_memory
from oracq.infrastructure.ir import Program, ValidationError

# Supported data types; files store source data types, and loading uniformly encodes them as unsigned words.
QRAM_DATA_TYPES = ("uint", "sint", "fixedpoint")


@dataclass(frozen=True)
class QramSegment:
    """One segment in a qram YAML file.

    Attributes:
        name: storage name, indexed by entry resource name at execution time.
        address_length: address bit width, 1..64, corresponding to the RIR QRAM address_width.
        word_length: data word bit width, 1..64, corresponding to the RIR QRAM data_width.
        type: data type, see QRAM_DATA_TYPES; data is interpreted under this type.
        data: dense source data array with the index as the address; cells above it are zero-padded when shorter than 2^address_length.
    """

    name: str
    address_length: int
    word_length: int
    type: str
    data: tuple[int, ...] | tuple[float, ...]


def _require_width(value: object, what: str) -> int:
    if type(value) is not int or not 1 <= value <= 64:
        raise ValidationError(f"{what} must be an integer in 1..64")
    return value


def _parse_words(name: str, word_length: int, data_type: str, words: object) -> tuple[int, ...] | tuple[float, ...]:
    limit = 1 << word_length
    if not isinstance(words, list):
        raise ValidationError(f"data of segment {name} must be an array")
    for value in words:
        if data_type == "uint":
            if type(value) is not int or not 0 <= value < limit:
                raise ValidationError(f"words of segment {name} must be unsigned integers within the word_length bit width")
        elif data_type == "sint":
            if type(value) is not int or not -(limit >> 1) <= value < limit >> 1:
                raise ValidationError(f"words of segment {name} must be signed integers within the word_length bit width")
        elif type(value) not in (int, float) or not 0 <= value < 1:
            raise ValidationError(f"fixedpoint values of segment {name} must be at least 0 and less than 1")
    return tuple(words)  # type: ignore[return-value]  # element types have been checked one by one per data_type.


def _encode_word(segment: QramSegment, value: int | float) -> int:
    """Encode a single source datum as an unsigned word (truncating toward zero like FixedFormat.encode)."""
    if segment.type == "uint":
        return int(value)
    if segment.type == "sint":
        return int(value) & ((1 << segment.word_length) - 1)
    return int(value * (1 << segment.word_length))


def _encode_words(segment: QramSegment) -> list[int]:
    return [_encode_word(segment, value) for value in segment.data]


def parse_qram_yaml(document: str) -> tuple[QramSegment, ...]:
    """Parse and validate one qram YAML memory definition text.

    Args:
        document: YAML text whose top level contains only the ``qram_segments`` list.

    Returns:
        tuple: segment tuple in order of appearance in the file; data keeps the source values, interpreted per type.

    Raises:
        ValidationError: the text is not valid YAML, the structure does not match, a segment name is missing or duplicated,
            a bit width is out of range, the data type is unsupported, or the source data falls outside the range
            set by type and word_length.
    """
    try:
        root = yaml.safe_load(document)
    except yaml.YAMLError as exc:
        raise ValidationError(f"QRAM memory definition is not valid YAML: {exc}") from exc
    if not isinstance(root, Mapping) or set(root) != {"qram_segments"}:
        raise ValidationError("the QRAM memory definition must contain qram_segments and nothing else")
    segments = root["qram_segments"]
    if not isinstance(segments, list):
        raise ValidationError("qram_segments must be a list of segments")
    fields = {"name", "address_length", "word_length", "type", "data"}
    result: list[QramSegment] = []
    seen: set[str] = set()
    for index, raw in enumerate(segments):
        where = f"segment {index}"
        if not isinstance(raw, Mapping) or set(raw) != fields:
            raise ValidationError(
                f"{where} must contain exactly name, address_length, word_length, type and data"
            )
        name = raw["name"]
        if not isinstance(name, str) or not name or name in seen:
            raise ValidationError("segment names must be nonempty strings unique within the file")
        seen.add(name)
        address_length = _require_width(raw["address_length"], f"address_length of segment {name}")
        word_length = _require_width(raw["word_length"], f"word_length of segment {name}")
        data_type = raw["type"]
        if data_type not in QRAM_DATA_TYPES:
            raise ValidationError(f"data type of segment {name} must be one of {QRAM_DATA_TYPES}")
        words = raw["data"]
        if isinstance(words, list) and len(words) > 1 << address_length:
            raise ValidationError(f"data length of segment {name} cannot exceed 2^address_length")
        data = _parse_words(name, word_length, data_type, words)
        result.append(QramSegment(name, address_length, word_length, data_type, data))
    return tuple(result)


def load_qram_yaml(path: str | Path) -> dict[str, list[int]]:
    """Load a qram YAML memory definition file into the memory mapping accepted by executors.

    Source data is encoded per segment type into unsigned words (sint takes the
    two's-complement bit pattern, fixedpoint multiplies by 2^word_length then
    truncates toward zero) and zero-padded to 2^address_length length; the resource
    name set and data ranges are cross-checked against program declarations by
    check_memory at execution time.

    Args:
        path: YAML file path.

    Returns:
        dict: mapping from resource names to dense word arrays, the same shape as
        the ``memory`` parameter of execution entry points such as ``simulate``.
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
    """Write the entry's QRAM binding data as qram YAML text.

    The data is first normalized and cross-checked by check_memory (sparse dicts
    densified, with zero cells and missing high cells both treated as zero), then
    segments are generated in entry resource declaration order; on the executor
    side the words are already unsigned integers within the bit width, so the
    segment type is always uint; the dense array avoids the string-key ambiguity
    of sparse dicts.

    Args:
        program: closed RIR program providing resource name and bit width declarations.
        memory: data provided by resource name; each item is a word sequence or a sparse ``address -> word`` dict,
            and ``None`` is treated as all zeros. The resource set must exactly match the entry declaration.

    Returns:
        str: qram YAML text, written to disk or hashed by the caller.
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
