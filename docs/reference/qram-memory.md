# QRAM memory definition (qram YAML)

**English** · <a href="../zh/reference/qram-memory.html">简体中文</a>

Memory data is not written into the program serialization text (see [RIR](rir.md)); it is bound separately as an execution input. This document defines the format of the binding file `*.qram.yaml`: one file describes a set of QRAM segments (the `qram_segments` list), cross-checked against the program by entry resource name at execution time.

```yaml
qram_segments:
  - name: values            # storage name, indexed by this name in code
    address_length: 2       # address bit width
    word_length: 3          # data word bit width
    type: uint              # data type: uint / sint / fixedpoint
    data: [1, 2, 4, 7]      # concrete data, dense array, index is the address
```

## Fields

| field | type | constraint | description |
|------|------|------|------|
| `qram_segments` | list | required, sole top-level key | a set of QRAM definitions |
| `name` | string | non-empty, unique within the file | storage name; matched exactly against entry resource names at execution |
| `address_length` | integer | 1..64 | address bit width; corresponds to `address_width` of the RIR QRAM |
| `word_length` | integer | 1..64 | data word bit width; corresponds to `data_width` of the RIR QRAM |
| `type` | string | `uint` / `sint` / `fixedpoint` | source data type of `data`; encoded as an unsigned word on load |
| `data` | array | see the table below | dense source data array; the index is the address |

Value ranges and encodings of the three data types:

| type | data element range | encoded on load as an unsigned word |
|------|--------------|---------------------|
| `uint` | unsigned integer, 0 .. 2^`word_length` − 1 | unchanged |
| `sint` | signed integer, −2^(`word_length`−1) .. 2^(`word_length`−1) − 1 | two's-complement bit pattern |
| `fixedpoint` | fixed-point number in [0, 1) | multiplied by 2^`word_length` then truncated toward zero (same convention as `FixedFormat.encode`) |

## Semantic rules

- At load time each source datum is validated locally against `word_length` (ranges in the table above); an out-of-range value is an error. The length of `data` must not exceed 2^`address_length`; when shorter, high-address cells are zero-padded (padding).
- Segment names must be unique within the file.
- At execution time the QRAM declarations inside the program are authoritative: {obj}`check_memory <oracq.infrastructure.execution.check_memory>` requires the segment resource-name set to match the entry declarations exactly, no more and no less; encoded words must fall within the program-declared bit widths, and an `address_length`/`word_length` disagreeing with the program declaration is not reported as a separate error.
- Segments produced by {obj}`dump_qram_yaml <oracq.infrastructure.qram_schema.dump_qram_yaml>` always have type `uint` (the executor-side words are already unsigned integers within the declared bit width).

## Python API

```python
from oracq import dump_qram_yaml, load_qram_yaml, simulate

memory = load_qram_yaml("memory.qram.yaml")  # -> {"values": [1, 2, 4, 7]}
state = simulate(program, memory)             # fully compatible with the existing memory argument

text = dump_qram_yaml(program, memory)        # emits text after canonical validation; the caller writes it to disk
```

{obj}`load_qram_yaml <oracq.infrastructure.qram_schema.load_qram_yaml>` returns a mapping from resource names to dense word arrays (already encoded per type and zero-padded to full length), shaped like the `memory` parameter of {obj}`simulate <oracq.infrastructure.execution.simulate>`, {obj}`run_pysparq <oracq.infrastructure.backends.pysparq.run_pysparq>`, and {obj}`run_originir <oracq.infrastructure.backends.originir.run_originir>`. `dump_qram_yaml` first cross-validates through `check_memory` and densifies the sparse dictionary, then emits segments in the order of the entry's resource declarations; when the entry declares no QRAM resource it emits an empty `qram_segments`.

## CLI

```bash
oracq run closed.rir.yaml --memory memory.qram.yaml
```

`--memory` accepts only this format. When the structure or data is invalid it reports a Chinese-language error with exit code 2 and never truncates silently. See [Command line](../manual/cli.md) for all subcommands.

## Schema

The structural constraints are in [qram-memory.schema.json](schemas/qram-memory.schema.json), including the `data` element type constraints branched on `type`. Cross-field constraints (array length versus `address_length`, word value range versus `word_length` and `type`, segment-name uniqueness) cannot be expressed in JSON Schema and are enforced by the Python loader.
