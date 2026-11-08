# Command line

**English** · <a href="../zh/manual/cli.html">简体中文</a>

The command line processes already saved RIR and also provides an entry point
for [math function compilation](math-functions.md#command-line-and-cases).
RIR input auto-detects both YAML and JSON text (the conventional extension is
`.rir.yaml`); when writing RIR text, the format is selected with
`--format yaml|json`, defaulting to `yaml`.

```bash
oracq validate program.rir.yaml
oracq inspect program.rir.yaml
oracq requirements program.rir.yaml
oracq emit program.rir.yaml -o program.originir
oracq emit program.rir.yaml --basis toffoli-u3-cz -o basis.originir
```

{obj}`validate <oracq.infrastructure.validation.validate>` checks the
structure and reports the open status. `inspect` prints the registers,
resources, and capability description of the entry. `requirements` lists the
unimplemented oracles reachable from the entry together with their call paths.

## Binding implementations

```bash
oracq bind open.rir.yaml --bindings bindings.json --report binding-report.json -o closed.rir.yaml
```

A binding manifest is a JSON file that maps slot names to the implementation's
RIR file and resource names:

```json
{
  "Function": {
    "program": "lookup.rir.yaml",
    "resources": {"table": "values"}
  }
}
```

`program` is resolved relative to the directory containing the binding
manifest, and the RIR file it points to may be YAML or JSON text. The
implementation's layout and alpha must match the slot. `--format` applies to
the RIR text written by {obj}`bind <oracq.infrastructure.linking.bind>`,
`canonicalize`, and `compile-function`; the programs parsed from the two
formats are identical field by field. For binding semantics and a runnable
example see [Tutorial: replacing oracles](../tutorials/oracle-binding.md).

## Open and closed resource analysis

```bash
oracq estimate open.rir.yaml --allow-open -o open-cost.json
oracq estimate closed.rir.yaml -o closed-cost.json
```

The open report keeps oracle call counts and unknown workspaces; known costs
must not be taken as the final total cost. For counting conventions and known
gaps see [Resource estimation](resource-estimation.md). On binding failure the
diagnostics named by `--report` are still written, and the command exits with
a non-zero status.

## Execution

```bash
oracq run closed.rir.yaml --memory memory.qram.yaml
oracq run closed.rir.yaml --memory memory.qram.yaml --backend pysparq
oracq run closed.rir.yaml --memory memory.qram.yaml --backend originir
```

The default is the [reference executor](backends.md#reference-executor). The
memory file is qram YAML (format in the reference's [QRAM memory
definitions](../reference/qram-memory.md)); it lists segments by entry
resource name, and unspecified cells are zero. A failed execution returns a
non-zero exit code; unfinished modules are never silently skipped.

## Math functions and QHAM

### Independent compiler studio

```bash
oracq gui
oracq gui --no-browser --port 8000
```

The installed package contains its own local graphical interface; no separate
demo-suite checkout or external web assets are needed. Edit a bounded pure
function, select fixed-point width, fraction bits, kernel degree and checkpoint
size, then compare the real compiler outputs against retained intermediates.
Both variants use exactly the same numeric configuration. Checkpointing can
reduce workspace while increasing gate counts; the table reports both.

The qubit peak comes from the actual strict backend layout, including its
additional scratch. Download the standard compressed RIR JSON, strict modular
gate netlist, math IR and provenance report for downstream tools. These are
logical costs before rotation synthesis, QRAM timing or fault-tolerant
scheduling; compiling alone does not certify numerical accuracy.

Each browser compilation runs in an isolated worker with a 60-second deadline.
The editor accepts loop-free sources of at most 128 AST nodes, width 2..16 and
kernel degree 1..8. The CLI below handles larger sources and static loops.
Source is interpreted by the restricted compiler frontend, never executed as
a Python script. The studio binds only to the loopback interface.

### Function and PDE command lines

```bash
oracq compile-function examples/math_functions.py --function pressure \
  --width 12 --fraction 6 --checkpoint-size 16 \
  --mir-output out/pressure.mir.json -o out/pressure.rir.yaml

python -m oracq.applications.qham --example burgers --order 2 --eta=-0.4
```

The old QHAM entry point `python -m oracq.qham` is kept for compatibility. It
and the new entry point call the same implementation; full usage is described
in [General QHAM automatic generation](qham.md).
