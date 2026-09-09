"已序列化 RIR 的验证、后端导出与小规模执行。"

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from pyqecclang import (
    Binding,
    Operation,
    bind,
    describe_oracle,
    dumps,
    export_originir,
    loads,
    run_originir,
    run_pysparq,
    simulate,
    unresolved,
)


def main():
    parser = argparse.ArgumentParser(prog="pyqecclang")
    parser.add_argument(
        "command",
        choices=[
            "validate",
            "inspect",
            "canonicalize",
            "emit",
            "run",
            "requirements",
            "bind",
            "compile-function",
        ],
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--memory", type=Path)
    parser.add_argument("--bindings", type=Path)
    parser.add_argument(
        "--backend", choices=["reference", "originir", "pysparq"], default="reference"
    )
    parser.add_argument("--basis", choices=["default", "toffoli-u3-cz"], default="default")
    parser.add_argument("--native-arithmetic", action="store_true")
    parser.add_argument("--native-cache", default="out/native-cache")
    parser.add_argument("--function")
    parser.add_argument("--width", type=int, default=12)
    parser.add_argument("--fraction", type=int, default=6)
    parser.add_argument("--degree", type=int, default=6)
    parser.add_argument("--constants", default="{}", help="生成期常量 JSON 对象")
    parser.add_argument("--inputs", help="输入类型 JSON 对象")
    parser.add_argument("--mir-output", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "compile-function":
            from pyqecclang.algorithms.arithmetic import FixedFormat
            from pyqecclang.infrastructure.mathfunc import Index, MathConfig, compile_function

            inputs = json.loads(args.inputs) if args.inputs else None
            if inputs is not None and not isinstance(inputs, dict):
                raise ValueError("--inputs 必须为 JSON 对象")
            if inputs is not None:
                if any(isinstance(v, dict) and set(v) != {"index"} for v in inputs.values()):
                    raise ValueError("索引类型需要 {index: 位宽}")
                inputs = {
                    k: Index(v["index"]) if isinstance(v, dict) else v for k, v in inputs.items()
                }
            constants = json.loads(args.constants)
            if not isinstance(constants, dict):
                raise ValueError("--constants 必须为 JSON 对象")
            compiled = compile_function(
                args.input.read_text(),
                entry=args.function,
                fmt=FixedFormat(args.width, args.fraction),
                config=MathConfig(degree=args.degree),
                constants=constants,
                inputs=inputs,
            )
            if args.mir_output:
                args.mir_output.write_text(compiled.math_ir.dumps())
            if args.output:
                args.output.write_text(dumps(compiled.program()))
            else:
                print(dumps(compiled.program()), end="")
            return
        program = loads(args.input.read_text())
        if args.command == "inspect":
            result = (
                json.dumps(
                    describe_oracle(Operation.from_program(program)).to_dict(),
                    ensure_ascii=False,
                    indent=2,
                )
                + "\n"
            )
        elif args.command == "validate":
            missing = unresolved(program)
            state = f"open ({len(missing)} unresolved oracles)" if missing else "closed"
            result = f"RIR {program.version}: {len(program.modules)} modules; valid; {state}\n"
        elif args.command == "requirements":
            result = (
                json.dumps(
                    [asdict(item) for item in unresolved(program)], ensure_ascii=False, indent=2
                )
                + "\n"
            )
        elif args.command == "bind":
            if args.bindings is None:
                raise ValueError("bind 需要 --bindings JSON 清单")
            manifest = json.loads(args.bindings.read_text())
            bindings = {}
            for name, item in manifest.items():
                path = args.bindings.parent / item["program"]
                implementation = loads(path.read_text())
                op = Operation(
                    implementation.main,
                    tuple(m for m in implementation.modules if m.name != implementation.entry),
                )
                bindings[name] = Binding(op, item.get("resources", {}))
            result = dumps(bind(program, bindings))
        elif args.command == "canonicalize":
            result = dumps(program)
        elif args.command == "emit":
            from pyqecclang.infrastructure.backends.basis import export_toffoli_u3_cz

            exporter = export_toffoli_u3_cz if args.basis == "toffoli-u3-cz" else export_originir
            result = exporter(program).text
        else:
            memory = json.loads(args.memory.read_text()) if args.memory else None
            if memory is not None:
                memory = {
                    name: {int(k): v for k, v in cells.items()}
                    if isinstance(cells, dict)
                    else cells
                    for name, cells in memory.items()
                }
            if args.backend == "originir":
                vector = run_originir(program, memory)
                amplitudes = {
                    str(i): [complex(v).real, complex(v).imag]
                    for i, v in enumerate(vector)
                    if abs(v) > 1e-14
                }
                result = (
                    json.dumps(
                        {"basis": "qubit_index_lsb_first", "amplitudes": amplitudes}, indent=2
                    )
                    + "\n"
                )
            else:
                runner = simulate if args.backend == "reference" else run_pysparq
                options, report = {}, {}
                if args.native_arithmetic:
                    if args.backend != "pysparq":
                        raise ValueError("--native-arithmetic 需要 --backend pysparq")
                    from pyqecclang.algorithms.arithmetic import arithmetic_native_registry

                    options = {
                        "native_registry": arithmetic_native_registry(
                            program, cache_dir=args.native_cache
                        ),
                        "report": report,
                    }
                state = runner(program, memory, **options)
                result = (
                    json.dumps(
                        {
                            "execution": report,
                            "registers": [r.name for r in state.registers],
                            "amplitudes": [
                                {"values": list(key), "real": value.real, "imag": value.imag}
                                for key, value in sorted(state.amplitudes.items())
                            ],
                        },
                        indent=2,
                    )
                    + "\n"
                )
        if args.output:
            args.output.write_text(result)
        else:
            print(result, end="")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"pyqecclang: {exc}\n")


if __name__ == "__main__":
    main()
