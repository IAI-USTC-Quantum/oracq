"已序列化 RIR 的验证、后端导出与小规模执行。"

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from oracq import (
    Binding,
    Operation,
    RegisterState,
    bind_with_report,
    describe_oracle,
    dumps,
    estimate_resources,
    export_originir,
    load_qram_yaml,
    loads,
    run_originir,
    run_pysparq,
    simulate,
    unresolved,
)


def main() -> None:
    """命令行入口：按子命令验证、检查、规范化、导出或执行已序列化的 RIR。

    子命令包括 validate、inspect、canonicalize、emit、run、requirements、
    bind 与 compile-function；输入与取值错误统一经 ``parser.exit`` 以退出码 2 报告。
    """
    parser = argparse.ArgumentParser(prog="oracq")
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
            "estimate",
        ],
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--format", choices=["yaml", "json"], default="yaml", help="RIR 输出文本格式")
    parser.add_argument("--memory", type=Path)
    parser.add_argument("--bindings", type=Path)
    parser.add_argument("--report", type=Path, help="绑定诊断报告的 JSON 路径")
    parser.add_argument("--allow-open", action="store_true", help="资源分析保留未实现 oracle 的调用台账")
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
            from oracq.algorithms.common.arithmetic import FixedFormat
            from oracq.infrastructure.mathfunc import Index, MathConfig, compile_function

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
                args.output.write_text(dumps(compiled.program(), format=args.format))
            else:
                print(dumps(compiled.program(), format=args.format), end="")
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
            # state 在 validate 分支为状态描述字符串，在 run 分支为执行终态。
            state: str | RegisterState = (
                f"open ({len(missing)} unresolved oracles)" if missing else "closed"
            )
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
            bindings: dict[str, Binding | Operation] = {}
            for name, item in manifest.items():
                path = args.bindings.parent / item["program"]
                implementation = loads(path.read_text())
                op = Operation(
                    implementation.main,
                    tuple(m for m in implementation.modules if m.name != implementation.entry),
                )
                bindings[name] = Binding(op, item.get("resources", {}))
            linked = bind_with_report(program, bindings)
            if args.report:
                args.report.write_text(json.dumps(linked.report.to_dict(), ensure_ascii=False, indent=2) + "\n")
            result = dumps(linked.require(), format=args.format)
        elif args.command == "estimate":
            result = json.dumps(
                estimate_resources(program, require_closed=not args.allow_open).to_dict(),
                ensure_ascii=False, indent=2,
            ) + "\n"
        elif args.command == "canonicalize":
            result = dumps(program, format=args.format)
        elif args.command == "emit":
            from oracq.infrastructure.backends.basis import export_toffoli_u3_cz

            exporter = export_toffoli_u3_cz if args.basis == "toffoli-u3-cz" else export_originir
            result = exporter(program).text
        else:
            memory = load_qram_yaml(args.memory) if args.memory else None
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
                options: dict[str, object]
                report: dict[str, int | list[str] | str]
                options, report = {}, {}
                if args.native_arithmetic:
                    if args.backend != "pysparq":
                        raise ValueError("--native-arithmetic 需要 --backend pysparq")
                    from oracq.algorithms.common.arithmetic import arithmetic_native_registry

                    options = {
                        "native_registry": arithmetic_native_registry(
                            program, cache_dir=args.native_cache
                        ),
                        "report": report,
                    }
                # runner 按 backend 动态选择，**options 为异构关键字包，无法静态验证。
                state = runner(program, memory, **options)  # type: ignore[arg-type]
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
        parser.exit(2, f"oracq: {exc}\n")


if __name__ == "__main__":
    main()
