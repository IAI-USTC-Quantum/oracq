"""Compose a linear ODE's operator and initial inputs, then inspect or execute."""

from __future__ import annotations

import argparse

import oracq


def make_problem(generator_kind: str, initial_kind: str) -> oracq.ODEProblem:
    """Choose operator and state access independently for the same physical system."""
    generator: oracq.GeneratorInput = oracq.MatrixInput(
        [[-0.5, 0.25], [0.25, -0.5]], encoding=generator_kind, angle_width=4,
    )
    initial: oracq.InitialInput = oracq.UniformInput(0.2) if initial_kind == "uniform" else oracq.ArrayInput(
        [0.2, 0.1], encoding=initial_kind, angle_width=4,
    )
    return oracq.ODEProblem(generator, initial, final_time=0.01,
                            components=("u", "v"), dissipative=True,
                            evidence="symmetric example generator has eigenvalues -0.25 and -0.75")


def main() -> None:
    """Prepare without native dependencies; execute only when a backend is selected."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", choices=("gates", "qram"), default="gates")
    parser.add_argument("--initial", choices=("uniform", "gates", "qram"), default="gates")
    parser.add_argument("--method", choices=("taylor", "lchs", "cbmd", "schrodingerization"), default="taylor")
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--backend", choices=("pysparq", "originir_ext"))
    args = parser.parse_args()
    instance = oracq.qode_solve(make_problem(args.generator, args.initial), method=args.method,
                               config=oracq.ODEConfig(taylor_degree=args.degree)).prepare()
    circuit = instance.circuit()
    print("Input contract:", instance.prepared.report.ok)
    print("RIR modules:", len(circuit.program.modules))
    print("Resources:", circuit.resource_estimate().to_dict())
    if args.backend:
        result = instance.run_pysparq() if args.backend == "pysparq" else instance.run_originir_ext()
        print("Physical success probability:", result.success_probability)
        print("Conditional physical amplitudes:", result.physical_amplitudes())
        if result.amplitude_scale is not None:
            print("Approximate physical values:", result.physical_values())


if __name__ == "__main__":
    main()
