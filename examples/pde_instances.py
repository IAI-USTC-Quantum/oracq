"""Compose PDE grids and initial inputs, inspect RIR, and optionally execute."""

from __future__ import annotations

import argparse

import oracq


def make_problem(grid_kind: str, initial_kind: str) -> oracq.PDEProblem:
    """Use either grid representation with any of the demonstrated initial inputs."""
    grid: oracq.PDEGrid
    if grid_kind == "uniform":
        grid = oracq.UniformGrid1D(0, 2, 2)
    else:
        derivative = (("x", 2),)
        layout = oracq.UnstructuredGridConfig(2, {derivative: 0}, value_scale=3)
        data = layout.encode_matrices({derivative: [[-2, 2], [2, -2]]}, word_length=4)
        bank = oracq.register_qram(data, address_length=2, word_length=4, name="mesh")
        grid = oracq.UnstructuredGrid(bank, layout)
    initial: oracq.InitialInput
    if initial_kind == "uniform":
        initial = oracq.UniformInput(0.2)
    else:
        initial = oracq.ArrayInput([0.2, 0.1], encoding=initial_kind, angle_width=4)
    return oracq.PDEProblem(type="heat", grid=grid, initial=initial, final_time=0.01)


def main() -> None:
    """Prepare and inspect by default; native execution is selected explicitly."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", choices=("uniform", "qram"), default="uniform")
    parser.add_argument("--initial", choices=("uniform", "gates", "qram"), default="gates")
    parser.add_argument("--backend", choices=("pysparq", "originir_ext"))
    args = parser.parse_args()
    instance = oracq.qpde_solve(make_problem(args.grid, args.initial), config=oracq.QHAMConfig(order=1, taylor_degree=1)).prepare()
    circuit = instance.circuit()
    print("Input contract:", instance.prepared.report.ok)
    print("RIR modules:", len(circuit.program.modules))
    print("Resources:", circuit.resource_estimate().to_dict())
    if args.backend:
        result = instance.run_pysparq() if args.backend == "pysparq" else instance.run_originir_ext()
        print("Physical success probability:", result.success_probability)
        print("Approximate physical values:", result.physical_values())


if __name__ == "__main__":
    main()
