"python -m oracq.qham: PDE JSON to general QCL derivation."

import argparse
import json
from pathlib import Path

from oracq.applications.qham.examples import example_pde
from oracq.applications.qham.linearization import Block, QHAMPlan
from oracq.applications.qham.pde import PolynomialPDE
from oracq.applications.qham.report import export_derivation


def main() -> None:
    """Parse the command line, derive the QHAM linearization from the input PDE, and export the report."""
    parser = argparse.ArgumentParser(description="Automatically derive a QHAM quantum-adapted linearization from a finite polynomial PDE")
    parser.add_argument("input", nargs="?", type=Path, help="PDE 0.1 JSON; a built-in example is selected when omitted")
    parser.add_argument(
        "--example",
        choices=["burgers", "kdv", "reaction", "coupled", "vector_burgers_2d"],
        default="burgers",
    )
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--eta", type=float, default=-0.4)
    parser.add_argument("--state-width", type=int, default=2)
    parser.add_argument("--max-blocks", type=int, default=256)
    parser.add_argument("--row", help="Query a single row, e.g. physical, one, 0,1")
    parser.add_argument("-o", "--output", type=Path, default=Path("out/qham-general/derivation"))
    args = parser.parse_args()
    try:
        if not 1 <= args.state_width <= 64 or args.max_blocks < 0:
            raise ValueError("state-width must be in the range 1 to 64 and max-blocks must be non-negative")
        pde = (
            PolynomialPDE.loads(args.input.read_text()) if args.input else example_pde(args.example)
        )
        row = None
        if args.row is not None:
            row = (
                Block("physical")
                if args.row == "physical"
                else Block(
                    "tensor", () if args.row == "one" else tuple(map(int, args.row.split(",")))
                )
            )
        result = export_derivation(
            QHAMPlan(pde, args.order),
            args.output,
            state_width=args.state_width,
            eta=args.eta,
            max_blocks=args.max_blocks,
            row=row,
        )
        print(json.dumps(result, ensure_ascii=False, indent=2))
    except (ValueError, OSError) as exc:
        parser.exit(2, str(exc) + "\n")


if __name__ == "__main__":
    main()
