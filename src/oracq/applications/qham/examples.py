"PDE examples for automated derivation; parameters and spatial discretization are kept separate."

from oracq.applications.qham.pde import Field, Known, PolynomialPDE


def example_pde(name: str) -> PolynomialPDE:
    """Return a built-in example PDE by name.

    The available examples are ``burgers``, ``kdv``, ``reaction``,
    ``coupled``, and ``vector_burgers_2d``; all are one-dimensional except
    ``vector_burgers_2d``, which uses the ``("x", "y")`` spatial axes.

    Args:
        name: Example name.

    Returns:
        PolynomialPDE: The PDE of the matching example; ``label`` is the
        example name.

    Raises:
        ValueError: The example name is unknown.
    """
    u, v = Field("u"), Field("v")
    if name == "burgers":
        return PolynomialPDE.from_equations(
            {"u": 0.1 * u.d("x", 2) - u * u.d("x") + Known("f")}, label=name
        )
    if name == "kdv":
        return PolynomialPDE.from_equations({"u": -u.d("x", 3) - 6 * u * u.d("x")}, label=name)
    if name == "reaction":
        return PolynomialPDE.from_equations(
            {"u": 0.1 * u.d("x", 2) + u - u**3 + Known("f")}, label=name
        )
    if name == "coupled":
        return PolynomialPDE.from_equations(
            {
                "u": -0.2 * u + 0.1 * v + 0.2 * u * v + Known("f"),
                "v": 0.1 * u - 0.3 * v - 0.1 * v * v,
            },
            label=name,
        )
    if name == "vector_burgers_2d":
        return PolynomialPDE.from_equations(
            {
                "u": 0.1 * (u.d("x", 2) + u.d("y", 2)) - u * u.d("x") - v * u.d("y"),
                "v": 0.1 * (v.d("x", 2) + v.d("y", 2)) - u * v.d("x") - v * v.d("y"),
            },
            axes=("x", "y"),
            label=name,
        )
    raise ValueError("unknown PDE example: " + name)
