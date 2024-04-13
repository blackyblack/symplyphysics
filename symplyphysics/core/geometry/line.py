from sympy import symbols, Point, Line, Expr, solve

def two_point_function(p1: Point, p2: Point, x: Expr) -> Expr:
    """Constructs a linear function of ``x`` using two points ``p1`` and ``p2``.
    
    Raises ``NotImplementedError`` if the line equation does not depend on ``y``.
    """

    xsym = x

    x, y = symbols("x y")
    line = Line(p1, p2)
    equation = line.equation(x=x, y=y)  # type: ignore[attr-defined]
    
    y_functions = solve(equation, y)

    if not y_functions:
        return NotImplemented

    return y_functions[0].subs(x, xsym)
