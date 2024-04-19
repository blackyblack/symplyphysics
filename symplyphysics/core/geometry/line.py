from sympy import symbols, Point2D, Line2D, Expr, solve


def two_point_function(p1: Point2D, p2: Point2D, x: Expr) -> Expr:
    """Constructs a linear function of ``x`` using two points ``p1`` and ``p2``.
    
    Raises ``ValueError`` if the line equation does not depend on ``y`` or if the points are not unique.
    """

    x_, y_ = symbols("x y")
    line = Line2D(p1, p2)
    equation = line.equation(x=x_, y=y_)

    y_functions = solve(equation, y_)

    if not y_functions:
        raise ValueError("Line equation does not depend on y.")

    return y_functions[0].subs(x_, x)
