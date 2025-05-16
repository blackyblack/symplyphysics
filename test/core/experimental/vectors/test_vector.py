from collections import defaultdict
from pytest import raises
from sympy import Symbol as SymSymbol, Function as SymFunction, S, Basic, Mul, Expr
from symplyphysics import units, dimensionless, symbols
from symplyphysics.core.dimensions import dimsys_SI
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.vectors import (
    is_vector_expr,
    into_terms,
    split_factor,
    VectorSymbol,
    VectorNorm as norm,
    VectorDot as dot,
    VectorCross as cross,
    VectorMixedProduct,
    VectorFunction,
    AppliedVectorFunction,
    vector_diff,
)
from symplyphysics.core.experimental.solvers import vector_equals


def test_is_vector_expr() -> None:
    assert is_vector_expr(0)
    assert is_vector_expr(S.Zero)
    assert not is_vector_expr(1)
    assert not is_vector_expr(S(-5))

    x = SymSymbol("x")
    y = SymSymbol("y")
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    assert not is_vector_expr(x)
    assert is_vector_expr(v)

    assert not is_vector_expr(x + y)
    assert not is_vector_expr(x + v)
    assert is_vector_expr(v + w)

    assert not is_vector_expr(x * y)
    assert is_vector_expr(v * x)
    assert is_vector_expr((v + w) * x)

    u = VectorFunction("u", nargs=1)
    assert is_vector_expr(u(x))
    assert is_vector_expr(u(v))


def as_combination(expr: Expr) -> dict[Expr, Expr]:
    combination: dict[Expr, Expr] = defaultdict(lambda: S.Zero)

    for term in into_terms(expr):
        vector, factor = split_factor(term)

        combination[vector] += factor

    return combination


def check_combination(expr: Expr, combination: dict[Expr, Expr]) -> bool:
    expr_combination = as_combination(expr)

    if set(expr_combination.keys()) != set(combination.keys()):
        return False

    return all(
        expr_equals(expr_combination[vector], factor) for vector, factor in combination.items())


def test_combination() -> None:
    assert check_combination(0, {})

    x = SymSymbol("x")
    y = SymSymbol("y")
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    assert check_combination(v, {v: 1})
    assert check_combination(x * v, {v: x})
    assert check_combination(v + w + 0, {v: 1, w: 1})
    assert check_combination(x * v + y * w, {v: x, w: y})

    u = VectorFunction("u", nargs=1)

    expr = (v + w * 2 + u(x) / x) * (x + y)
    assert check_combination(expr, {v: x + y, w: 2 * (x + y), u(x): (x + y) / x})

    with raises(ValueError):
        tuple(into_terms(1))

    with raises(ValueError):
        tuple(into_terms(x))

    with raises(ValueError):
        tuple(into_terms(x + v))

    with raises(ValueError):
        tuple(into_terms(x * v + 1))


def test_init() -> None:
    name = "F"
    dim = units.force
    latex = "\\mathbf{F}"

    force = VectorSymbol(name, dim, display_latex=latex)
    assert force.display_name == name
    assert force.display_latex == latex
    assert force.dimension == dim

    dimensionless_force = VectorSymbol(name, display_latex=latex)
    assert dimensionless_force.display_name == name
    assert dimensionless_force.display_latex == latex
    assert dimensionless_force.dimension == dimensionless


def test_equality() -> None:
    a = VectorSymbol("a")
    b = VectorSymbol("b")
    assert a != b
    assert a == a  # pylint: disable=comparison-with-itself
    assert b == b  # pylint: disable=comparison-with-itself
    assert a != 0
    assert b != 0

    c1 = VectorSymbol("c")
    c2 = VectorSymbol("c")
    assert c1 != c2


def test_vector_scaling() -> None:
    force = VectorSymbol("F", units.force)
    scale = SymSymbol("k", real=True)

    scaled_force = force * scale
    assert isinstance(scaled_force, Mul)
    vector1, scale1 = split_factor(scaled_force)
    assert vector1 == force
    assert scale1 == scale

    doubly_scaled_force = force * scale * scale
    assert isinstance(doubly_scaled_force, Mul)
    vector2, scale2 = split_factor(doubly_scaled_force)
    assert vector2 == force
    assert expr_equals(scale2, scale**2)

    scaled_zero = 0 * scale
    assert scaled_zero == 0

    scaled_by_zero = force * 0
    assert scaled_by_zero == 0

    negated_force = -force
    assert isinstance(negated_force, Mul)
    vector3, scale3 = split_factor(negated_force)
    assert vector3 == force
    assert scale3 == -1

    positive_force = +1 * force
    assert positive_force == force

    # .subs

    new_scaled_force = scaled_force.subs(scale, 2)
    assert new_scaled_force == force * 2

    acceleration = VectorSymbol("a", units.acceleration)
    scaled_acceleration = scaled_force.subs(force, acceleration)
    assert scaled_acceleration == acceleration * scale


def test_vector_norm() -> None:
    # VectorScale

    force = VectorSymbol("F", units.force)
    real_scale = SymSymbol("k", real=True)
    pos_scale = SymSymbol("m", positive=True)

    assert expr_equals(norm(force * real_scale), norm(force) * abs(real_scale))
    assert expr_equals(norm(force * pos_scale), norm(force) * pos_scale)
    assert expr_equals(norm(-force), norm(force))

    # VectorAdd

    v1 = VectorSymbol("v_1")
    v2 = VectorSymbol("v_2")
    assert expr_equals(norm(v1 + v2), norm(v1 + v2, evaluate=False))
    assert expr_equals(norm(v1 + v1), norm(v1) * 2)
    assert expr_equals(norm(0 - v2), norm(v2))

    # .subs

    assert norm(force * real_scale).subs(real_scale, 3) == norm(force) * 3

    assert norm(force).subs(force, v1) == norm(v1)


def test_vector_add() -> None:
    force_1 = VectorSymbol("F_1", units.force)
    force_2 = VectorSymbol("F_2", units.force)
    force_3 = VectorSymbol("F_3", units.force)

    assert force_1 + 0 == force_1
    assert 0 + force_1 == force_1

    sum_12 = force_1 + force_2
    assert set(sum_12.args) == {force_1, force_2}
    assert vector_equals(sum_12 + 0, sum_12)
    assert vector_equals(sum_12, sum_12)
    assert sum_12 - sum_12 == 0

    sum_123 = force_1 + force_2 + force_3
    assert set(sum_123.args) == {force_1, force_2, force_3}
    assert sum_123 + sum_123 == sum_123 * 2

    assert force_1 + force_2 - force_1 == force_2
    assert force_1 + force_2 - force_2 - force_1 == 0
    assert force_1 + force_1 == force_1 * 2
    assert force_1 - force_2 + (force_1 + force_3) * 2 == force_1 * 3 - force_2 + force_3 * 2

    mass = symbols.mass
    acceleration = VectorSymbol("a", units.acceleration)
    assert set((force_1 - acceleration * mass).args) == {force_1, acceleration * mass * -1}

    # NOTE: dimensions of sub-expressions are not yet checked
    _ = force_1 + acceleration

    # .subs

    assert sum_123.subs(force_3, force_1) == force_1 * 2 + force_2


def test_vector_dot() -> None:
    f1 = VectorSymbol("F_1", units.force)
    f2 = VectorSymbol("F_2", units.force)
    v1 = VectorSymbol("v_1", units.velocity)
    v2 = VectorSymbol("v_2", units.velocity)

    assert expr_equals(dot(0, 0), 0)

    assert expr_equals(dot(f1, f1), norm(f1)**2)
    assert expr_equals(dot(f1, f2), dot(f1, f2))
    assert expr_equals(dot(0, f1), 0)
    assert expr_equals(dot(f1, 0), 0)

    assert expr_equals(dot(f1 + f2, f1), norm(f1)**2 + dot(f1, f2))
    assert expr_equals(dot(f1 + f2, f2 + f1), norm(f1)**2 + 2 * dot(f1, f2) + norm(f2)**2)
    assert expr_equals(dot(f1 + f2, f1 - f2), norm(f1)**2 - norm(f2)**2)
    assert expr_equals(dot(f1 + f2, 0), 0)
    assert expr_equals(dot(0, f1 + f2), 0)

    assert expr_equals(
        dot(f1 + f2 * 2, -v1 + v2),
        -1 * dot(f1, v1) - 2 * dot(f2, v1) + dot(f1, v2) + 2 * dot(f2, v2),
    )


def test_vector_cross() -> None:
    f1 = VectorSymbol("F_1", units.force)
    f2 = VectorSymbol("F_2", units.force)
    v1 = VectorSymbol("v_1", units.velocity)
    v2 = VectorSymbol("v_2", units.velocity)

    assert vector_equals(cross(f1, f1), 0)
    assert not vector_equals(cross(f1, f2), 0)
    assert vector_equals(cross(0, 0), 0)
    assert vector_equals(cross(0, f1), 0)
    assert vector_equals(cross(f1, 0), 0)

    assert vector_equals(cross(f1, f2), cross(f1, f2))
    assert vector_equals(cross(f1, f2), -1 * cross(f2, f1))

    assert vector_equals(cross(f1, f1 + f2), cross(f1, f2))
    assert vector_equals(cross(f1 + f2, f1 + f2), 0)

    assert vector_equals(
        cross(f1 + f2 * 2, -v1 + v2),
        -1 * cross(f1, v1) + cross(f1, v2) - 2 * cross(f2, v1) + 2 * cross(f2, v2),
    )


def test_vector_mixed() -> None:
    a = VectorSymbol("a")
    b = VectorSymbol("b")
    c = VectorSymbol("c")

    # even permutation
    assert expr_equals(VectorMixedProduct(a, b, c), VectorMixedProduct(b, c, a))
    assert expr_equals(VectorMixedProduct(a, b, c), VectorMixedProduct(c, a, b))

    # odd permutation
    assert expr_equals(VectorMixedProduct(a, b, c), -1 * VectorMixedProduct(b, a, c))
    assert expr_equals(VectorMixedProduct(a, b, c), -1 * VectorMixedProduct(a, c, b))
    assert expr_equals(VectorMixedProduct(a, b, c), -1 * VectorMixedProduct(c, b, a))

    # repeating arguments
    assert expr_equals(VectorMixedProduct(a, a, c), 0)
    assert expr_equals(VectorMixedProduct(a, b, a), 0)
    assert expr_equals(VectorMixedProduct(a, b, b), 0)
    assert expr_equals(VectorMixedProduct(a + b, a + b, c), 0)
    assert expr_equals(VectorMixedProduct(a, b + c, b + c), 0)
    assert expr_equals(VectorMixedProduct(a + c, b, a + c), 0)

    assert expr_equals(VectorMixedProduct(0, b, c), 0)
    assert expr_equals(VectorMixedProduct(0, 0, c), 0)
    assert expr_equals(VectorMixedProduct(a, 0, c), 0)
    assert expr_equals(VectorMixedProduct(a, 0, 0), 0)
    assert expr_equals(VectorMixedProduct(0, 0, 0), 0)
    assert expr_equals(VectorMixedProduct(0, b, 0), 0)
    assert expr_equals(VectorMixedProduct(a, b, 0), 0)

    assert expr_equals(VectorMixedProduct(a, b, c), dot(a, cross(b, c)))

    assert expr_equals(VectorMixedProduct(a - b, b - c, c - a), 0)


def test_vector_function() -> None:  # pylint: disable=too-many-statements
    a = SymSymbol("a")
    v = VectorSymbol("v")

    f = VectorFunction("f", arguments=(a,), dimension=units.length)
    assert isinstance(f(a), AppliedVectorFunction)
    assert isinstance(f(v), AppliedVectorFunction)
    assert f(a) == f(a)
    assert f(a) != f(v)
    with raises(TypeError):
        f()
    with raises(TypeError):
        f(a, v)
    with raises(TypeError):
        f(a, v, v, a)
    assert f == f  # pylint: disable=comparison-with-itself
    assert f.arguments == (a,)
    assert dimsys_SI.equivalent_dims(f.dimension, units.length)
    assert f.display_name == "f"
    assert f.display_latex == "\\mathbf{f}"

    g = VectorFunction("g")
    assert isinstance(g(), AppliedVectorFunction)
    assert isinstance(g(a), AppliedVectorFunction)
    assert isinstance(g(v), AppliedVectorFunction)
    assert isinstance(g(a, a), AppliedVectorFunction)
    assert isinstance(g(v, v), AppliedVectorFunction)
    assert isinstance(g(a, v), AppliedVectorFunction)
    assert isinstance(g(v, a), AppliedVectorFunction)
    assert isinstance(g(a, a, f(a)), AppliedVectorFunction)
    assert isinstance(g(a, v, a, v, v, v), AppliedVectorFunction)
    assert g != f
    assert g() == g()
    assert g() != g(a)
    assert g(a) == g(a)
    assert g(a) != g(v)
    assert g(a, v) == g(a, v)
    assert g(a, a, a) == g(a, a, a)
    assert g(a, v, a) != g(a, v, v)
    assert g.arguments is None
    assert dimsys_SI.is_dimensionless(g.dimension)

    h = VectorFunction("h", nargs=3)
    b = SymSymbol("b")
    c = SymSymbol("c")
    assert isinstance(h(a, b, c), AppliedVectorFunction)
    assert isinstance(h(v, b, c), AppliedVectorFunction)
    assert isinstance(h(v, v, c), AppliedVectorFunction)
    assert isinstance(h(v, v, v), AppliedVectorFunction)
    with raises(TypeError):
        h()
    with raises(TypeError):
        h(v)
    with raises(TypeError):
        h(v, c)
    with raises(TypeError):
        h(a, b, c, v)

    sigma = VectorFunction(
        "sigma",
        dimension=units.area,
        display_latex="\\mathbf{\\sigma}",
        nargs=1,
    )
    assert sigma.display_name == "sigma"
    assert sigma.dimension == units.area
    assert sigma.display_latex == "\\mathbf{\\sigma}"
    assert dimsys_SI.is_dimensionless(h.dimension)

    q = VectorFunction(dimension=units.current)
    assert q.display_name.startswith("FUN")
    assert dimsys_SI.equivalent_dims(q.dimension, units.current)
    assert q.display_latex

    # Arguments must be applied functions
    with raises(ValueError):
        q(q)

    w = SymFunction("w")
    with raises(ValueError):
        q(w)

    # Note that `sympy.Function` doesn't check for unapplied `VectorFunction` instances
    _ = w(q)  # pylint: disable=not-callable


def test_vector_derivative() -> None:
    x = SymSymbol("x")
    y = SymSymbol("y")

    v = VectorSymbol("v")
    w = VectorSymbol("w")

    f = VectorFunction("f")
    g = SymFunction("g")
    h = VectorFunction("h", nargs=1)

    # 0
    assert vector_equals(vector_diff(0, x), 0)

    # VectorSymbol
    assert vector_equals(vector_diff(v, x), 0)

    # VectorScale

    assert vector_equals(
        vector_diff(v * x, x),
        v,
    )

    assert vector_equals(
        vector_diff(v * x + w * (1 - x), x),
        v - w,
    )

    # VectorFunction

    assert vector_equals(
        vector_diff(f(), x),
        0,
    )

    assert vector_equals(
        vector_diff(f(x), x),
        vector_diff(f(x), x, evaluate=False),
    )

    assert vector_equals(
        vector_diff(f(y), x),
        0,
    )

    assert vector_equals(
        vector_diff(f(x, y), x),
        vector_diff(f(x, y), x, evaluate=False),
    )

    assert vector_equals(
        vector_diff(f(x, y), x),
        vector_diff(f(x, y), x, evaluate=False),
    )

    assert vector_equals(
        vector_diff(f(x) * x, x),
        f(x) + vector_diff(f(x), x) * x,
    )

    assert vector_equals(
        vector_diff(f(x) + h(x), x),
        vector_diff(f(x), x) + vector_diff(h(x), x),
    )

    # Requires a vector analogue of `sympy.Subs`
    with raises(NotImplementedError):
        # should evaluate to `Subs(vector_diff(f(y), y), y, vector_diff(g(x), x))`
        vector_diff(f(g(x)), x)  # pylint: disable=not-callable

    with raises(NotImplementedError):
        vector_diff(f(x, x), x)

    with raises(NotImplementedError):
        # should evaluate to `dot(Jacobian(f)(h), vector_diff(h(x)))`
        vector_diff(f(h(x)), x)

    # first argument is not a vector
    with raises(ValueError):
        vector_diff(Basic(), x)

    # first argument is not a vector
    with raises(ValueError):
        vector_diff(S.One, x)

    # differentiation symbol is a constant
    with raises(ValueError):
        vector_diff(f(x), 1)

    # differentiation symbol is not an expression
    with raises(TypeError):
        vector_diff(f(x), Basic())

    # differentiation symbol is not an expression, but is sympify-able
    with raises(ValueError):
        vector_diff(f(x), "x")
