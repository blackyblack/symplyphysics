from pytest import raises
from sympy import Basic, Symbol as SymSymbol, pi
from symplyphysics import units, dimensionless, symbols, assert_equal
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.experimental.vectors import VectorSymbol, ZERO, norm, VectorScale, Scale


def test_init() -> None:
    name = "F"
    dim = units.force
    latex = "\\mathbf{F}"

    force = VectorSymbol(name, dim, display_latex=latex)
    assert force.display_name == name
    assert force.display_latex == latex
    assert force.dimension == dim
    assert force.norm is None
    assert not force.is_zero

    dimensionless_force = VectorSymbol(name, display_latex=latex)
    assert dimensionless_force.display_name == name
    assert dimensionless_force.display_latex == latex
    assert dimensionless_force.dimension == dimensionless
    assert force.norm is None
    assert not force.is_zero

    force_with_norm = VectorSymbol(name, dim, norm=symbols.force, display_latex=latex)
    assert force_with_norm.norm == symbols.force
    assert not force_with_norm.is_zero

    one_newton_force = VectorSymbol(name, dim, norm=1 * units.newton, display_latex=latex)
    assert_equal(one_newton_force.norm, 1 * units.newton)
    assert not one_newton_force.is_zero

    assert VectorSymbol("a", norm=0) is ZERO

    # Correct dimension of norm
    VectorSymbol("S", units.area, norm=1 * units.meter**2)
    VectorSymbol("a", norm=3)

    # Dimension of norm is different from dimension of vector
    with raises(UnitsError):
        VectorSymbol("S", units.area, norm=1)
    with raises(UnitsError):
        VectorSymbol("A", norm=1 * units.meter**2)

    # If norm is a number, it cannot be negative
    with raises(ValueError):
        VectorSymbol("k", norm=-1)

    # This also applies to more complex expressions
    negative1 = SymSymbol("n_1", negative=True)
    negative2 = SymSymbol("n_2", negative=True)
    with raises(ValueError):
        VectorSymbol("n", norm=negative1 + negative2)

    # Norm must be an Expr
    with raises(TypeError):
        VectorSymbol("a", norm=Basic())


def test_equality() -> None:
    a = VectorSymbol("a")
    b = VectorSymbol("b")
    assert a != b
    assert a == a  # pylint: disable=comparison-with-itself
    assert b == b  # pylint: disable=comparison-with-itself
    assert a != ZERO
    assert b != ZERO

    c1 = VectorSymbol("c")
    c2 = VectorSymbol("c")
    assert c1 != c2

    assert ZERO == ZERO  # pylint: disable=comparison-with-itself


def test_zero_vector() -> None:
    assert norm(ZERO) == 0
    assert ZERO.is_zero


def test_vector_scaling() -> None:
    force = VectorSymbol("F", units.force)
    scale = SymSymbol("k", real=True)

    scaled_force = force * scale
    assert isinstance(scaled_force, VectorScale)
    assert scaled_force.vector == force
    assert scaled_force.scale == scale

    doubly_scaled_force = force * scale * scale
    assert isinstance(doubly_scaled_force, VectorScale)
    assert doubly_scaled_force.vector == force
    assert expr_equals(doubly_scaled_force.scale, scale**2)

    scaled_zero = ZERO * scale
    assert scaled_zero == ZERO

    scaled_by_zero = force * 0
    assert scaled_by_zero == ZERO

    assert ZERO * 0 == ZERO

    negated_force = -force
    assert isinstance(negated_force, VectorScale)
    assert negated_force.vector == force  # pylint: disable=no-member
    assert negated_force.scale == -1  # pylint: disable=no-member

    positive_force = +force
    assert positive_force == force

    # .subs

    new_scaled_force = scaled_force.subs(scale, 2)
    assert new_scaled_force == force * 2

    acceleration = VectorSymbol("a", units.acceleration)
    scaled_acceleration = scaled_force.subs(force, acceleration)
    assert scaled_acceleration == acceleration * scale

    # .evalf

    # TODO: for some reason this doesn't produce the intended result
    scaled_vector: VectorScale = force * pi
    with raises(AssertionError):
        assert scaled_vector.evalf(3) == force * pi.evalf(3)


def test_vector_norm() -> None:
    # VectorSymbol

    displacement = VectorSymbol("s", units.length, norm=3 * units.meter)
    assert_equal(norm(displacement), 3 * units.meter)

    # VectorScale

    force = VectorSymbol("F", units.force)
    real_scale = SymSymbol("k", real=True)
    pos_scale = SymSymbol("m", positive=True)

    assert expr_equals(norm(force * real_scale), norm(force) * abs(real_scale))
    assert expr_equals(norm(force * pos_scale), norm(force) * pos_scale)
    assert expr_equals(norm(-force), norm(force))

    unit = VectorSymbol("a", norm=1)
    assert expr_equals(norm(unit), 1)
    assert expr_equals(norm(unit * real_scale), abs(real_scale))
    assert expr_equals(norm(unit * pos_scale), pos_scale)

    # VectorAdd

    v1 = VectorSymbol("v_1")
    v2 = VectorSymbol("v_2")
    assert norm(v1 + v2).args[0] == v1 + v2
    assert norm(v1 + v1) == norm(v1) * 2
    assert norm(ZERO - v2) == norm(v2)

    # .subs

    assert norm(force * real_scale).subs(real_scale, 3) == norm(force) * 3

    # NOTE: We have to force `.doit` after substitution since `__new__` doesn't handle it
    assert norm(force).subs(force, unit).doit() == norm(unit)

    # .evalf

    # NOTE: Both sides differ slightly, (LHS - RHS) / norm(F) = -4.8e-4
    with raises(AssertionError):
        assert expr_equals(norm(force * pi).evalf(3), norm(force) * pi.evalf(3))


def test_vector_add() -> None:
    force_1 = VectorSymbol("F_1", units.force)
    force_2 = VectorSymbol("F_2", units.force)
    force_3 = VectorSymbol("F_3", units.force)

    assert ZERO + ZERO == ZERO
    assert force_1 + ZERO == force_1
    assert ZERO + force_1 == force_1

    sum_12 = force_1 + force_2
    assert set(sum_12.args) == {force_1, force_2}
    assert sum_12 + ZERO == sum_12
    assert ZERO + sum_12 == sum_12

    sum_123 = force_1 + force_2 + force_3
    assert set(sum_123.args) == {force_1, force_2, force_3}
    assert sum_123 + sum_123 == sum_123 * 2

    assert force_1 + force_2 - force_1 == force_2
    assert force_1 + force_2 - force_2 - force_1 == ZERO
    assert force_1 + force_1 == force_1 * 2
    assert force_1 - force_2 + (force_1 + force_3) * 2 == force_1 * 3 - force_2 + force_3 * 2

    mass = symbols.mass
    acceleration = VectorSymbol("a", units.acceleration)
    assert set((force_1 - Scale(mass) * acceleration).args) == {force_1, acceleration * mass * -1}

    # NOTE: dimensions of sub-expressions are not yet checked
    _ = force_1 + acceleration

    # .subs

    # NOTE: We have to force `.doit` after substitution since `__new__` doesn't handle it
    assert sum_123.subs(force_3, force_1).doit() == force_1 * 2 + force_2

    # .evalf

    # NOTE: for some reason, `pi` in the LHS simply gets evaluated with maximum pricision
    with raises(AssertionError):
        assert (force_1 * pi + force_2).evalf(3) == force_1 * pi.evalf(3) + force_2
