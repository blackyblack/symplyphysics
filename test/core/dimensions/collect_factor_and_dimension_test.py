from pytest import raises
from sympy import evaluate, Min, sin, pi, Derivative, sympify
from sympy.physics.units.systems.si import dimsys_SI
from symplyphysics import symbols, units, Quantity, assert_equal, dimensionless
from symplyphysics.core.dimensions import collect_factor_and_dimension


def test_quantity() -> None:
    expr = Quantity(4.4e4 * units.meter)
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 4.4e4)
    assert dimsys_SI.equivalent_dims(dim, units.length)

    # the factor is not always the SI value of the quantity
    expr = Quantity(1e-3 * units.kilogram)
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 1)
    assert dimsys_SI.equivalent_dims(dim, units.mass)


def test_mul() -> None:
    force_qty = Quantity(4 * units.newton)
    speed_qty = Quantity(5e-3 * units.meter / units.second)
    expr = force_qty * speed_qty
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 20)
    assert dimsys_SI.equivalent_dims(dim, units.power)


def test_add() -> None:
    first_time = Quantity(4 * units.second)
    second_time = Quantity(10 * units.second)
    third_time = Quantity(6 * units.second)
    mass = Quantity(units.kilogram)

    expr = first_time + second_time + third_time
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 20)
    assert dimsys_SI.equivalent_dims(dim, units.time)

    expr = first_time + second_time + mass
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_abs() -> None:
    with evaluate(False):
        expr = abs(Quantity(-4.4 * units.coulomb))
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 4.4)
    assert dimsys_SI.equivalent_dims(dim, units.charge)


def test_min() -> None:
    first_time = Quantity(4 * units.second)
    second_time = Quantity(10 * units.second)
    third_time = Quantity(6 * units.second)
    mass = Quantity(units.kilogram)

    with evaluate(False):
        expr = Min(first_time, second_time, third_time)
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 4)
    assert dimsys_SI.equivalent_dims(dim, units.time)

    with evaluate(False):
        expr = Min(first_time, second_time, mass)
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_function() -> None:
    dimensionless_qty = Quantity(4.4 * units.degree)
    with evaluate(False):
        expr = sin(dimensionless_qty)
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, sin(4.4 * pi / 180))
    assert dimsys_SI.equivalent_dims(dim, dimensionless)

    dimensionful_qty = Quantity(-2 * units.meter)
    with evaluate(False):
        expr = sin(dimensionful_qty)
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_derivative() -> None:
    with evaluate(False):
        expr = Derivative(Quantity(3), symbols.length)
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_rest() -> None:
    expr = sympify(4)
    factor, dim = collect_factor_and_dimension(expr)
    assert_equal(factor, 4)
    assert dimsys_SI.is_dimensionless(dim)
