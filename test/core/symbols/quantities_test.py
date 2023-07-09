from pytest import raises
from sympy import Derivative, cos, pi
from symplyphysics import (units, Quantity, SI, dimensionless)
from symplyphysics.core.symbols.quantities import collect_factor_and_dimension

# Test collect_factor_and_dimension


def test_basic_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    expr = a
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 10
    assert SI.get_dimension_system().equivalent_dims(d, units.speed)


def test_prefixed_collect_factor():
    a = Quantity(10 * units.kilometer)
    expr = a
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 10000
    assert SI.get_dimension_system().equivalent_dims(d, units.length)


def test_dimensionless_collect_factor():
    a = Quantity(10)
    expr = a
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 10
    assert SI.get_dimension_system().equivalent_dims(d, dimensionless)


def test_sum_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.meter / units.second)
    expr = a + b
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 15
    assert SI.get_dimension_system().equivalent_dims(d, units.speed)

    expr = b + a
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 15
    assert SI.get_dimension_system().equivalent_dims(d, units.speed)


def test_invalid_sum_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.meter)
    expr = a + b
    with raises(ValueError):
        collect_factor_and_dimension(expr)
    expr = b + a
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_mul_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(5 * units.second)
    expr = a * b
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 50
    assert SI.get_dimension_system().equivalent_dims(d, units.length)

    expr = b * a
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 50
    assert SI.get_dimension_system().equivalent_dims(d, units.length)


def test_pow_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(2)
    expr = a**b
    (f, d) = collect_factor_and_dimension(expr)
    assert f == 100
    assert SI.get_dimension_system().equivalent_dims(d, units.speed**2)


def test_invalid_pow_collect_factor():
    a = Quantity(10 * units.meter / units.second)
    b = Quantity(2)
    expr = b**a
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_dimension_collect_factor():
    (f, d) = collect_factor_and_dimension(units.length)
    assert f == 1
    assert SI.get_dimension_system().equivalent_dims(d, units.length)


def test_invalid_derivative_collect_factor():
    a = Quantity(10 * units.second)
    expr = Derivative(a**2, a)
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_function_collect_factor():
    a = Quantity(pi)
    expr = cos(a)
    # Make sure it is unevaluated
    assert expr != -1
    (f, d) = collect_factor_and_dimension(expr)
    assert f == -1
    assert SI.get_dimension_system().equivalent_dims(d, dimensionless)


def test_invalid_function_collect_factor():
    a = Quantity(pi * units.second)
    expr = cos(a)
    # Make sure it is unevaluated
    assert expr != -1
    with raises(ValueError):
        collect_factor_and_dimension(expr)


def test_dimensionless_expr_function_collect_factor():
    a = Quantity(pi * units.second)
    b = Quantity(1 * units.second)
    expr = cos(a / b)
    # Make sure it is unevaluated
    assert expr != -1
    (f, d) = collect_factor_and_dimension(expr)
    assert f == -1
    assert SI.get_dimension_system().equivalent_dims(d, dimensionless)
