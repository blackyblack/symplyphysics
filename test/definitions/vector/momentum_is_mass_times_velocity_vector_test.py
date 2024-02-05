from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import cos, pi, sin
from symplyphysics import (
    units,
    convert_to,
    Quantity,
    SI,
    errors,
    QuantityVector,
)
from symplyphysics.definitions.vector import momentum_is_mass_times_velocity_vector as momentum_def

# Description
## A particle of mass m = 0.2 kg is moving in space, its velocity vector being (-1, 0, 2) m/s.
## The vector of its linear momentum is (-0.2, 0, 0.4) kg*m/s.


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(0.2 * units.kilogram)
    v = QuantityVector([
        Quantity(-1.0 * units.meter / units.second),
        Quantity(0.0 * units.meter / units.second),
        Quantity(2.0 * units.meter / units.second),
    ])
    p = QuantityVector([
        Quantity(-0.2 * units.kilogram * units.meter / units.second),
        Quantity(0.0 * units.kilogram * units.meter / units.second),
        Quantity(0.4 * units.kilogram * units.meter / units.second),
    ])
    Args = namedtuple("Args", "m v p")
    return Args(m=m, v=v, p=p)


def test_momentum_definition(test_args):
    result = momentum_def.calculate_momentum(test_args.m, test_args.v)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.momentum)
    for result_component, correct_component in zip(result.components, test_args.p.components):
        result_value = convert_to(result_component, units.kilogram * units.meter / units.second)
        correct_value = convert_to(correct_component, units.kilogram * units.meter / units.second)
        assert result_value == approx(correct_value, 1e-3)


def test_momentum_bad_mass(test_args):
    m_bad_scalar = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(m_bad_scalar, test_args.v)
    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)


def test_momentum_bad_velocity(test_args):
    v_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, v_bad_vector)

    v_scalar = Quantity(1.0 * units.meter / units.second)
    with raises(AttributeError):
        momentum_def.calculate_momentum(test_args.m, v_scalar)

    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, [100, 100])


def test_velocity_law(test_args):
    result = momentum_def.calculate_velocity(test_args.m, test_args.p)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    for result_component, correct_component in zip(result.components, test_args.v.components):
        result_value = convert_to(result_component, units.meter / units.second)
        correct_value = convert_to(correct_component, units.meter / units.second)
        assert result_value == approx(correct_value, 1e-3)


def test_velocity_bad_mass(test_args):
    m_bad_scalar = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        momentum_def.calculate_velocity(m_bad_scalar, test_args.p)
    with raises(TypeError):
        momentum_def.calculate_velocity(100, test_args.p)


def test_velocity_bad_momentum(test_args):
    p_bad_vector = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(1.0 * units.meter),
    ])
    with raises(errors.UnitsError):
        momentum_def.calculate_velocity(test_args.m, p_bad_vector)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(AttributeError):
        momentum_def.calculate_velocity(test_args.m, p_scalar)

    with raises(TypeError):
        momentum_def.calculate_velocity(test_args.m, 100)
    with raises(TypeError):
        momentum_def.calculate_velocity(test_args.m, [100])
