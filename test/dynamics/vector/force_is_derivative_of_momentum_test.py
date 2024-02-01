from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    errors,
    convert_to,
    Quantity,
    SI,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import force_is_derivative_of_momentum as force_momentum_law

# Description
## A force is acting upon a body. At one moment of time, its momentum amounts to (-1, 2, 3) kg*m/s.
## After a second, it was (3, 4, 5) kg*m/s. The force acting on the body amounts to (4, 2, 2) N.


@fixture(name="test_args")
def test_args_fixture():
    p0 = QuantityVector([
        Quantity(-1.0 * units.kilogram * units.meter / units.second),
        Quantity(2.0 * units.kilogram * units.meter / units.second),
        Quantity(3.0 * units.kilogram * units.meter / units.second),
    ])
    p1 = QuantityVector([
        Quantity(3.0 * units.kilogram * units.meter / units.second),
        Quantity(4.0 * units.kilogram * units.meter / units.second),
        Quantity(5.0 * units.kilogram * units.meter / units.second),
    ])
    dt = Quantity(1.0 * units.second)
    Args = namedtuple("Args", "p0 p1 dt")
    return Args(p0=p0, p1=p1, dt=dt)


def test_basic_law(test_args):
    result_force = force_momentum_law.calculate_force(test_args.p0, test_args.p1, test_args.dt)
    assert len(result_force.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result_force.dimension, units.force)
    for result_component, correct_value in zip(result_force.components, (4, 2, 2)):
        result_value = convert_to(result_component, units.newton).evalf(3)
        assert result_value == approx(correct_value, 1e-3)


def test_bad_momenta(test_args):
    p_bad_vector = QuantityVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ])
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(p_bad_vector, test_args.p1, test_args.dt)
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(test_args.p0, p_bad_vector, test_args.dt)

    p_scalar = Quantity(1.0 * units.kilogram * units.meter / units.second)
    with raises(AttributeError):
        force_momentum_law.calculate_force(p_scalar, test_args.p1, test_args.dt)
    with raises(AttributeError):
        force_momentum_law.calculate_force(test_args.p0, p_scalar, test_args.dt)

    with raises(TypeError):
        force_momentum_law.calculate_force(100, test_args.p1, test_args.dt)
    with raises(TypeError):
        force_momentum_law.calculate_force([100], test_args.p1, test_args.dt)
    with raises(TypeError):
        force_momentum_law.calculate_force(test_args.p0, 100, test_args.dt)
    with raises(TypeError):
        force_momentum_law.calculate_force(test_args.p0, [100], test_args.dt)


def test_bad_time(test_args):
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        force_momentum_law.calculate_force(test_args.p0, test_args.p1, tb)
    with raises(TypeError):
        force_momentum_law.calculate_force(test_args.p0, test_args.p1, 100)
