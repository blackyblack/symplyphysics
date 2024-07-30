from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    errors,
    QuantityVector,
)
from symplyphysics.laws.dynamics.vector import (
    kinetic_energy_via_angular_momentum_and_angular_velocity as law,
)

Args = namedtuple("Args", "l w")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    l_unit = units.kilogram * units.meter**2 / units.second
    l = QuantityVector([3 * l_unit, -2 * l_unit, 4 * l_unit])

    w_unit = units.radian / units.second
    w = QuantityVector([1 * w_unit, 1 * w_unit, 2 * w_unit])

    return Args(l=l, w=w)


def test_law(test_args: Args) -> None:
    result = law.calculate_kinetic_energy(test_args.l, test_args.w)
    assert_equal(result, 4.5 * units.joule)


def test_bad_angular_momentum(test_args: Args) -> None:
    lb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy(lb_vector, test_args.w)

    lb_scalar = units.kilogram * units.meter**2 / units.second
    with raises(AttributeError):
        law.calculate_kinetic_energy(lb_scalar, test_args.w)

    with raises(TypeError):
        law.calculate_kinetic_energy(100, test_args.w)
    with raises(TypeError):
        law.calculate_kinetic_energy([100], test_args.w)


def test_bad_angular_velocity(test_args: Args) -> None:
    wb_vector = QuantityVector([units.coulomb])
    with raises(errors.UnitsError):
        law.calculate_kinetic_energy(test_args.l, wb_vector)

    wb_scalar = units.radian / units.second
    with raises(AttributeError):
        law.calculate_kinetic_energy(test_args.l, wb_scalar)

    with raises(TypeError):
        law.calculate_kinetic_energy(test_args.l, 100)
    with raises(TypeError):
        law.calculate_kinetic_energy(test_args.l, [100])
