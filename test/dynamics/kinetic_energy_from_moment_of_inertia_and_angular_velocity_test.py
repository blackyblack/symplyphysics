from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import kinetic_energy_from_moment_of_inertia_and_angular_velocity as kinetic_energy_law

# Description
## If we have a stone with 2kg*m^2 of moment of inertia spinning with 3rad/s angular velocity, this stone should bear 9 joules of kinetic energy.

Args = namedtuple("Args", ["I", "w"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    I = Quantity(2 * units.kilogram * units.meter**2)
    w = Quantity(3 * units.radian / units.second)
    return Args(I=I, w=w)


def test_basic_energy(test_args: Args) -> None:
    result = kinetic_energy_law.calculate_energy(test_args.I, test_args.w)
    assert_equal(result, 9 * units.joule)


def test_bad_inertia_moment(test_args: Args) -> None:
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(Ib, test_args.w)
    with raises(TypeError):
        kinetic_energy_law.calculate_energy(100, test_args.w)


def test_bad_velocity(test_args: Args) -> None:
    wb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(test_args.I, wb)
    with raises(TypeError):
        kinetic_energy_law.calculate_energy(test_args.I, 100)
