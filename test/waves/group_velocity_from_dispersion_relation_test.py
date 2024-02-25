from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.waves import group_velocity_from_dispersion_relation as group_velocity_law

# Description
## A wave packet is propagating in space. For wavenumber k = 0.1 rad/m the angular frequency
## of the wave is 1.0 rad/s. For k = 0.2 rad/m the angular frequency is 1.41 rad/s. Then the
## group velocity of the wave packet is 4.10 m/s.

Args = namedtuple("Args", "w0 w1 k0 k1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    w0 = Quantity(1.00 * units.radian / units.second)
    w1 = Quantity(1.41 * units.radian / units.second)
    k0 = Quantity(0.10 * units.radian / units.meter)
    k1 = Quantity(0.20 * units.radian / units.meter)
    return Args(w0=w0, w1=w1, k0=k0, k1=k1)


def test_law(test_args: Args) -> None:
    result = group_velocity_law.calculate_group_velocity(test_args.w0, test_args.w1, test_args.k0, test_args.k1)
    assert_equal(result, 4.10 * units.meter / units.second)


def test_bad_frequency(test_args: Args) -> None:
    wb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        group_velocity_law.calculate_group_velocity(wb, test_args.w1, test_args.k0, test_args.k1)
    with raises(TypeError):
        group_velocity_law.calculate_group_velocity(100, test_args.w1, test_args.k0, test_args.k1)
    with raises(errors.UnitsError):
        group_velocity_law.calculate_group_velocity(test_args.w0, wb, test_args.k0, test_args.k1)
    with raises(TypeError):
        group_velocity_law.calculate_group_velocity(test_args.w0, 100, test_args.k0, test_args.k1)


def test_bad_wavenumber(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        group_velocity_law.calculate_group_velocity(test_args.w0, test_args.w1, kb, test_args.k1)
    with raises(TypeError):
        group_velocity_law.calculate_group_velocity(test_args.w0, test_args.w1, 100, test_args.k1)
    with raises(errors.UnitsError):
        group_velocity_law.calculate_group_velocity(test_args.w0, test_args.w1, test_args.k0, kb)
    with raises(TypeError):
        group_velocity_law.calculate_group_velocity(test_args.w0, test_args.w1, test_args.k0, 100)
