from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.thermodynamics import (
    fractional_volume_change_via_small_temperature_change as law,)

Args = namedtuple("Args", "a dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(3.0e-5 / units.kelvin)
    dt = Quantity(0.4 * units.kelvin)
    return Args(a=a, dt=dt)


def test_law(test_args: Args) -> None:
    result = law.calculate_fractional_volume_change(test_args.a, test_args.dt)
    assert_equal(result, 1.2e-5)


def test_bad_expansion_coefficient(test_args: Args) -> None:
    ab = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_fractional_volume_change(ab, test_args.dt)
    with raises(TypeError):
        law.calculate_fractional_volume_change(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_fractional_volume_change(test_args.a, tb)
    with raises(TypeError):
        law.calculate_fractional_volume_change(test_args.a, 100)
