from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import isentropic_speed_of_sound as speed_law

# Description
## When gas density increased by 1 mg/L, the pressure in the gas increased by around 110 Pa.
## The speed of sound in the gas is 332 m/s.

Args = namedtuple("Args", "drho dp")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    drho = Quantity(1 * units.milligram / units.liter)
    dp = Quantity(110 * units.pascal)
    return Args(drho=drho, dp=dp)


def test_law(test_args: Args) -> None:
    result = speed_law.calculate_speed_of_sound(test_args.drho, test_args.dp)
    assert_equal(result, 332 * units.meter / units.second, tolerance=2e-3)


def test_bad_density(test_args: Args) -> None:
    rhob = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed_of_sound(rhob, test_args.dp)
    with raises(TypeError):
        speed_law.calculate_speed_of_sound(100, test_args.dp)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        speed_law.calculate_speed_of_sound(test_args.drho, pb)
    with raises(TypeError):
        speed_law.calculate_speed_of_sound(test_args.drho, 100)
