from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves import speed_of_sound_via_bulk_modulus_and_density as sound_speed_law

# Description
## In air at 20 degrees Celsius (bulk modulus B = 142 kPa, density rho = 1.204 kg/m**3)
## and atmospheric pressure, the speed of sound is approximately 343 m/s.

Args = namedtuple("Args", "b rho")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(142 * prefixes.kilo * units.pascal)
    rho = Quantity(1.204 * units.kilogram / units.meter**3)
    return Args(b=b, rho=rho)


def test_law(test_args: Args) -> None:
    result = sound_speed_law.calculate_speed_of_sound(test_args.b, test_args.rho)
    assert_equal(result, 343 * units.meter / units.second, tolerance=2e-3)


def test_bad_bulk_modulus(test_args: Args) -> None:
    b_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        sound_speed_law.calculate_speed_of_sound(b_bad, test_args.rho)
    with raises(TypeError):
        sound_speed_law.calculate_speed_of_sound(100, test_args.rho)


def test_bad_density(test_args: Args) -> None:
    rho_bad = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        sound_speed_law.calculate_speed_of_sound(test_args.b, rho_bad)
    with raises(TypeError):
        sound_speed_law.calculate_speed_of_sound(test_args.b, 100)
