from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    wavelength_of_standing_wave_in_string_with_fixed_ends as standing_wave_law,
)

# Description
## A harmonic standing wave is occurring in a string of length l = 30 cm. The wavelength
## of the third harmonic is 20 cm.

Args = namedtuple("Args", "n l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 3
    l = Quantity(30.0 * units.centimeter)
    return Args(n=n, l=l)


def test_law(test_args: Args) -> None:
    result = standing_wave_law.calculate_wavelength(test_args.n, test_args.l)
    assert_equal(result, 20.0 * units.centimeter)


def test_bad_integer_factor(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_wavelength(test_args.n, nb)


def test_bad_string_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        standing_wave_law.calculate_wavelength(test_args.n, lb)
    with raises(TypeError):
        standing_wave_law.calculate_wavelength(test_args.n, 100)
