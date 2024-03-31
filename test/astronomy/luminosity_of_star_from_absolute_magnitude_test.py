from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import luminosity_of_star_from_absolute_magnitude as luminosity_law

# Description
## The absolute magnitude of red dwarfs is approximately 17. Then the luminosity value will be 0.00001584.
## https://учисьучись.рф/materials/shkolnaya-programma/astronomy/vidimayaiabsolyutnayazvezdnayavelichinasveti/

Args = namedtuple("Args", ["absolute_magnitude"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    absolute_magnitude = 17

    return Args(absolute_magnitude=absolute_magnitude)


def test_basic_luminosity(test_args: Args) -> None:
    result = luminosity_law.calculate_luminosity(test_args.absolute_magnitude)
    assert_equal(result, 0.00001584)


def test_bad_absolute_magnitude() -> None:
    absolute_magnitude = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        luminosity_law.calculate_luminosity(absolute_magnitude)
