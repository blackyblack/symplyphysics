from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import (
    corrected_planet_period_squared_is_proportional_to_cube_of_semimajor_axis as law,
)
from symplyphysics.quantities import solar_mass

Args = namedtuple("Args", "a ms mp")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(units.astronomical_unit)
    ms = solar_mass
    mp = solar_mass / 4
    return Args(a=a, ms=ms, mp=mp)


def test_law(test_args: Args) -> None:
    result = law.calculate_rotation_period(test_args.a, test_args.ms, test_args.mp)
    assert_equal(result, 0.89 * units.year, tolerance=5e-3)


def test_bad_length(test_args: Args) -> None:
    ab = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_rotation_period(ab, test_args.ms, test_args.mp)
    with raises(TypeError):
        law.calculate_rotation_period(100, test_args.ms, test_args.mp)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_rotation_period(test_args.a, mb, test_args.mp)
    with raises(TypeError):
        law.calculate_rotation_period(test_args.a, 100, test_args.mp)
    with raises(errors.UnitsError):
        law.calculate_rotation_period(test_args.a, test_args.ms, mb)
    with raises(TypeError):
        law.calculate_rotation_period(test_args.a, test_args.ms, 100)
