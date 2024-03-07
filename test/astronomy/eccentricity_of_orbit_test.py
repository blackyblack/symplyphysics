from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.astronomy import eccentricity_of_orbit as eccentricity_law

# Description
## Let the minor axis be 0.5 meter and the major axis be 1 meter. Then the eccentricity of the ellipse is 0.866.
## https://www.symbolab.com/solver/ellipse-function-eccentricity-calculator/eccentricity%201%5Ccdot%20x%5E%7B2%7D%20%2B%204%5Ccdot%20y%5E%7B2%7D%3D1?or=input

Args = namedtuple("Args", ["small_semi_axis", "large_semi_axis"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    small_semi_axis = Quantity(1 * units.meter)
    large_semi_axis = Quantity(2 * units.meter)

    return Args(small_semi_axis=small_semi_axis, large_semi_axis=large_semi_axis)


def test_basic_eccentricity(test_args: Args) -> None:
    result = eccentricity_law.calculate_eccentricity(test_args.small_semi_axis, test_args.large_semi_axis)
    assert_equal(result, 0.866)


def test_bad_semi_axis(test_args: Args) -> None:
    bad_semi_axis = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        eccentricity_law.calculate_eccentricity(bad_semi_axis, test_args.large_semi_axis)
    with raises(TypeError):
        eccentricity_law.calculate_eccentricity(100, test_args.large_semi_axis)
    with raises(errors.UnitsError):
        eccentricity_law.calculate_eccentricity(test_args.small_semi_axis, bad_semi_axis)
    with raises(TypeError):
        eccentricity_law.calculate_eccentricity(test_args.small_semi_axis, 100)
