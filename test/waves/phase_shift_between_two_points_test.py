from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.waves import phase_shift_between_two_points as phase_difference_law

# Description
## The distance to the first point is 10 meter, the distance to the second point is 16 meter, and the wavelength is 12 meter.
## Then the phase difference of the oscillations of the two points will be equal to pi.
## https://www.bog5.in.ua/problems/sav/koleban/probl%20kol24.html

Args = namedtuple("Args", ["distance_first_point", "distance_second_point", "wavelength"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    distance_first_point = Quantity(10 * units.meter)
    distance_second_point = Quantity(16 * units.meter)
    wavelength = Quantity(12 * units.meter)

    return Args(distance_first_point=distance_first_point,
        distance_second_point=distance_second_point,
        wavelength=wavelength)


def test_basic_phase_difference(test_args: Args) -> None:
    result = phase_difference_law.calculate_phase_difference(test_args.distance_first_point,
        test_args.distance_second_point, test_args.wavelength)
    assert_equal(result, pi)


def test_bad_distances_and_wavelength(test_args: Args) -> None:
    bad_distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        phase_difference_law.calculate_phase_difference(bad_distance,
            test_args.distance_second_point, test_args.wavelength)
    with raises(TypeError):
        phase_difference_law.calculate_phase_difference(100, test_args.distance_second_point,
            test_args.wavelength)
    with raises(errors.UnitsError):
        phase_difference_law.calculate_phase_difference(test_args.distance_first_point,
            bad_distance, test_args.wavelength)
    with raises(TypeError):
        phase_difference_law.calculate_phase_difference(test_args.distance_first_point, 100,
            test_args.wavelength)
    with raises(errors.UnitsError):
        phase_difference_law.calculate_phase_difference(test_args.distance_first_point,
            test_args.distance_second_point, bad_distance)
    with raises(TypeError):
        phase_difference_law.calculate_phase_difference(test_args.distance_first_point,
            test_args.distance_second_point, 100)
