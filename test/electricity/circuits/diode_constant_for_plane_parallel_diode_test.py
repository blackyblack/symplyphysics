from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.circuits import diode_constant_for_plane_parallel_diode as constant_law

# Description
## The area of the electrodes is 20 [centimeter^2]. The distance between the electrodes is 1 centimeter.
## Then the diode constant will be equal to 46.7 [microampere / volt^(3 / 2)].

Args = namedtuple("Args", ["electrode_area", "distance_between_electrodes"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    electrode_area = Quantity(20 * units.centimeter**2)
    distance_between_electrodes = Quantity(1 * units.centimeter)

    return Args(electrode_area=electrode_area,
        distance_between_electrodes=distance_between_electrodes)


def test_basic_diode_constant(test_args: Args) -> None:
    result = constant_law.calculate_diode_constant(
        test_args.electrode_area, test_args.distance_between_electrodes)
    assert_equal(result, 46.7 * prefixes.micro * units.ampere / units.volt**(3 / 2))


def test_bad_electrode_area(test_args: Args) -> None:
    electrode_area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_diode_constant(electrode_area,
            test_args.distance_between_electrodes)
    with raises(TypeError):
        constant_law.calculate_diode_constant(100,
            test_args.distance_between_electrodes)


def test_bad_distance_between_electrodes(test_args: Args) -> None:
    distance_between_electrodes = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        constant_law.calculate_diode_constant(
            test_args.electrode_area, distance_between_electrodes)
    with raises(TypeError):
        constant_law.calculate_diode_constant(
            test_args.electrode_area, 100)
