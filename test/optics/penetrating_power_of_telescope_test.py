from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.optics import penetrating_power_of_telescope as power_law

# Description
## With a lens diameter of 50 millimeters, the penetrating force of the telescope will be 11.
## https://college.ru/astronomy/course/content/chapter2/tsol3.html

Args = namedtuple("Args", ["lens_diameter"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    lens_diameter = Quantity(50 * units.millimeter)

    return Args(lens_diameter=lens_diameter)


def test_basic_speed(test_args: Args) -> None:
    result = power_law.calculate_penetrating_power(test_args.lens_diameter)
    assert_equal(result, 11)


def test_bad_lens_diameter() -> None:
    lens_diameter = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_law.calculate_penetrating_power(lens_diameter)
    with raises(TypeError):
        power_law.calculate_penetrating_power(100)
