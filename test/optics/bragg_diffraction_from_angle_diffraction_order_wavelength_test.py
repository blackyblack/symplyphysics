from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.optics import bragg_diffraction_from_angle_diffraction_order_wavelength as distance_law

# Description
## Let the diffraction order be 1, the wavelength is 700 nanometer, and the sliding angle is
## 45 degree (pi / 4 radian). Then the distance between the crystal planes will be 494.9 nanometer.
## https://www.indigomath.ru//raschety/ZfLMDK.html

Args = namedtuple("Args", ["diffraction_order", "wavelength", "angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    diffraction_order = 1
    wavelength = Quantity(700 * units.nanometer)
    angle = pi / 4
    return Args(diffraction_order=diffraction_order, wavelength=wavelength, angle=angle)


def test_basic_distance(test_args: Args) -> None:
    result = distance_law.calculate_distance(test_args.diffraction_order, test_args.wavelength,
        test_args.angle)
    assert_equal(result, 494.9 * units.nanometer)


def test_bad_diffraction_order(test_args: Args) -> None:
    diffraction_order = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(diffraction_order, test_args.wavelength, test_args.angle)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(test_args.diffraction_order, wavelength, test_args.angle)
    with raises(TypeError):
        distance_law.calculate_distance(test_args.diffraction_order, 100, test_args.angle)


def test_bad_angle(test_args: Args) -> None:
    angle = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(test_args.diffraction_order, test_args.wavelength, angle)
