from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.optics import bragg_diffraction_from_angle_diffraction_order_wavelength as distance_law

# Description
## Let the diffraction order be 1, the wavelength is 700 nanometer, and the sliding angle is
## 45 degree (pi / 4 radian). Then the distance between the crystal planes will be 494 nanometer.
## https://www.indigomath.ru//raschety/ZfLMDK.html


@fixture(name="test_args")
def test_args_fixture():
    diffraction_order = 1
    wavelength = Quantity(700 * units.nanometer)
    angle = pi / 4
    Args = namedtuple("Args", ["diffraction_order", "wavelength", "angle"])
    return Args(diffraction_order=diffraction_order,
        wavelength=wavelength,
        angle=angle)


def test_basic_distance(test_args):
    result = distance_law.calculate_distance(test_args.diffraction_order, test_args.wavelength,
        test_args.angle)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_power = convert_to(result, units.nanometer).evalf(2)
    assert result_power == approx(494, 0.1)


def test_bad_diffraction_order(test_args):
    diffraction_order = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(diffraction_order, test_args.wavelength, test_args.angle)
    with raises(TypeError):
        distance_law.calculate_distance(True, test_args.wavelength, test_args.angle)


def test_bad_wavelength(test_args):
    wavelength = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(test_args.diffraction_order, wavelength, test_args.angle)
    with raises(TypeError):
        distance_law.calculate_distance(test_args.diffraction_order, 100, test_args.angle)


def test_bad_angle(test_args):
    angle = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        distance_law.calculate_distance(test_args.diffraction_order, test_args.wavelength, angle)
    with raises(AttributeError):
        distance_law.calculate_distance(test_args.diffraction_order, test_args.wavelength, True)
