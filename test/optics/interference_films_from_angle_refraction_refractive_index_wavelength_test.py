from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.optics import interference_films_from_angle_refraction_refractive_index_wavelength as thickness_law

# Description
## Let the wavelength be 1500 nanometer, the order of interference is 1, the angle of refraction is
## 60 degree (pi / 3 radian), and the refractive index is 5. Then the film thickness will be 300 nanometer.
## https://www.indigomath.ru//raschety/uHv0Ia.html


@fixture(name="test_args")
def test_args_fixture():
    wavelength = Quantity(1500 * units.nanometer)
    double_order_interference = 1
    angle_refraction = pi / 3
    refractive_index = 5

    Args = namedtuple("Args", ["wavelength", "double_order_interference", "angle_refraction", "refractive_index"])
    return Args(wavelength=wavelength, double_order_interference=double_order_interference, refractive_index=refractive_index, angle_refraction=angle_refraction)


def test_basic_thickness(test_args):
    result = thickness_law.calculate_film_thickness(test_args.wavelength, test_args.double_order_interference, test_args.refractive_index, test_args.angle_refraction)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.nanometer).evalf(5)
    assert result == approx(300, rel=0.1)


def test_bad_wavelength(test_args):
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(wavelength, test_args.double_order_interference, test_args.refractive_index, test_args.angle_refraction)
    with raises(TypeError):
        thickness_law.calculate_film_thickness(100, test_args.double_order_interference, test_args.refractive_index, test_args.angle_refraction)


def test_bad_double_order_interference(test_args):
    double_order_interference = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, double_order_interference, test_args.refractive_index, test_args.angle_refraction)
    with raises(ValueError):
        thickness_law.calculate_film_thickness(test_args.wavelength, -1, test_args.refractive_index, test_args.angle_refraction)


def test_bad_refractive_index(test_args):
    refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, test_args.double_order_interference, refractive_index, test_args.angle_refraction)


def test_bad_angle_refraction(test_args):
    angle_refraction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, test_args.double_order_interference, test_args.refractive_index, angle_refraction)
