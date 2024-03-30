from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.optics import radius_newton_ring_from_wavelength_order_radius_curvature_refractive_index as radius_law

# Description
## Let the wavelength be 589 nanometer, the order of interference is 6,
## radius of curvature is 4.99 meter and the refractive index is 1.
## Then the radius of the dark Newtons's ring will be 4.2 millimeter.
## https://uchi.ru/otvety/questions/diametr-shestogo-temnogo-koltsa-nyutona-okazalsya-ravnim-d6-84mm-ustanovka-dlya-polucheni?utm_referrer=https%3a%2f%2fwww.google.ru%2f


@fixture(name="test_args")
def test_args_fixture():
    wavelength = Quantity(589 * units.nanometer)
    double_order_interference = 6
    radius_curvature = Quantity(4.99 * units.meter)
    refractive_index = 1

    Args = namedtuple("Args", ["wavelength", "double_order_interference", "radius_curvature", "refractive_index"])
    return Args(wavelength=wavelength, double_order_interference=double_order_interference, refractive_index=refractive_index, radius_curvature=radius_curvature)


def test_basic_radius(test_args):
    result = radius_law.calculate_radius(test_args.wavelength, test_args.double_order_interference, test_args.refractive_index, test_args.radius_curvature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.millimeter).evalf(5)
    assert result == approx(4.2, rel=0.1)


def test_bad_wavelength(test_args):
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(wavelength, test_args.double_order_interference, test_args.refractive_index, test_args.radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(100, test_args.double_order_interference, test_args.refractive_index, test_args.radius_curvature)


def test_bad_double_order_interference(test_args):
    double_order_interference = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, double_order_interference, test_args.refractive_index, test_args.radius_curvature)
    with raises(ValueError):
        radius_law.calculate_radius(test_args.wavelength, -1, test_args.refractive_index, test_args.radius_curvature)


def test_bad_refractive_index(test_args):
    refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.double_order_interference, refractive_index, test_args.radius_curvature)


def test_bad_radius_curvature(test_args):
    radius_curvature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.double_order_interference, test_args.refractive_index, radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.wavelength, test_args.double_order_interference, test_args.refractive_index, 100)
