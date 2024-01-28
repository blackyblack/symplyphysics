from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import pi
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.optics import radius_newton_ring_from_wavelength_order_radius_curvature_refractive_index as radius_law

# Description
## Let the wavelength be 1500 nanometer, the order of interference is 3,
## radius of curvature is 10 centimeter and the refractive index is 5.
## Then the film thickness will be 300 nanometer.
## https://www.indigomath.ru//raschety/XdKWfj.html


@fixture(name="test_args")
def test_args_fixture():
    wavelength = Quantity(1500 * units.nanometer)
    order_interference = 3
    radius_curvature = Quantity(10 * units.centimeter)
    refractive_index = 5

    Args = namedtuple("Args", ["wavelength", "order_interference", "radius_curvature", "refractive_index"])
    return Args(wavelength=wavelength, order_interference=order_interference, refractive_index=refractive_index, radius_curvature=radius_curvature)


def test_basic_radius(test_args):
    result = radius_law.calculate_radius(test_args.wavelength, test_args.order_interference, test_args.refractive_index, test_args.radius_curvature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result = convert_to(result, units.millimeter).evalf(5)
    assert result == approx(0.3, rel=0.1)


def test_bad_wavelength(test_args):
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(wavelength, test_args.order_interference, test_args.refractive_index, test_args.radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(100, test_args.order_interference, test_args.refractive_index, test_args.radius_curvature)


def test_bad_order_interference(test_args):
    order_interference = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, order_interference, test_args.refractive_index, test_args.radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.wavelength, True, test_args.refractive_index, test_args.radius_curvature)


def test_bad_refractive_index(test_args):
    refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_interference, refractive_index, test_args.radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_interference, True, test_args.radius_curvature)


def test_bad_radius_curvature(test_args):
    radius_curvature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_interference, test_args.refractive_index, radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_interference, test_args.refractive_index, 100)
