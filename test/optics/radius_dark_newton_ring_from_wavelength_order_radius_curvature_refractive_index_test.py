from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (units, Quantity, errors, assert_equal)
from symplyphysics.laws.optics import radius_dark_newton_ring_from_wavelength_order_radius_curvature_refractive_index as radius_law

# Description
## Let the wavelength be 589 nanometer, the order of ring is 6,
## radius of curvature is 4.99 meter and the refractive index is 1.
## Then the radius of the dark Newtons's ring will be 4.2 millimeter.
## https://uchi.ru/otvety/questions/diametr-shestogo-temnogo-koltsa-nyutona-okazalsya-ravnim-d6-84mm-ustanovka-dlya-polucheni?utm_referrer=https%3a%2f%2fwww.google.ru%2f

Args = namedtuple("Args", ["wavelength", "order_of_ring", "radius_curvature", "refractive_index_between_lens_plate"])

@fixture(name="test_args")
def test_args_fixture() -> Args:
    wavelength = Quantity(589 * units.nanometer)
    order_of_ring = 6
    radius_curvature = Quantity(4.99 * units.meter)
    refractive_index_between_lens_plate = 1
    
    return Args(wavelength=wavelength, order_of_ring=order_of_ring, refractive_index_between_lens_plate=refractive_index_between_lens_plate, radius_curvature=radius_curvature)


def test_basic_radius(test_args: Args) -> None:
    result = radius_law.calculate_radius(test_args.wavelength, test_args.order_of_ring, test_args.refractive_index_between_lens_plate, test_args.radius_curvature)
    assert_equal(result, 4.2 * units.millimeter)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(wavelength, test_args.order_of_ring, test_args.refractive_index_between_lens_plate, test_args.radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(100, test_args.order_of_ring, test_args.refractive_index_between_lens_plate, test_args.radius_curvature)


def test_bad_order_of_ring(test_args: Args) -> None:
    order_of_ring = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, order_of_ring, test_args.refractive_index_between_lens_plate, test_args.radius_curvature)
    with raises(ValueError):
        radius_law.calculate_radius(test_args.wavelength, -1, test_args.refractive_index_between_lens_plate, test_args.radius_curvature)


def test_bad_refractive_index_between_lens_plate(test_args: Args) -> None:
    refractive_index_between_lens_plate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_of_ring, refractive_index_between_lens_plate, test_args.radius_curvature)


def test_bad_radius_curvature(test_args: Args) -> None:
    radius_curvature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_of_ring, test_args.refractive_index_between_lens_plate, radius_curvature)
    with raises(TypeError):
        radius_law.calculate_radius(test_args.wavelength, test_args.order_of_ring, test_args.refractive_index_between_lens_plate, 100)
