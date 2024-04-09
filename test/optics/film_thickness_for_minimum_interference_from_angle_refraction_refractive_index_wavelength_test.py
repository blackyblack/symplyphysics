from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (units, Quantity, errors, assert_equal)
from symplyphysics.laws.optics import film_thickness_for_minimum_interference_from_angle_refraction_refractive_index_wavelength as thickness_law

# Description
## Let the wavelength be 400 nanometer, the order of interference is 4, the angle of refraction is
## 22 degree (0.12 * pi radian), and the refractive index is 1.33. Then the film thickness will be 647 nanometer.
## https://exir.ru/5/resh/5_80.htm

Args = namedtuple("Args", ["wavelength", "order_interference", "angle_refraction", "refractive_index"])

@fixture(name="test_args")
def test_args_fixture() -> Args:
    wavelength = Quantity(400 * units.nanometer)
    order_interference = 4
    angle_refraction = 0.12 * pi
    refractive_index = 1.33

    return Args(wavelength=wavelength, order_interference=order_interference, refractive_index=refractive_index, angle_refraction=angle_refraction)


def test_basic_thickness(test_args: Args) -> None:
    result = thickness_law.calculate_film_thickness(test_args.wavelength, test_args.order_interference, test_args.refractive_index, test_args.angle_refraction)
    assert_equal(result, 647 * units.nanometer)


def test_bad_wavelength(test_args: Args) -> None:
    wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(wavelength, test_args.order_interference, test_args.refractive_index, test_args.angle_refraction)
    with raises(TypeError):
        thickness_law.calculate_film_thickness(100, test_args.order_interference, test_args.refractive_index, test_args.angle_refraction)


def test_bad_order_interference(test_args: Args) -> None:
    order_interference = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, order_interference, test_args.refractive_index, test_args.angle_refraction)
    with raises(ValueError):
        thickness_law.calculate_film_thickness(test_args.wavelength, -1, test_args.refractive_index, test_args.angle_refraction)


def test_bad_refractive_index(test_args: Args) -> None:
    refractive_index = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, test_args.order_interference, refractive_index, test_args.angle_refraction)


def test_bad_angle_refraction(test_args: Args) -> None:
    angle_refraction = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        thickness_law.calculate_film_thickness(test_args.wavelength, test_args.order_interference, test_args.refractive_index, angle_refraction)
