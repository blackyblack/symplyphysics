from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal)

from symplyphysics.laws.hydro import jurin_equation

## Source of numbers: https://scask.ru/e_book_phis.php?id=62

Args = namedtuple("Args", [
    "surface_tension_coefficient",
    "angle",
    "density_of_liquid",
    "radius",
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_tension_coefficient = Quantity(0.073 * units.newton / units.meter)
    angle = 0
    density_of_liquid = Quantity(1000 * units.kilogram / units.meter**3)
    radius = Quantity(5e-7 * units.meter)
    return Args(
        surface_tension_coefficient=surface_tension_coefficient,
        angle=angle,
        density_of_liquid=density_of_liquid,
        radius=radius,
    )


def test_basic_law(test_args: Args) -> None:
    result = jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle,
        test_args.density_of_liquid, test_args.radius)
    assert_equal(result, 29.77 * units.meter)


def test_bad_surface_tension_coefficient(test_args: Args) -> None:
    bad_surface_tension_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        jurin_equation.calculate_height(bad_surface_tension_coefficient, test_args.angle,
            test_args.density_of_liquid, test_args.radius)
    with raises(TypeError):
        jurin_equation.calculate_height(100, test_args.angle, test_args.density_of_liquid,
            test_args.radius)


def test_bad_angle(test_args: Args) -> None:
    bad_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, bad_angle,
            test_args.density_of_liquid, test_args.radius)


def test_bad_density_of_liquid(test_args: Args) -> None:
    bad_density_of_liquid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle,
            bad_density_of_liquid, test_args.radius)
    with raises(TypeError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle, 100,
            test_args.radius)


def test_bad_radius(test_args: Args) -> None:
    bad_radius = Quantity(1 * units.coulomb)
    bad_radius_2 = Quantity(40 * units.meter)
    with raises(errors.UnitsError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle,
            test_args.density_of_liquid, bad_radius)
    with raises(TypeError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle,
            test_args.density_of_liquid, 100)
    with raises(ValueError):
        jurin_equation.calculate_height(test_args.surface_tension_coefficient, test_args.angle,
            test_args.density_of_liquid, bad_radius_2)
