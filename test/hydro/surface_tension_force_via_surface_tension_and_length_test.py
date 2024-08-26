from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.hydro import surface_tension_force_via_surface_tension_and_length as force_law

# Description
## Let the liquid be water with a surface tension coefficient equal to 7.286e-2 newton.
## Then, with a contour length of 3 centimeter, the surface tension force will be equal to 2185.8 micronewton.
## https://www.indigomath.ru//raschety/BP2MOO.html

Args = namedtuple("Args", ["surface_coefficient", "contour_length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_coefficient = Quantity(7.286e-2 * units.newton / units.meter)
    contour_length = Quantity(3 * units.centimeter)
    return Args(surface_coefficient=surface_coefficient, contour_length=contour_length)


def test_basic_force(test_args: Args) -> None:
    result = force_law.calculate_force(test_args.surface_coefficient, test_args.contour_length)
    assert_equal(result, 2185.8 * prefixes.micro * units.newton)


def test_bad_surface_coefficient(test_args: Args) -> None:
    surface_coefficient = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        force_law.calculate_force(surface_coefficient, test_args.contour_length)
    with raises(TypeError):
        force_law.calculate_force(100, test_args.contour_length)


def test_bad_contour_length(test_args: Args) -> None:
    contour_length = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.surface_coefficient, contour_length)
    with raises(TypeError):
        force_law.calculate_force(test_args.surface_coefficient, 100)
