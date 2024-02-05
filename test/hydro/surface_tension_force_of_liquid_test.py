from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, units, SI, convert_to, Quantity, errors, prefixes)
from symplyphysics.laws.hydro import surface_tension_force_of_liquid as force_law

# Description
## Let the liquid be water with a surface tension coefficient equal to 7.286e-2 newton.
## Then, with a contour length of 3 centimeter, the surface tension force will be equal to 2185.8 micronewton.
## https://www.indigomath.ru//raschety/BP2MOO.html


@fixture(name="test_args")
def test_args_fixture():
    surface_coefficient = Quantity(7.286e-2 * units.newton / units.meter)
    contour_length = Quantity(3 * units.centimeter)

    Args = namedtuple("Args", ["surface_coefficient", "contour_length"])
    return Args(surface_coefficient=surface_coefficient, contour_length=contour_length)


def test_basic_force(test_args):
    result = force_law.calculate_force(test_args.surface_coefficient, test_args.contour_length)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result = convert_to(result, prefixes.micro * units.newton).evalf(5)
    assert_approx(result, 2185.8)


def test_bad_surface_coefficient(test_args):
    surface_coefficient = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        force_law.calculate_force(surface_coefficient, test_args.contour_length)
    with raises(TypeError):
        force_law.calculate_force(100, test_args.contour_length)


def test_bad_contour_length(test_args):
    contour_length = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        force_law.calculate_force(test_args.surface_coefficient, contour_length)
    with raises(TypeError):
        force_law.calculate_force(test_args.surface_coefficient, 100)
