from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro import excessive_pressure_under_curved_surface_of_bubble as pressure_law

# Description
## The coefficient of surface tension of the bubble is 0.04 [newton / meter], the radius of the bubble is 7 centimeter.
## Then the excessive pressure under the curved surface of the bubble will be 2.286 pascal.
## http://ru.solverbook.com/spravochnik/molekulyarnaya-fizika-i-termodinamika/davlenie-pod-iskrivlennoj-poverxnostyu-zhidkosti/

Args = namedtuple("Args", ["surface_tension_of_the_liquid", "radius_of_bubble"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    surface_tension_of_the_liquid = Quantity(0.04 * units.newton / units.meter)
    radius_of_bubble = Quantity(7 * units.centimeter)
    return Args(surface_tension_of_the_liquid=surface_tension_of_the_liquid,
        radius_of_bubble=radius_of_bubble)


def test_basic_excessive_pressure(test_args: Args) -> None:
    result = pressure_law.calculate_excessive_pressure(test_args.surface_tension_of_the_liquid,
        test_args.radius_of_bubble)
    assert_equal(result, 2.286 * units.pascal)


def test_bad_surface_tension_of_the_liquid(test_args: Args) -> None:
    surface_tension_of_the_liquid = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_excessive_pressure(surface_tension_of_the_liquid,
            test_args.radius_of_bubble)
    with raises(TypeError):
        pressure_law.calculate_excessive_pressure(100, test_args.radius_of_bubble)


def test_bad_radius_of_bubble(test_args: Args) -> None:
    radius_of_bubble = Quantity(10 * units.coulomb)
    with raises(errors.UnitsError):
        pressure_law.calculate_excessive_pressure(test_args.surface_tension_of_the_liquid,
            radius_of_bubble)
    with raises(TypeError):
        pressure_law.calculate_excessive_pressure(test_args.surface_tension_of_the_liquid, 100)
