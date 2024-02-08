from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.hydro import volume_flux_is_constant as continuity_equation

# Description
## A fluid moves steadily through a pipe of varying cross-sectional area.
## In one section of the pipe, the cross-sectional area is 1 m^2, and the speed of the flow is 6 m/s.
## In another section of the pipe, the area is 2 m^2, and the speed of the flow is 3 m/s.
## This law should confirm the validity of these figures.


@fixture(name="test_args")
def test_args_fixture():
    area_before = Quantity(1 * units.meter**2)
    speed_before = Quantity(6 * units.meter / units.second)
    area_after = Quantity(2 * units.meter**2)
    Args = namedtuple("Args", ["area_before", "speed_before", "area_after"])
    return Args(
        area_before=area_before,
        speed_before=speed_before,
        area_after=area_after,
    )


def test_basic_flow_speed(test_args):
    result = continuity_equation.calculate_fluid_speed(test_args.area_before,
        test_args.speed_before, test_args.area_after)
    assert_equal(result, 3 * units.meter / units.second)


def test_bad_tube_area_before(test_args):
    bad_area_before = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_fluid_speed(bad_area_before, test_args.speed_before,
            test_args.area_after)
    with raises(TypeError):
        continuity_equation.calculate_fluid_speed(1, test_args.speed_before, test_args.area_after)


def test_bad_fluid_speed_before(test_args):
    bad_fluid_speed_before = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_fluid_speed(test_args.area_before, bad_fluid_speed_before,
            test_args.area_after)
    with raises(TypeError):
        continuity_equation.calculate_fluid_speed(test_args.area_before, 1, test_args.area_after)


def test_bad_tube_area_after(test_args):
    bad_area_after = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_fluid_speed(test_args.area_before, test_args.speed_before,
            bad_area_after)
    with raises(TypeError):
        continuity_equation.calculate_fluid_speed(test_args.area_before, test_args.speed_before, 1)
