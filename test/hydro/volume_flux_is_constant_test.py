from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.hydro import volume_flux_is_constant as continuity_equation

# Description
## A fluid moves steadily through a pipe of varying cross-sectional area.
## In one section of the pipe, the cross-sectional area is 1 m^2, and the speed of the flow is 6 m/s.
## In another section of the pipe, the area is 2 m^2, and the speed of the flow is 3 m/s.
## This law should confirm the validity of these figures.


@fixture(name="test_args")
def test_args_fixture():
    area_start = Quantity(1 * units.meter**2)
    v_start = Quantity(6 * units.meter / units.second)
    area_end = Quantity(2 * units.meter**2)
    Args = namedtuple("Args", ["area_start", "v_start", "area_end"])
    return Args(
        area_start=area_start,
        v_start=v_start,
        area_end=area_end,
    )


def test_basic_flow_speed(test_args):
    result = continuity_equation.calculate_flow_speed(
        test_args.area_start, test_args.v_start, test_args.area_end
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity)
    result_velocity = convert_to(result, units.meter / units.second).evalf(3)
    assert result_velocity == approx(3, 0.001)


def test_bad_area_start(test_args):
    bad_area_start = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_flow_speed(bad_area_start, test_args.v_start, test_args.area_end)
    with raises(TypeError):
        continuity_equation.calculate_flow_speed(1, test_args.v_start, test_args.area_end)


def test_bad_v_start(test_args):
    bad_v_start = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_flow_speed(test_args.area_start, bad_v_start, test_args.area_end)
    with raises(TypeError):
        continuity_equation.calculate_flow_speed(test_args.area_start, 1, test_args.area_end)


def test_bad_area_end(test_args):
    bad_area_end = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        continuity_equation.calculate_flow_speed(test_args.area_start, test_args.v_start, bad_area_end)
    with raises(TypeError):
        continuity_equation.calculate_flow_speed(test_args.area_start, test_args.v_start, 1)
