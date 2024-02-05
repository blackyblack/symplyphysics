from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.hydro import pressure_from_force_and_area as pressure


@fixture(name="test_args")
def test_args_fixture():
    force = Quantity(588 * units.newton)
    area = Quantity(0.002 * units.meter**2)
    Args = namedtuple("Args", ["area", "force"])
    return Args(
        area=area,
        force=force,
    )


def test_basic_pressure(test_args):
    result = pressure.calculate_pressure(test_args.force, test_args.area)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(5)
    assert_approx(result_pressure, 294000)


def test_bad_force(test_args):
    ag = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure.calculate_pressure(ag, test_args.area)
    with raises(errors.UnitsError):
        pressure.calculate_pressure(test_args.force, ag)
    with raises(TypeError):
        pressure.calculate_pressure(100, test_args.area)
    with raises(TypeError):
        pressure.calculate_pressure(test_args.force, 100)
