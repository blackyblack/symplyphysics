from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
import symplyphysics.laws.hydro.inner_pressure_of_fluid as inner_pressure


@fixture(name="test_args")
def test_args_fixture():
    P_static = Quantity(10 * units.pascal)
    P_dynamic = Quantity(2 * units.pascal)
    P_hydrostatic = Quantity(3 * units.pascal)
    Args = namedtuple("Args", "P_static, P_dynamic, P_hydrostatic")
    return Args(
        P_static=P_static,
        P_dynamic=P_dynamic,
        P_hydrostatic=P_hydrostatic,
    )


def test_basic_law(test_args):
    result = inner_pressure.calculate_inner_pressure(
        test_args.P_static,
        test_args.P_dynamic,
        test_args.P_hydrostatic,
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_inner_pressure = convert_to(result, units.pascal).evalf(5)
    assert result_inner_pressure == approx(15, 1e-5)


def test_bad_static_pressure(test_args):
    P_static_bad = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(P_static_bad, test_args.P_dynamic, test_args.P_hydrostatic)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(100, test_args.P_dynamic, test_args.P_hydrostatic)


def test_bad_dynamic_pressure(test_args):
    P_dynamic_bad = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(test_args.P_static, P_dynamic_bad, test_args.P_hydrostatic)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(test_args.P_static, 100, test_args.P_hydrostatic)
    

def test_bad_hydrostatic_pressure(test_args):
    P_hydrostatic_bad = Quantity(1 * units.kilogram)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(test_args.P_static, test_args.P_dynamic, P_hydrostatic_bad)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(test_args.P_static, test_args.P_dynamic, 100)
