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
    p_static = Quantity(10 * units.pascal)
    p_dynamic = Quantity(2 * units.pascal)
    p_hydrostatic = Quantity(3 * units.pascal)
    Args = namedtuple("Args", "p_static, p_dynamic, p_hydrostatic")
    return Args(
        p_static=p_static,
        p_dynamic=p_dynamic,
        p_hydrostatic=p_hydrostatic,
    )


def test_basic_law(test_args):
    result = inner_pressure.calculate_inner_pressure(
        test_args.p_static,
        test_args.p_dynamic,
        test_args.p_hydrostatic,
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_inner_pressure = convert_to(result, units.pascal).evalf(5)
    assert result_inner_pressure == approx(15, 1e-5)


def test_bad_static_pressure(test_args):
    p_bad = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(p_bad, test_args.p_dynamic, test_args.p_hydrostatic)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(100, test_args.p_dynamic, test_args.p_hydrostatic)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(test_args.p_static, p_bad, test_args.p_hydrostatic)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(test_args.p_static, 100, test_args.p_hydrostatic)
    with raises(errors.UnitsError):
        inner_pressure.calculate_inner_pressure(test_args.p_static, test_args.p_dynamic, p_bad)
    with raises(TypeError):
        inner_pressure.calculate_inner_pressure(test_args.p_static, test_args.p_dynamic, 100)
