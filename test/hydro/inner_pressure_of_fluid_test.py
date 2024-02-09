from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
import symplyphysics.laws.hydro.inner_pressure_of_fluid as inner_pressure

Args = namedtuple("Args", "p_static, p_dynamic, p_hydrostatic")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p_static = Quantity(10 * units.pascal)
    p_dynamic = Quantity(2 * units.pascal)
    p_hydrostatic = Quantity(3 * units.pascal)
    return Args(
        p_static=p_static,
        p_dynamic=p_dynamic,
        p_hydrostatic=p_hydrostatic,
    )


def test_basic_law(test_args: Args) -> None:
    result = inner_pressure.calculate_inner_pressure(
        test_args.p_static,
        test_args.p_dynamic,
        test_args.p_hydrostatic,
    )
    assert_equal(result, 15 * units.pascal)


def test_bad_static_pressure(test_args: Args) -> None:
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
