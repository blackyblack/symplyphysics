from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import work_done_by_gas_at_constant_pressure

Args = namedtuple("Args", ["p", "v_1", "v_2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(2.52 * units.pascal)
    v_1 = Quantity(1 * units.meter**3)
    v_2 = Quantity(2 * units.meter**3)
    return Args(p=p, v_1=v_1, v_2=v_2)


def test_gas_work(test_args: Args) -> None:
    result = work_done_by_gas_at_constant_pressure.calculate_work(test_args.p, test_args.v_1,
        test_args.v_2)
    assert_equal(result, 2.52 * units.joule)


def test_bad_pressure(test_args: Args) -> None:
    bp = Quantity(2 * units.coulomb)
    with raises(errors.UnitsError):
        work_done_by_gas_at_constant_pressure.calculate_work(bp, test_args.v_1, test_args.v_2)
    with raises(TypeError):
        work_done_by_gas_at_constant_pressure.calculate_work(2, test_args.v_1, test_args.v_2)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_done_by_gas_at_constant_pressure.calculate_work(test_args.p, vb, test_args.v_2)
    with raises(errors.UnitsError):
        work_done_by_gas_at_constant_pressure.calculate_work(test_args.p, test_args.v_1, vb)
    with raises(TypeError):
        work_done_by_gas_at_constant_pressure.calculate_work(test_args.p, 20, test_args.v_2)
    with raises(TypeError):
        work_done_by_gas_at_constant_pressure.calculate_work(test_args.p, test_args.v_1, 20)
