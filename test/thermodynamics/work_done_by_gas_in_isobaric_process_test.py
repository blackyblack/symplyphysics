from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    convert_to,
    SI,
)
from symplyphysics.laws.thermodynamics import work_done_by_gas_in_isobaric_process


@fixture(name="test_args")
def test_args_fixture():
    p = Quantity(2.52 * units.pascal)
    v_1 = Quantity(1 * units.meter**3)
    v_2 = Quantity(2 * units.meter**3)
    Args = namedtuple("Args", ["p", "v_1", "v_2"])
    return Args(p=p, v_1=v_1, v_2=v_2)


def test_gas_work(test_args):
    result = work_done_by_gas_in_isobaric_process.calculate_work(test_args.p, test_args.v_1, test_args.v_2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).evalf(4)
    assert result_work == approx(2.52, 0.01)


def test_bad_pressue(test_args):
    bp = Quantity(2 * units.coulomb)
    with raises(errors.UnitsError):
        work_done_by_gas_in_isobaric_process.calculate_work(
            bp, test_args.v_1, test_args.v_2)
    with raises(errors.UnitsError):
        work_done_by_gas_in_isobaric_process.calculate_work(
            2, test_args.v_1, test_args.v_2)


def test_bad_volume(test_args):
    vb = Quantity(1 * units.coulomb)
    vb2 = Quantity(2 * units.coulomb)
    with raises(errors.UnitsError):
        work_done_by_gas_in_isobaric_process.calculate_work(
            test_args.p, vb, vb2)
    with raises(errors.UnitsError):
        work_done_by_gas_in_isobaric_process.calculate_work(
            test_args.p, 20, 30)
