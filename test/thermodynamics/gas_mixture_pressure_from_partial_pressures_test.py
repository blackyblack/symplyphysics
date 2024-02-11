from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.thermodynamics import gas_mixture_pressure_from_partial_pressures as gas_mixture_pressure

Args = namedtuple("Args", ["p1", "p2", "p3"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p1 = Quantity(200 * units.pascal)
    p2 = Quantity(100 * units.pascal)
    p3 = Quantity(400 * units.pascal)
    return Args(p1=p1, p2=p2, p3=p3)


def test_basic_gas_mixture_pressure(test_args: Args) -> None:
    result = gas_mixture_pressure.calculate_total_pressure(
        [test_args.p1, test_args.p2, test_args.p3])
    assert_equal(result, 700 * units.pascal)


def test_array_with_bad_element(test_args: Args) -> None:
    bad_element = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        gas_mixture_pressure.calculate_total_pressure([test_args.p1, bad_element])
    with raises(TypeError):
        gas_mixture_pressure.calculate_total_pressure([test_args.p1, 100])
    with raises(errors.UnitsError):
        gas_mixture_pressure.calculate_total_pressure([bad_element, test_args.p2])
    with raises(TypeError):
        gas_mixture_pressure.calculate_total_pressure([100, test_args.p2])
    with raises(errors.UnitsError):
        gas_mixture_pressure.calculate_total_pressure([bad_element, bad_element])
    with raises(TypeError):
        gas_mixture_pressure.calculate_total_pressure([100, 100])
    with raises(TypeError):
        gas_mixture_pressure.calculate_total_pressure(test_args.p1)
