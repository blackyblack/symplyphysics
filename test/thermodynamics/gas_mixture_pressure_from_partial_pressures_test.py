from collections import namedtuple
from pytest import approx, fixture, raises
from sympy import S
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.thermodynamics import gas_mixture_pressure_from_partial_pressures as gas_mixture_pressure


@fixture(name="test_args")
def test_args_fixture():
    p1 = Quantity(200 * units.pascal)
    p2 = Quantity(100 * units.pascal)
    p3 = Quantity(400 * units.pascal)
    Args = namedtuple("Args", ["p1", "p2", "p3"])
    return Args(
        p1=p1,
        p2=p2,
        p3=p3
    )


def test_basic_gas_mixture_pressure(test_args):
    result = gas_mixture_pressure.calculate_total_pressure([test_args.p1, test_args.p2, test_args.p3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(3)
    assert result_pressure == approx(700, 0.001)


def test_array_with_bad_element(test_args):
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
