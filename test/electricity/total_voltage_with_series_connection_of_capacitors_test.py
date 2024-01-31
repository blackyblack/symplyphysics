from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.electricity import total_voltage_with_series_connection_of_capacitors as total_voltage


@fixture(name="test_args")
def test_args_fixture():
    u1 = Quantity(2 * units.volt)
    u2 = Quantity(5 * units.volt)
    u3 = Quantity(4 * units.volt)
    Args = namedtuple("Args", ["u1", "u2", "u3"])
    return Args(
        u1=u1,
        u2=u2,
        u3=u3
    )


def test_basic_total_voltage(test_args):
    result = total_voltage.calculate_total_voltage([test_args.u1, test_args.u2, test_args.u3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result_voltage = convert_to(result, units.volt).evalf(3)
    assert result_voltage == approx(11, 0.001)


def test_array_with_bad_element(test_args):
    bad_element = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        total_voltage.calculate_total_voltage([test_args.u1, bad_element])
    with raises(TypeError):
        total_voltage.calculate_total_voltage([test_args.u1, 100])
    with raises(errors.UnitsError):
        total_voltage.calculate_total_voltage([bad_element, test_args.u2])
    with raises(TypeError):
        total_voltage.calculate_total_voltage([100, test_args.u2])
    with raises(errors.UnitsError):
        total_voltage.calculate_total_voltage([bad_element, bad_element])
    with raises(TypeError):
        total_voltage.calculate_total_voltage([100, 100])
    with raises(TypeError):
        total_voltage.calculate_total_voltage(test_args.u1)
