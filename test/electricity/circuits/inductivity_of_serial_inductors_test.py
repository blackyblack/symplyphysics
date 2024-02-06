from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import inductivity_of_serial_inductors as serial_inductor

# Description
## Assert we have two inductors with 1H and 2H inductances connected in series.
## Resulting inductance should be sum of 1 and 2 - 3H.


@fixture(name="test_args")
def test_args_fixture():
    L1 = Quantity(1 * units.henry)
    L2 = Quantity(2 * units.henry)
    Args = namedtuple("Args", ["L1", "L2"])
    return Args(L1=L1, L2=L2)


def test_basic_inductivity(test_args):
    result = serial_inductor.calculate_serial_inductance([test_args.L1, test_args.L2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.inductance)
    result_inductance = convert_to(result, units.henry).evalf(3)
    assert_approx(result_inductance, 3)


def test_three_inductors_array(test_args):
    L3 = Quantity(3 * units.henry)
    result = serial_inductor.calculate_serial_inductance([test_args.L1, test_args.L2, L3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.inductance)
    result_inductance = convert_to(result, units.henry).evalf(3)
    assert_approx(result_inductance, 6)


def test_bad_inductivity(test_args):
    Lb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        serial_inductor.calculate_serial_inductance([Lb, test_args.L2])
    with raises(TypeError):
        serial_inductor.calculate_serial_inductance([100, test_args.L2])
    with raises(errors.UnitsError):
        serial_inductor.calculate_serial_inductance([test_args.L1, Lb])
    with raises(TypeError):
        serial_inductor.calculate_serial_inductance([test_args.L1, 100])
    with raises(errors.UnitsError):
        serial_inductor.calculate_serial_inductance([Lb, Lb])
    with raises(TypeError):
        serial_inductor.calculate_serial_inductance([100, 100])
    with raises(TypeError):
        serial_inductor.calculate_serial_inductance(test_args.L1)
