from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import inductance_in_serial_connection as serial_inductor

# Description
## Assert we have two inductors with 1H and 2H inductances connected in series.
## Resulting inductance should be sum of 1 and 2 - 3H.

Args = namedtuple("Args", ["L1", "L2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    L1 = Quantity(1 * units.henry)
    L2 = Quantity(2 * units.henry)
    return Args(L1=L1, L2=L2)


def test_basic_inductivity(test_args: Args) -> None:
    result = serial_inductor.calculate_serial_inductance([test_args.L1, test_args.L2])
    assert_equal(result, 3 * units.henry)


def test_three_inductors_array(test_args: Args) -> None:
    L3 = Quantity(3 * units.henry)
    result = serial_inductor.calculate_serial_inductance([test_args.L1, test_args.L2, L3])
    assert_equal(result, 6 * units.henry)


def test_bad_inductivity(test_args: Args) -> None:
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
