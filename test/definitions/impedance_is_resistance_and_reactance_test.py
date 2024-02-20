from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import impedance_is_resistance_and_reactance as impedance_def

# Description
## If the object has 10 Ohm resistance in series with capacitor with reactance of -1.59 Ohm,
## the impedance magnitude of the circuit is 10.12 Ohm.
## Taken from https://wiraelectrical.com/impedance-and-admittance/

Args = namedtuple("Args", ["R", "X"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    R = Quantity(10 * units.ohm)
    X = Quantity(-1.59 * units.ohm)
    return Args(R=R, X=X)


def test_basic_impedance(test_args: Args) -> None:
    result = impedance_def.calculate_impedance_magnitude(test_args.R, test_args.X)
    assert_equal(result, 10.12 * units.ohm)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        impedance_def.calculate_impedance_magnitude(Rb, test_args.X)
    with raises(TypeError):
        impedance_def.calculate_impedance_magnitude(100, test_args.X)


def test_bad_reactance(test_args: Args) -> None:
    Xb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        impedance_def.calculate_impedance_magnitude(test_args.R, Xb)
    with raises(TypeError):
        impedance_def.calculate_impedance_magnitude(test_args.R, 100)
