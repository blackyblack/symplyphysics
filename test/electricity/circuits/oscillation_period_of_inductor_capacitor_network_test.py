from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import oscillation_period_of_inductor_capacitor_network as lc

# Description
## Assert we have a capacitor with 1 farad capacitance and 1 henry inductor.
## According to Tomson's formula oscillation period of this circuit should be 6.28 seconds

Args = namedtuple("Args", ["L", "C"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    L = Quantity(1 * units.henry)
    C = Quantity(1 * units.farad)
    return Args(L=L, C=C)


def test_basic_period(test_args: Args) -> None:
    result = lc.calculate_oscillation_period(test_args.L, test_args.C)
    assert_equal(result, 6.28 * units.second)


def test_bad_inductance(test_args: Args) -> None:
    Lb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        lc.calculate_oscillation_period(Lb, test_args.C)
    with raises(TypeError):
        lc.calculate_oscillation_period(100, test_args.C)


def test_bad_capacity(test_args: Args) -> None:
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        lc.calculate_oscillation_period(test_args.L, Cb)
    with raises(TypeError):
        lc.calculate_oscillation_period(test_args.L, 100)
