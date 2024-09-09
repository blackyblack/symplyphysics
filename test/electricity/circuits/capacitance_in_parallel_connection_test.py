from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import capacitance_in_parallel_connection as parallel_capacitor

# Description
## Assert we are having two capacitors with 2 and 3 Farads of capacitance.
## Resulting capacitance should be 5 Farads.

Args = namedtuple("Args", ["C1", "C2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    C1 = Quantity(2 * units.farad)
    C2 = Quantity(3 * units.farad)
    return Args(C1=C1, C2=C2)


def test_basic_capacitance(test_args: Args) -> None:
    result = parallel_capacitor.calculate_parallel_capacitance([test_args.C1, test_args.C2])
    assert_equal(result, 5 * units.farad)


def test_three_capacitors_array(test_args: Args) -> None:
    C3 = Quantity(5 * units.farad)
    result = parallel_capacitor.calculate_parallel_capacitance([test_args.C1, test_args.C2, C3])
    assert_equal(result, 10 * units.farad)


def test_bad_capacitance(test_args: Args) -> None:
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_capacitor.calculate_parallel_capacitance([Cb, test_args.C2])
    with raises(TypeError):
        parallel_capacitor.calculate_parallel_capacitance([100, test_args.C2])
    with raises(errors.UnitsError):
        parallel_capacitor.calculate_parallel_capacitance([test_args.C1, Cb])
    with raises(TypeError):
        parallel_capacitor.calculate_parallel_capacitance([test_args.C1, 100])
    with raises(errors.UnitsError):
        parallel_capacitor.calculate_parallel_capacitance([Cb, Cb])
    with raises(TypeError):
        parallel_capacitor.calculate_parallel_capacitance([100, 100])
    with raises(TypeError):
        parallel_capacitor.calculate_parallel_capacitance(test_args.C1)
