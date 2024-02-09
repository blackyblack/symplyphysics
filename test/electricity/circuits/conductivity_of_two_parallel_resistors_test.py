from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity.circuits import conductivity_of_two_parallel_resistors as parallel_resistor

# Description
## Assert we have two resistors with 2 Ohm and 4 Ohm impedances.
## Accordind to calculator (https://www.chipdip.ru/calc/parallel-resistors) resulting resistance should be 1.333 Ohm.

Args = namedtuple("Args", ["R1", "R2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    R1 = Quantity(2 * units.ohm)
    R2 = Quantity(4 * units.ohm)
    return Args(R1=R1, R2=R2)


def test_basic_resistance(test_args: Args) -> None:
    result = parallel_resistor.calculate_resistance(test_args.R1, test_args.R2)
    assert_equal(result, 1.333 * units.ohm)


def test_bad_resistance(test_args: Args) -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_resistance(Rb, test_args.R2)
    with raises(TypeError):
        parallel_resistor.calculate_resistance(100, test_args.R2)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_resistance(test_args.R1, Rb)
    with raises(TypeError):
        parallel_resistor.calculate_resistance(test_args.R1, 100)
