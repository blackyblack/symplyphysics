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
from symplyphysics.laws.electricity.circuits import resistivity_of_parallel_resistors as parallel_resistor

# Description
## Assert we have three resistors with 2 Ohm, 2 Ohm, and 4 Ohm resistances in parallel.
## Resulting resistance should be 0.8 Ohm.
## If we add another resistor of 12 Ohm, the resulting resistance should be 0.75 Ohm.


@fixture(name="test_args")
def test_args_fixture():
    R1 = Quantity(2 * units.ohm)
    R2 = Quantity(2 * units.ohm)
    R3 = Quantity(4 * units.ohm)
    Args = namedtuple("Args", ["R1", "R2", "R3"])
    return Args(R1=R1, R2=R2, R3=R3)


def test_basic_resistance(test_args):
    result = parallel_resistor.calculate_parallel_resistance(
        [test_args.R1, test_args.R2, test_args.R3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_resistance = convert_to(result, units.ohm).evalf(3)
    assert_approx(result_resistance, 0.8)


def test_four_resistors_array(test_args):
    R4 = Quantity(12 * units.ohm)
    result = parallel_resistor.calculate_parallel_resistance(
        [test_args.R1, test_args.R2, test_args.R3, R4])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_resistance = convert_to(result, units.ohm).evalf(3)
    assert_approx(result_resistance, 0.75)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_parallel_resistance([Rb, test_args.R2])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_resistance([100, test_args.R2])
    with raises(errors.UnitsError):
        parallel_resistor.calculate_parallel_resistance([test_args.R1, Rb])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_resistance([test_args.R1, 100])
    with raises(TypeError):
        parallel_resistor.calculate_parallel_resistance(test_args.R1)
