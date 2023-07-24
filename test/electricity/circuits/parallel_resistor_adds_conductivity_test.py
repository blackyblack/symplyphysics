from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import parallel_resistor_adds_conductivity as parallel_resistor

# Description
## Assert we have two resistors with 2 Ohm and 4 Ohm impedances.
## Accordind to calculator (https://www.chipdip.ru/calc/parallel-resistors) resulting resistance should be 1.333 Ohm.


@fixture(name="test_args")
def test_args_fixture():
    R1 = Quantity(2 * units.ohm)
    R2 = Quantity(4 * units.ohm)
    Args = namedtuple("Args", ["R1", "R2"])
    return Args(R1=R1, R2=R2)


def test_basic_resistance(test_args):
    result = parallel_resistor.calculate_resistance(test_args.R1, test_args.R2)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.impedance)
    result_resistance = convert_to(result, units.ohm).evalf(3)
    assert result_resistance == approx(1.333, 0.001)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_resistance(Rb, test_args.R2)
    with raises(TypeError):
        parallel_resistor.calculate_resistance(100, test_args.R2)
    with raises(errors.UnitsError):
        parallel_resistor.calculate_resistance(test_args.R1, Rb)
    with raises(TypeError):
        parallel_resistor.calculate_resistance(test_args.R1, 100)
    
