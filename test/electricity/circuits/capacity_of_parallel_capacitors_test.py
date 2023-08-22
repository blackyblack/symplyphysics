from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity.circuits import capacity_of_parallel_capacitors as parallel_capacitor

# Description
## Assert we are having two capacitors with 2 and 3 Farads of capacitance.
## Resulting capacitance should be 5 Farads.


@fixture(name="test_args")
def test_args_fixture():
    C1 = Quantity(2 * units.farad)
    C2 = Quantity(3 * units.farad)
    Args = namedtuple("Args", ["C1", "C2"])
    return Args(C1=C1, C2=C2)


def test_basic_capacitance(test_args):
    result = parallel_capacitor.calculate_parallel_capacitance([test_args.C1, test_args.C2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)
    result_conductance = convert_to(result, units.farad).evalf(3)
    assert result_conductance == approx(5, 0.001)


def test_three_capacitors_array(test_args):
    C3 = Quantity(5 * units.farad)
    result = parallel_capacitor.calculate_parallel_capacitance([test_args.C1, test_args.C2, C3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)
    result_capacitance = convert_to(result, units.farad).evalf(3)
    assert result_capacitance == approx(10, 0.01)


def test_bad_capacitance(test_args):
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
