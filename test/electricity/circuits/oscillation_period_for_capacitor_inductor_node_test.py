# Description
## Assert we have a capacitor with 1 farad capacitance and 1 henry inductor.
## Accordind to Tomson's formula oscillation period of this circuit should be 6.28 seconds

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.electricity.circuits import oscillation_period_for_capacitor_inductor_node as lc

@fixture
def test_args():
    L = Quantity(units.inductance, 1 * units.henry)
    C = Quantity(units.capacitance, 1 * units.farad)
    Args = namedtuple("Args", ["L", "C"])
    return Args(L=L, C=C)

def test_basic_period(test_args):
    result = lc.calculate_oscillation_period(test_args.L, test_args.C)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.time)
    result_voltage = convert_to(result, units.second).subs(units.second, 1).evalf(2)
    assert result_voltage == approx(6.28, 0.01)

def test_bad_inductance(test_args):
    Lb = Quantity(units.length)
    with raises(errors.UnitsError):
        lc.calculate_oscillation_period(Lb, test_args.C)
    with raises(TypeError):
        lc.calculate_oscillation_period(100, test_args.C)

def test_bad_capacity(test_args):
    Cb = Quantity(units.length)
    with raises(errors.UnitsError):
        lc.calculate_oscillation_period(test_args.L, Cb)
    with raises(TypeError):
        lc.calculate_oscillation_period(test_args.L, 100)
