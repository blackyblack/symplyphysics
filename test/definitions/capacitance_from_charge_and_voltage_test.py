# Description
## Assume we have a capacitor with 6 Coulombs of charge in it and having 3 volts on it's terminals.
## According to the definition of Capacitance this capacitor should have C = 2 Farads.

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors, solve
)
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_def

@fixture
def test_args():
    C = units.Quantity('C')
    SI.set_quantity_dimension(C, units.capacitance)    
    SI.set_quantity_scale_factor(C, 2 * units.farad)
    Q = units.Quantity('Q')
    SI.set_quantity_dimension(Q, units.charge)
    SI.set_quantity_scale_factor(Q, 6 * units.coulomb)
    U = units.Quantity('U')
    SI.set_quantity_dimension(U, units.voltage)
    SI.set_quantity_scale_factor(U, 3 * units.volts)

    Args = namedtuple('Args', ['C', 'Q', 'U'])
    return Args(C = C, Q = Q, U = U)

def test_basic_capacitance(test_args):
    result = capacitance_def.calculate_capacitance(test_args.Q, test_args.U)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)
    result_cap = convert_to(result, units.farad).subs(units.farad, 1).evalf(4)
    assert result_cap == approx(2, 0.01)


def test_bad_charge(test_args):
    Qb = units.Quantity('Qb')
    SI.set_quantity_dimension(Qb, units.length)
    SI.set_quantity_scale_factor(Qb, 1 * units.meter)

    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(Qb, test_args.U)

    with raises(TypeError):
        capacitance_def.calculate_capacitance(100, test_args.U)

def test_bad_voltage(test_args):
    Vb = units.Quantity('Vb')
    SI.set_quantity_dimension(Vb, units.length)
    SI.set_quantity_scale_factor(Vb, 1 * units.meter)

    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(test_args.Q, Vb)

    with raises(TypeError):
        capacitance_def.calculate_capacitance(test_args.Q, 100)
