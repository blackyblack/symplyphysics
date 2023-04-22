# Description
## Assume we have a capacitor with 6 Coulombs of charge in it and having 3 volts on it's terminals.
## According to the definition of Capacitance this capacitor should have C = 2 Farads.

from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import capacitance_from_charge_and_voltage as capacitance_def


@fixture
def test_args():
    C = Quantity(2 * units.farad)
    Q = Quantity(6 * units.coulomb)
    U = Quantity(3 * units.volts)
    Args = namedtuple("Args", ["C", "Q", "U"])
    return Args(C=C, Q=Q, U=U)


def test_basic_capacitance(test_args):
    result = capacitance_def.calculate_capacitance(test_args.Q, test_args.U)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.capacitance)
    result_cap = convert_to(result, capacitance_def.definition_units_SI).subs(units.farad,
        1).evalf(4)
    assert result_cap == approx(2, 0.01)


def test_bad_charge(test_args):
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(Qb, test_args.U)
    with raises(TypeError):
        capacitance_def.calculate_capacitance(100, test_args.U)


def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitance_def.calculate_capacitance(test_args.Q, Vb)
    with raises(TypeError):
        capacitance_def.calculate_capacitance(test_args.Q, 100)
