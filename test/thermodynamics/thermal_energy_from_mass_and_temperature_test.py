from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as amount_energy

# How much energy does it take to heat some volume of water weighing 0.5 kilograms то 50 Celsius degree?
# Specific heat capacity of water is 4200 J/kg*K , ignore losses.
@fixture
def test_args():
    C = units.Quantity('C')
    SI.set_quantity_dimension(C, units.energy / (units.mass * units.temperature))
    SI.set_quantity_scale_factor(C, 4200 * units.joule / (units.kilogram * units.kelvin))
    m = units.Quantity('m')
    SI.set_quantity_dimension(m, units.mass)
    SI.set_quantity_scale_factor(m, 0.5 * units.kilogram)
    t1 = units.Quantity('t1')
    SI.set_quantity_dimension(t1, units.temperature)
    SI.set_quantity_scale_factor(t1, 273 * units.kelvin)
    t2 = units.Quantity('t2')
    SI.set_quantity_dimension(t2, units.temperature)
    SI.set_quantity_scale_factor(t2, 323 * units.kelvin)
    Args = namedtuple('Args', ['C', 'm', 't1', 't2'])
    return Args(C=C, m=m, t1=t1, t2=t2)

def test_basic_amount(test_args):
    result = amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, test_args.t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(7)
    assert result_energy == approx(105000.1, 0.000001)


def test_bad_specific_heat(test_args):
    bC = units.Quantity('bC')
    SI.set_quantity_dimension(bC, units.mass)
    SI.set_quantity_scale_factor(bC, 1 * units.kilogram)
    with raises(errors.UnitsError):
         amount_energy.calculate_amount_energy(bC, test_args.m, test_args.t2, test_args.t1)
    with raises(TypeError):
         amount_energy.calculate_amount_energy(100, test_args.m, test_args.t2, test_args.t1)

def test_bad_body_mass(test_args):
    bm = units.Quantity('bm')
    SI.set_quantity_dimension(bm, units.temperature)
    SI.set_quantity_scale_factor(bm, 1 * units.kelvin)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, bm, test_args.t2, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, 100, test_args.t2, test_args.t1)

def test_bad_temperature(test_args):
    bt = units.Quantity('bt')
    SI.set_quantity_dimension(bt, units.mass)
    SI.set_quantity_scale_factor(bt, 1 * units.kilogram)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, bt, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, 100, test_args.t1)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, bt)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, 100)