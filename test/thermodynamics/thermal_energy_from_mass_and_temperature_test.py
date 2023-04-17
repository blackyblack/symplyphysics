from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.thermodynamics import thermal_energy_from_mass_and_temperature as amount_energy

# How much energy does it take to heat some volume of water weighing 0.5 kilograms то 50 Celsius degree?
# Specific heat capacity of water is 4200 J/kg*K , ignore losses.
@fixture
def test_args():
    C = Quantity(units.energy / (units.mass * units.temperature), 4200 * units.joule / (units.kilogram * units.kelvin))
    m = Quantity(units.mass, 0.5 * units.kilogram)
    t1 = Quantity(units.temperature, 273 * units.kelvin)
    t2 = Quantity(units.temperature, 323 * units.kelvin)
    Args = namedtuple("Args", ["C", "m", "t1", "t2"])
    return Args(C=C, m=m, t1=t1, t2=t2)

def test_basic_amount(test_args):
    result = amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, test_args.t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(7)
    assert result_energy == approx(105000.1, 0.000001)

def test_bad_specific_heat(test_args):
    Cb = Quantity(units.mass)
    with raises(errors.UnitsError):
         amount_energy.calculate_amount_energy(Cb, test_args.m, test_args.t2, test_args.t1)
    with raises(TypeError):
         amount_energy.calculate_amount_energy(100, test_args.m, test_args.t2, test_args.t1)

def test_bad_body_mass(test_args):
    mb = Quantity(units.temperature)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, mb, test_args.t2, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, 100, test_args.t2, test_args.t1)

def test_bad_temperature(test_args):
    tb = Quantity(units.mass)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, tb, test_args.t1)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, 100, test_args.t1)
    with raises(errors.UnitsError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, tb)
    with raises(TypeError):
        amount_energy.calculate_amount_energy(test_args.C, test_args.m, test_args.t2, 100)
