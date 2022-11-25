from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

#How much energy can a household electric kettle produce to heat water in one minute
#at 220 volts and a heater resistance of 36 ohms?
@fixture
def test_args():
    U = units.Quantity('U')
    SI.set_quantity_dimension(U, units.voltage)
    SI.set_quantity_scale_factor(U, 220 * units.volt)
    t = units.Quantity('t')
    SI.set_quantity_dimension(t, units.time)
    SI.set_quantity_scale_factor(t, 60 * units.second)
    R = units.Quantity('R')
    SI.set_quantity_dimension(R, units.impedance)
    SI.set_quantity_scale_factor(R, 36 * units.ohm)
    Args = namedtuple('Args', ['U', 't', 'R'])
    return Args(U=U, t=t, R=R)

def test_basic_amount(test_args):
    result = joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(6)
    assert result_energy == approx(80666.6, 0.000001)

def test_bad_voltage(test_args):
    bU = units.Quantity('bU')
    SI.set_quantity_dimension(bU, units.time)
    SI.set_quantity_scale_factor(bU, 1 * units.second)
    with raises(errors.UnitsError):
         joule_lenz_law.calculate_amount_energy(bU, test_args.t, test_args.R)
    with raises(TypeError):
         joule_lenz_law.calculate_amount_energy(100, test_args.t, test_args.R)

def test_bad_time(test_args):
    bt = units.Quantity('bt')
    SI.set_quantity_dimension(bt, units.voltage)
    SI.set_quantity_scale_factor(bt, 1 * units.volt)
    with raises(errors.UnitsError):
         joule_lenz_law.calculate_amount_energy(test_args.U, bt, test_args.R)
    with raises(TypeError):
         joule_lenz_law.calculate_amount_energy(test_args.U, 100, test_args.R)

def test_bad_resistance(test_args):
    bR = units.Quantity('bR')
    SI.set_quantity_dimension(bR, units.time)
    SI.set_quantity_scale_factor(bR, 1 * units.second)
    with raises(errors.UnitsError):
         joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, bR)
    with raises(TypeError):
         joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, 100)