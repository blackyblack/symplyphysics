from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.dynamics import kinetic_energy_from_moment_of_inertia_and_angular_velocity as kinetic_energy_law
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type

# Description
## If we have a stone with 2kg*m^2 of moment of inertia spinning with 3rad/s angular velocity, this stone should bear 9 joules of kinetic energy.

@fixture
def test_args():
    I = units.Quantity('I')
    SI.set_quantity_dimension(I, units.mass * units.length**2)
    SI.set_quantity_scale_factor(I, 2 * units.kilogram * units.meter**2)

    w = units.Quantity('w')
    SI.set_quantity_dimension(w, angle_type / units.time)
    SI.set_quantity_scale_factor(w, 3 * units.radian / units.second)

    Args = namedtuple('Args', ['I', 'w'])
    return Args(I = I, w = w)

def test_basic_energy(test_args):
    result = kinetic_energy_law.calculate_energy(test_args.I, test_args.w)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_energy == approx(9.0, 0.01)

def test_bad_inertia_moment(test_args):
    Ib = units.Quantity('Ib')
    SI.set_quantity_dimension(Ib, units.length)
    SI.set_quantity_scale_factor(Ib, 1 * units.meter)

    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(Ib, test_args.w)

    with raises(TypeError):
        kinetic_energy_law.calculate_energy(100, test_args.w)

def test_bad_velocity(test_args):
    wb = units.Quantity('wb')
    SI.set_quantity_dimension(wb, units.length)
    SI.set_quantity_scale_factor(wb, 1 * units.meter)

    with raises(errors.UnitsError):
        kinetic_energy_law.calculate_energy(test_args.I, wb)

    with raises(TypeError):
        kinetic_energy_law.calculate_energy(test_args.I, 100)