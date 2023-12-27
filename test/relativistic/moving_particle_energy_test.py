import pytest
from symplyphysics.laws.relativistic.moving_particle_energy import calculate_particle_energy

from pytest import approx, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from sympy import Abs, sqrt
from sympy.physics.units import speed_of_light


def test_basic():
    e = calculate_particle_energy(Quantity(1 * units.kilogram), 
                                  Quantity(0.5 * speed_of_light))
    assert SI.get_dimension_system().equivalent_dims(e.dimension, units.energy)
    e_nondim = convert_to(e, 
                          units.kilogram * units.meter ** 2 / units.second ** 2).as_coeff_Mul()[0]
    c = convert_to(speed_of_light, units.meter / units.second).as_coeff_Mul()[0]
    check = c ** 2 * sqrt(1 / (1 - 0.5 ** 2))
    assert Abs(e_nondim - check) < 1e-16
    

def test_bad_mass():
    mb = Quantity(1 * units.coulomb)
    v = Quantity(0.5 * speed_of_light)
    with raises(errors.UnitsError):
        calculate_particle_energy(mb, v)
    with raises(TypeError):
        calculate_particle_energy(100, v)
