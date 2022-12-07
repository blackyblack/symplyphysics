from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.gravity import free_fall_acceleration_from_height as free_fall_acceleration

# G - universal gravity constant  6.672e-11 N*m^2/kg^2
# M - Earth mass constant         5.976e+24 kg
# R - Earth radius constant       6.371e+6 m
@fixture
def test_args():
    earth_mass = units.Quantity("earth_mass")
    SI.set_quantity_dimension(earth_mass, units.mass)
    SI.set_quantity_scale_factor(earth_mass, 5.976e+24 * units.kilogram)
    earth_radius = units.Quantity("earth_radius")
    SI.set_quantity_dimension(earth_radius, units.length)
    SI.set_quantity_scale_factor(earth_radius, 6.371e+6 * units.meter)
    height_from_surface = units.Quantity("height_from_surface")
    SI.set_quantity_dimension(height_from_surface, units.length)
    SI.set_quantity_scale_factor(height_from_surface, 1 * units.meter)
    Args = namedtuple('Args', ["earth_mass", "earth_radius", "height_from_surface"])
    return Args(earth_mass=earth_mass, earth_radius=earth_radius, height_from_surface=height_from_surface)

def test_basic_acceleration(test_args):
    result = free_fall_acceleration.calculate_acceleration(test_args.earth_mass, test_args.earth_radius, test_args.height_from_surface)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.acceleration)
    result_acceleration = convert_to(result, units.meter / units.second**2).subs(units.meter / units.second**2, 1).evalf(6)
    assert result_acceleration == approx(9.82316, 0.005)

def test_bad_earth_mass(test_args):
    emb = units.Quantity('emb')
    SI.set_quantity_dimension(emb, units.length)
    SI.set_quantity_scale_factor(emb, 1 * units.meter)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(emb, test_args.earth_radius, test_args.height_from_surface)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(100, test_args.earth_radius, test_args.height_from_surface)

def test_bad_earth_radius(test_args):
    erb = units.Quantity('erb')
    SI.set_quantity_dimension(erb, units.mass)
    SI.set_quantity_scale_factor(erb, 1 * units.kilogram)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, erb, test_args.height_from_surface)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, 100, test_args.height_from_surface)

def test_bad_height(test_args):
    hb = units.Quantity('hb')
    SI.set_quantity_dimension(hb, units.mass)
    SI.set_quantity_scale_factor(hb, 1 * units.kilogram)
    with raises(errors.UnitsError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, test_args.earth_radius, hb)
    with raises(TypeError):
        free_fall_acceleration.calculate_acceleration(test_args.earth_mass, test_args.earth_radius, 100)