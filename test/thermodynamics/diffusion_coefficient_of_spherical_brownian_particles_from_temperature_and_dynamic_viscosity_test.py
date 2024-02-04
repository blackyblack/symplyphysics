from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to
)

from symplyphysics.laws.thermodynamics import diffusion_coefficient_of_spherical_brownian_particles_from_temperature_and_dynamic_viscosity as diffusion_coefficient


@fixture(name="test_args")
def test_args_fixture():
    temperature = Quantity(298 * units.kelvin)
    particle_radius = Quantity(1e-8 * units.meter)
    dynamic_viscosity = Quantity(1.9e-7 * units.pascal * units.second)
    Args = namedtuple("Args", ["particle_radius", "dynamic_viscosity", "temperature"])
    return Args(
        particle_radius=particle_radius,
        dynamic_viscosity=dynamic_viscosity,
        temperature=temperature
    )


def test_basic_law(test_args):
    result = diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, test_args.dynamic_viscosity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.area / units.time)
    result_diffusion_coefficient = convert_to(result, units.meter**2 / units.second).evalf(3)
    assert result_diffusion_coefficient == approx(1.15e-7, abs=1e-4)


def test_bad_temperature(test_args):
    bad_temperature = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(bad_temperature, test_args.particle_radius, test_args.dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(100, test_args.particle_radius, test_args.dynamic_viscosity)


def test_bad_particle_radius(test_args):
    bad_particle_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, bad_particle_radius, test_args.dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, 100, test_args.dynamic_viscosity)


def test_bad_dynamic_viscosity(test_args):
    bad_dynamic_viscosity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, bad_dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, 100)
