from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    prefixes,
    units,
    Quantity,
    assert_equal
)

from symplyphysics.laws.thermodynamics import diffusion_coefficient_of_spherical_brownian_particles_from_temperature_and_dynamic_viscosity as diffusion_coefficient

# Using as reference https://studfile.net/preview/6307326/page:2/
# Reference calculations are very inaccurate. Recalculated manually.

Args = namedtuple("Args", ["particle_radius", "dynamic_viscosity", "temperature"])

@fixture(name="test_args")
def test_args_fixture() -> Args:
    temperature = Quantity(300 * units.kelvin)
    particle_radius = Quantity(3 * units.angstrom)
    # 0.01 poise
    dynamic_viscosity = Quantity(0.01 * units.gram / units.second / units.centimeter)
    return Args(
        particle_radius=particle_radius,
        dynamic_viscosity=dynamic_viscosity,
        temperature=temperature
    )


def test_basic_law(test_args: Args) -> None:
    result = diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, test_args.dynamic_viscosity)
    assert_equal(result, 7.32e-6 * (prefixes.centi * units.meter)**2 / units.second)


def test_bad_temperature(test_args: Args) -> None:
    bad_temperature = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(bad_temperature, test_args.particle_radius, test_args.dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(100, test_args.particle_radius, test_args.dynamic_viscosity)


def test_bad_particle_radius(test_args: Args) -> None:
    bad_particle_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, bad_particle_radius, test_args.dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, 100, test_args.dynamic_viscosity)


def test_bad_dynamic_viscosity(test_args: Args) -> None:
    bad_dynamic_viscosity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, bad_dynamic_viscosity)
    with raises(TypeError):
        diffusion_coefficient.calculate_diffusion_coefficient(test_args.temperature, test_args.particle_radius, 100)
