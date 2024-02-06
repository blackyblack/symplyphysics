from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, errors, units, Quantity)

from symplyphysics.laws.thermodynamics import diffusion_flux_from_diffusion_coefficient_and_concentration_gradient as diffusion_flux


@fixture(name="test_args")
def test_args_fixture():
    diffusion_coefficient = Quantity(0.001 * units.meter**2 / units.second)
    concentration_start = Quantity(10 * units.mole / units.meter**3)
    concentration_end = Quantity(15 * units.mole / units.meter**3)
    position = Quantity(5 * units.meter)
    Args = namedtuple("Args",
        ["diffusion_coefficient", "concentration_start", "concentration_end", "position"])
    return Args(diffusion_coefficient=diffusion_coefficient,
        concentration_start=concentration_start,
        concentration_end=concentration_end,
        position=position)


def test_basic_law(test_args):
    result = diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient,
        test_args.concentration_start, test_args.concentration_end, test_args.position)
    assert_equal(result, -0.001 * units.mole / (units.meter**2 * units.second))


def test_diffusion_coefficient(test_args):
    bad_diffusion_coefficient = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        diffusion_flux.calculate_diffusion_flux(bad_diffusion_coefficient,
            test_args.concentration_start, test_args.concentration_end, test_args.position)
    with raises(TypeError):
        diffusion_flux.calculate_diffusion_flux(100, test_args.concentration_start,
            test_args.concentration_end, test_args.position)


def test_bad_concentration(test_args):
    bad_concentration = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient, bad_concentration,
            test_args.concentration_end, test_args.position)
    with raises(TypeError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient, 100,
            test_args.concentration_end, test_args.position)
    with raises(errors.UnitsError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient,
            test_args.concentration_start, bad_concentration, test_args.position)
    with raises(TypeError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient,
            test_args.concentration_start, 100, test_args.position)


def test_bad_position(test_args):
    bad_position = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient,
            test_args.concentration_start, test_args.concentration_end, bad_position)
    with raises(TypeError):
        diffusion_flux.calculate_diffusion_flux(test_args.diffusion_coefficient,
            test_args.concentration_start, test_args.concentration_end, 100)
