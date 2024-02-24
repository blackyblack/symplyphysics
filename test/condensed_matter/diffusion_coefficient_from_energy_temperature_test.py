from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.condensed_matter import diffusion_coefficient_from_energy_temperature as diffusion_law

# Description
## Consider a dopant of phosphorus in silicon. The diffusion constant will be 3.85 [centimeter^2 / second],
## the activation energy will be 3.66 electronvolt. At a temperature of 1000 kelvin, the diffusion
## coefficient will be 1.38e-18 [centimeter^2 / second].
## https://elib.gsu.by/bitstream/123456789/5368/1/Лаб%20№4%20Моделирование%20процесса%20диффузии.pdf

Args = namedtuple("Args", ["energy", "diffusion_constant", "temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    energy = Quantity(3.66 * units.electronvolt)
    diffusion_constant = Quantity(3.85 * (units.centimeter**2 / units.second))
    temperature = Quantity(1000 * units.kelvin)

    return Args(energy=energy,
        diffusion_constant=diffusion_constant,
        temperature=temperature)


def test_basic_diffusion_coefficient(test_args: Args) -> None:
    result = diffusion_law.calculate_diffusion_coefficient(test_args.energy,
        test_args.diffusion_constant, test_args.temperature)
    assert_equal(result, 1.38e-18 * (units.centimeter**2 / units.second))


def test_bad_energy(test_args: Args) -> None:
    energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_law.calculate_diffusion_coefficient(energy, test_args.diffusion_constant,
            test_args.temperature)
    with raises(TypeError):
        diffusion_law.calculate_diffusion_coefficient(100, test_args.diffusion_constant,
            test_args.temperature)


def test_bad_diffusion_constant(test_args: Args) -> None:
    diffusion_constant = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_law.calculate_diffusion_coefficient(test_args.energy, diffusion_constant,
            test_args.temperature)
    with raises(TypeError):
        diffusion_law.calculate_diffusion_coefficient(test_args.energy, 100,
            test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        diffusion_law.calculate_diffusion_coefficient(test_args.energy,
            test_args.diffusion_constant, temperature)
    with raises(TypeError):
        diffusion_law.calculate_diffusion_coefficient(test_args.energy,
            test_args.diffusion_constant, 100)
