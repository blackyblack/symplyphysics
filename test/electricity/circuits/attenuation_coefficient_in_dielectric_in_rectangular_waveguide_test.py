from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (errors, units, Quantity, assert_equal,)

from symplyphysics.laws.electricity.circuits import attenuation_coefficient_in_dielectric_in_rectangular_waveguide as coefficient_law

## The resistance of the waveguide is 306.45 ohm. The resistance of the material filling the waveguide is equal to 254.167 ohm.
## the tangent of the dielectric loss angle is 1e-4. The wavelength is equal to 10 millimeter.
## The attenuation coefficient will be 0.0379 [1 / meter].
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["resistance_of_waveguide", "resistance_of_medium", "wavelength", "tangent_dielectric_loss_angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resistance_of_waveguide = Quantity(306.45 * units.ohm)
    resistance_of_medium = Quantity(254.167 * units.ohm)
    wavelength = Quantity(10 * units.millimeter)
    tangent_dielectric_loss_angle = 1e-4
    return Args(resistance_of_waveguide=resistance_of_waveguide, resistance_of_medium=resistance_of_medium, wavelength=wavelength, tangent_dielectric_loss_angle=tangent_dielectric_loss_angle)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, test_args.resistance_of_medium, test_args.wavelength,
        test_args.tangent_dielectric_loss_angle)
    assert_equal(result, 0.0379 * (1 / units.meter))


def test_bad_resistance(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_resistance, test_args.resistance_of_medium, test_args.wavelength, test_args.tangent_dielectric_loss_angle)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(100, test_args.resistance_of_medium, test_args.wavelength, test_args.tangent_dielectric_loss_angle)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, bad_resistance, test_args.wavelength, test_args.tangent_dielectric_loss_angle)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, 100, test_args.wavelength, test_args.tangent_dielectric_loss_angle)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, test_args.resistance_of_medium, bad_wavelength, test_args.tangent_dielectric_loss_angle)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, test_args.resistance_of_medium, 100, test_args.tangent_dielectric_loss_angle)


def test_bad_tangent_dielectric_loss_angle(test_args: Args) -> None:
    tangent_dielectric_loss_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.resistance_of_waveguide, test_args.resistance_of_medium, test_args.wavelength, tangent_dielectric_loss_angle)
