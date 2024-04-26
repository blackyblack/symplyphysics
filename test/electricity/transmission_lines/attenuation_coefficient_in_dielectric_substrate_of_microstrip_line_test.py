from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines import attenuation_coefficient_in_dielectric_substrate_of_microstrip_line as coefficient_law

## The relative permittivity of the microstrip line dielectric is 4. The effective permittivity is 2.5. The wavelength of the
## signal is 10 millimeter, and the tangent of the dielectric loss angle is 1e-4. Then the loss value will be 0.345 [dB / meter].

Args = namedtuple("Args", [
    "relative_permittivity", "effective_permittivity", "wavelength", "tangent_dielectric_loss_angle"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 4
    effective_permittivity = 2.5
    wavelength = Quantity(10 * units.millimeter)
    tangent_dielectric_loss_angle = 1e-4
    return Args(relative_permittivity=relative_permittivity,
        effective_permittivity=effective_permittivity,
        wavelength=wavelength,
        tangent_dielectric_loss_angle=tangent_dielectric_loss_angle)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(
        test_args.relative_permittivity, test_args.effective_permittivity, test_args.wavelength,
        test_args.tangent_dielectric_loss_angle)
    assert_equal(result, 0.345 * (1 / units.meter))


def test_bad_permittivities(test_args: Args) -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_permittivity,
            test_args.effective_permittivity, test_args.wavelength,
            test_args.tangent_dielectric_loss_angle)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            bad_permittivity, test_args.wavelength, test_args.tangent_dielectric_loss_angle)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.effective_permittivity, bad_wavelength, test_args.tangent_dielectric_loss_angle)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.effective_permittivity, 100, test_args.tangent_dielectric_loss_angle)


def test_bad_tangent_dielectric_loss_angle(test_args: Args) -> None:
    tangent_dielectric_loss_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.effective_permittivity, test_args.wavelength, tangent_dielectric_loss_angle)
