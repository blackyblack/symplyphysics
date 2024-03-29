from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (errors, units, Quantity, assert_equal,)

from symplyphysics.laws.electricity.waveguides import phase_velocity_of_wave_in_rectangular_waveguide as velocity_law

## The critical wavelength is 17.9 millimeters.The wavelength is 10 millimeters.
## The relative permittivity of the dielectric is 2.2, the relative permeability of the dielectric is 1.
## Then the phase velocity in the waveguide will be 243695071.57 meter per second.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["relative_permittivity", "relative_permeability", "wavelength", "critical_wavelength"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    relative_permeability = 1
    wavelength = Quantity(10 * units.millimeter)
    critical_wavelength = Quantity(17.9 * units.millimeter)
    return Args(relative_permittivity=relative_permittivity, relative_permeability=relative_permeability, wavelength=wavelength, critical_wavelength=critical_wavelength)


def test_basic_phase_velocity(test_args: Args) -> None:
    result = velocity_law.calculate_phase_velocity(test_args.relative_permittivity, test_args.relative_permeability, test_args.wavelength,
        test_args.critical_wavelength)
    assert_equal(result, 243695071.57 * (units.meter / units.second))


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_phase_velocity(relative_permittivity, test_args.relative_permeability, test_args.wavelength, test_args.critical_wavelength)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_phase_velocity(test_args.relative_permittivity, relative_permeability, test_args.wavelength, test_args.critical_wavelength)


def test_bad_wavelength(test_args: Args) -> None:
    bad_wavelength = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_phase_velocity(test_args.relative_permittivity, test_args.relative_permeability, bad_wavelength, test_args.critical_wavelength)
    with raises(TypeError):
        velocity_law.calculate_phase_velocity(test_args.relative_permittivity, test_args.relative_permeability, 100, test_args.critical_wavelength)
    with raises(errors.UnitsError):
        velocity_law.calculate_phase_velocity(test_args.relative_permittivity, test_args.relative_permeability, test_args.wavelength, bad_wavelength)
    with raises(TypeError):
        velocity_law.calculate_phase_velocity(test_args.relative_permittivity, test_args.relative_permeability, test_args.wavelength, 100)
