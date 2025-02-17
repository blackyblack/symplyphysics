from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.condensed_matter import concentration_of_intrinsic_charge_carriers as concentration_law

# Description
## It is known that the concentration of intrinsic charge carriers in silicon is equal to 1e10 [1 / cm^3].
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html

Args = namedtuple("Args", [
    "density_of_states_in_conduction_band", "density_of_states_in_valence_band", "temperature",
    "band_gap"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    density_of_states_in_conduction_band = Quantity(3.2e19 * (1 / units.centimeter**3))
    density_of_states_in_valence_band = Quantity(1.8e19 * (1 / units.centimeter**3))
    temperature = Quantity(300 * units.kelvin)
    band_gap = Quantity(1.12 * units.electronvolt)

    return Args(density_of_states_in_conduction_band=density_of_states_in_conduction_band,
        density_of_states_in_valence_band=density_of_states_in_valence_band,
        band_gap=band_gap,
        temperature=temperature)


def test_basic_charge_carriers_concentration(test_args: Args) -> None:
    result = concentration_law.calculate_concentration(
        test_args.density_of_states_in_conduction_band, test_args.density_of_states_in_valence_band,
        test_args.band_gap, test_args.temperature)
    # NOTE: high tolerance due to very high inaccuracy of available examples
    assert_equal(result, 1e10 / units.centimeter**3, relative_tolerance=0.1)


def test_bad_density_of_states_in_conduction_band(test_args: Args) -> None:
    density_of_states_in_conduction_band = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(100, test_args.density_of_states_in_valence_band,
            test_args.band_gap, test_args.temperature)


def test_bad_density_of_states_in_valence_band(test_args: Args) -> None:
    density_of_states_in_valence_band = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            density_of_states_in_valence_band, test_args.band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            100, test_args.band_gap, test_args.temperature)


def test_bad_band_gap(test_args: Args) -> None:
    band_gap = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, 100, test_args.temperature)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, 100)
