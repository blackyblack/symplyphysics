from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.condensed_matter import concentration_of_intrinsic_charge_carriers as concentration_law

# Description
## It is known that the concentration of intrinsic charge carriers in silicon is equal to 1e10 [1 / cm^3].
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html


@fixture(name="test_args")
def test_args_fixture():
    density_of_states_in_conduction_band = Quantity(3.2e19 * (1 / units.centimeter**3))
    density_of_states_in_valence_band = Quantity(1.8e19 * (1 / units.centimeter**3))
    temperature = Quantity(300 * units.kelvin)
    band_gap = Quantity(1.12 * units.electronvolt)

    Args = namedtuple("Args", [
        "density_of_states_in_conduction_band", "density_of_states_in_valence_band", "temperature",
        "band_gap"
    ])
    return Args(density_of_states_in_conduction_band=density_of_states_in_conduction_band,
        density_of_states_in_valence_band=density_of_states_in_valence_band,
        band_gap=band_gap,
        temperature=temperature)


def test_basic_charge_carriers_concentration(test_args):
    result = concentration_law.calculate_concentration(
        test_args.density_of_states_in_conduction_band, test_args.density_of_states_in_valence_band,
        test_args.band_gap, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.volume)
    result = convert_to(result, 1 / units.centimeter**3).evalf(5)
    assert result == approx(1e10, rel=0.1)


def test_bad_density_of_states_in_conduction_band(test_args):
    density_of_states_in_conduction_band = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(100, test_args.density_of_states_in_valence_band,
            test_args.band_gap, test_args.temperature)


def test_bad_density_of_states_in_valence_band(test_args):
    density_of_states_in_valence_band = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            density_of_states_in_valence_band, test_args.band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            100, test_args.band_gap, test_args.temperature)


def test_bad_band_gap(test_args):
    band_gap = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, band_gap, test_args.temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, 100, test_args.temperature)


def test_bad_temperature(test_args):
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, temperature)
    with raises(TypeError):
        concentration_law.calculate_concentration(test_args.density_of_states_in_conduction_band,
            test_args.density_of_states_in_valence_band, test_args.band_gap, 100)
