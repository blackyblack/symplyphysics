from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.conditions.electricity.maxwell_equations import derivative_of_magnetic_induction_in_time_is_rotor_of_electric_intensity as induction

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
    assert_equal(result, 1e10 / units.centimeter**3, tolerance=0.1)