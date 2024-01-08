from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units,
    SI,
    convert_to,
    Quantity,
)
from symplyphysics.laws.condensed_matter import concentration_of_intrinsic_charge_carriers as concentration_law

# Description
## It is known that the concentration of intrinsic charge carriers in silicon is equal to 1e10 [1 / cm^3].
## https://www.ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html


@fixture(name="test_args")
def test_args_fixture():
    density_conductivity = Quantity(3.2e19 * (1 / units.centimeter**3))
    density_valence = Quantity(1.8e19 * (1 / units.centimeter**3))
    temperature = Quantity(300 * units.kelvin)
    band_gap = Quantity(1.12 * 1.6e-19 * units.joule)

    Args = namedtuple("Args", ["density_conductivity", "density_valence", "temperature", "band_gap"])
    return Args(
        density_conductivity=density_conductivity,
        density_valence=density_valence,
        band_gap=band_gap, temperature=temperature)


def test_basic_concentration_intrinsic(test_args):
    result = concentration_law.calculate_concentration(test_args.density_conductivity,
        test_args.density_valence, test_args.band_gap, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length**3)
    result = convert_to(result, 1 / units.centimeter**3).evalf(5)
    assert result == approx(1e10, rel=0.1)
