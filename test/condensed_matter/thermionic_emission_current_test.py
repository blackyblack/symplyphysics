from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    units,
    SI,
    convert_to,
    Quantity,
    dimensionless
)
from symplyphysics.laws.condensed_matter import thermionic_emission_current as emission_law

# Description
## It is known that for strontium, the current density of thermionic emission at
## a temperature of 300 kelvin is equal to 0.53e4 [A / m^2].
## https://saecanet.com/Calculation_data/000034_000170_richardson-dushman.html


@fixture(name="test_args")
def test_args_fixture():
    thermodynamic_work = Quantity(2.59 * 1.6e-19 * units.joule)
    electron_mass_ratio = 1
    temperature = Quantity(300 * units.kelvin)

    Args = namedtuple("Args", ["thermodynamic_work", "electron_mass_ratio", "temperature"])
    return Args(
        thermodynamic_work=thermodynamic_work,
        electron_mass_ratio=electron_mass_ratio,
        temperature=temperature)


def test_basic_thermionic_current(test_args):
    result = emission_law.calculate_current(test_args.thermodynamic_work,
        test_args.electron_mass_ratio, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current / units.length**2)
    result = convert_to(result, units.A / units.meter**2).evalf(5)
    assert result == approx(0.33e-32, rel=0.01)
