from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    SI,
    convert_to,
    Quantity,
    errors
)
from symplyphysics.laws.condensed_matter import thermionic_emission_current as emission_law

# Description
## It is known that for strontium, the current density of thermionic emission at
## a temperature of 300 kelvin is equal to 0.33e-32 [A / m^2].
## https://saecanet.com/Calculation_data/000034_000170_richardson-dushman.html


@fixture(name="test_args")
def test_args_fixture():
    thermodynamic_work = Quantity(2.59 * units.electronvolt)
    temperature = Quantity(300 * units.kelvin)

    Args = namedtuple("Args", ["thermodynamic_work", "temperature"])
    return Args(
        thermodynamic_work=thermodynamic_work,
        temperature=temperature)


def test_basic_thermionic_current(test_args):
    result = emission_law.calculate_current(test_args.thermodynamic_work, test_args.temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current / units.area)
    result = convert_to(result, units.ampere / units.meter**2).evalf(5)
    assert result == approx(0.33e-32, rel=0.01)


def test_bad_thermodynamic_work(test_args):
    thermodynamic_work = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        emission_law.calculate_current(thermodynamic_work, test_args.temperature)
    with raises(TypeError):
        emission_law.calculate_current(100, test_args.temperature)


def test_bad_temperature(test_args):
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        emission_law.calculate_current(test_args.thermodynamic_work, temperature)
    with raises(TypeError):
        emission_law.calculate_current(test_args.thermodynamic_work, 100)
