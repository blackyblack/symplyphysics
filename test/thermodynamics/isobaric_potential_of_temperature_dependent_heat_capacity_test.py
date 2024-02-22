from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.thermodynamics import isobaric_potential_of_temperature_dependent_heat_capacity as potential_law

# Description
## The decomposition of barium carbonate by the reaction: BaCO3 = BaO + CO2.
## At a temperature of 298 kelvin, the thermal effect of the reaction is 552866 joule per mole, the entropy change is 155 [joule / (mole * kelvin)],
## The isobaric heat capacity is -8.79 [joule / (mole * kelvin)].
## Then the isobaric potential is 506676 joule per mole.
## https://studfile.net/preview/5797036/page:4/

Args = namedtuple("Args", ["thermal_effect", "entropy", "temperature", "heat_capacity",
                           "coefficient_capacity_1", "coefficient_capacity_2", "coefficient_capacity_3"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    thermal_effect = Quantity(552866 * units.joule / units.mole)
    entropy = Quantity(155 * (units.joule / units.mole / units.kelvin))
    temperature = Quantity(298 * units.kelvin)
    heat_capacity = Quantity(-8.79 * (units.joule / units.mole / units.kelvin))
    coefficient_capacity_1 = Quantity(1 * (units.joule / units.mole / units.kelvin))
    coefficient_capacity_2 =Quantity(2 * (units.joule / units.mole / units.kelvin**2))
    coefficient_capacity_3 = Quantity(3 * (units.joule * units.kelvin / units.mole))

    return Args(thermal_effect=thermal_effect,
        entropy=entropy,
        temperature=temperature,
        heat_capacity=heat_capacity,
        coefficient_capacity_1=coefficient_capacity_1,
        coefficient_capacity_2=coefficient_capacity_2,
        coefficient_capacity_3=coefficient_capacity_3,)


def test_basic_isobaric_potential(test_args: Args) -> None:
    result = potential_law.calculate_isobaric_potential(test_args.thermal_effect,
        test_args.entropy, test_args.temperature, test_args.heat_capacity,
        test_args.coefficient_capacity_1,
        test_args.coefficient_capacity_2,
        test_args.coefficient_capacity_3,)
    assert_equal(result, 506676 * units.joule / units.mole)


def test_bad_thermal_effect(test_args: Args) -> None:
    thermal_effect = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(thermal_effect, test_args.entropy,
            test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(100, test_args.entropy,
            test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)


def test_bad_entropy(test_args: Args) -> None:
    entropy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect, entropy,
            test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect, 100,
            test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)


def test_bad_temperature(test_args: Args) -> None:
    temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, 100, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)


def test_bad_heat_capacity(test_args: Args) -> None:
    heat_capacity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, 100, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)


def test_bad_coefficients_capacity(test_args: Args) -> None:
    coefficient_capacity_1 = Quantity(1 * units.coulomb)
    coefficient_capacity_2 = Quantity(1 * units.coulomb)
    coefficient_capacity_3 = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, coefficient_capacity_1,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, 100,
            test_args.coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            coefficient_capacity_2, test_args.coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            100, test_args.coefficient_capacity_3,)
    with raises(errors.UnitsError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, coefficient_capacity_3,)
    with raises(TypeError):
        potential_law.calculate_isobaric_potential(test_args.thermal_effect,
            test_args.entropy, test_args.temperature, test_args.heat_capacity, test_args.coefficient_capacity_1,
            test_args.coefficient_capacity_2, 100,)
