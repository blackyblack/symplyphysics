from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import energy_two_charges_from_distance_and_relative_permittivity as energy_law

# Description
## It is known that with a charge of q1 equal to 0.25 coulomb and a charge of q2 equal to 4 coulomb,
## interaction energy is 4.49e9 joule at a distance of 2 meters between charges and relative_permittivity equal to 1.
## https://matematika-club.ru/potencialnaya-ehnergiya-zaryada-q-kalkulyator-onlajn

Args = namedtuple("Args", ["relative_permittivity", "distance", "charge_1", "charge_2"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 1
    distance = Quantity(2 * units.meter)
    charge_1 = Quantity(0.25 * units.coulomb)
    charge_2 = Quantity(4 * units.coulomb)
    return Args(relative_permittivity=relative_permittivity,
        distance=distance,
        charge_1=charge_1,
        charge_2=charge_2)


def test_basic_energy(test_args: Args) -> None:
    result = energy_law.calculate_energy(test_args.relative_permittivity, test_args.distance,
        test_args.charge_1, test_args.charge_2)
    assert_equal(result, 4.49e9 * units.joule)


def test_bad_relative_permittivity(test_args: Args) -> None:
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(relative_permittivity, test_args.distance, test_args.charge_1,
            test_args.charge_2)
    with raises(TypeError):
        energy_law.calculate_energy(True, test_args.distance, test_args.charge_1,
            test_args.charge_2)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.relative_permittivity, distance, test_args.charge_1,
            test_args.charge_2)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.relative_permittivity, 100, test_args.charge_1,
            test_args.charge_2)


def test_bad_charges(test_args: Args) -> None:
    charge_1 = Quantity(1 * units.meter)
    charge_2 = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        energy_law.calculate_energy(test_args.relative_permittivity, test_args.distance, charge_1,
            charge_2)
    with raises(TypeError):
        energy_law.calculate_energy(test_args.relative_permittivity, test_args.distance, 100, 100)
