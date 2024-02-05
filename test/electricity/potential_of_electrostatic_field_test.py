from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_approx, units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import potential_of_electrostatic_field as potential_law

# Description
## It is known that with a potential energy of 10 joules and a charge of 5.5 coulomb,
## the electrostatic potential is 1.818 volts.
## https://matematika-club.ru/potencial-ehlektrostaticheskogo-polya-onlajn-kalkulyator


@fixture(name="test_args")
def test_args_fixture():
    potential_energy = Quantity(10 * units.joule)
    charge = Quantity(5.5 * units.coulomb)

    Args = namedtuple("Args", ["potential_energy", "charge"])
    return Args(potential_energy=potential_energy, charge=charge)


def test_basic_electrostatic_potential(test_args):
    result = potential_law.calculate_potential(test_args.potential_energy, test_args.charge)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage)
    result = convert_to(result, units.volt).evalf(5)
    assert_approx(result, 1.818)


def test_bad_potential_energy(test_args):
    potential_energy = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        potential_law.calculate_potential(potential_energy, test_args.charge)
    with raises(TypeError):
        potential_law.calculate_potential(100, test_args.charge)


def test_bad_charge(test_args):
    charge = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        potential_law.calculate_potential(test_args.potential_energy, charge)
    with raises(TypeError):
        potential_law.calculate_potential(test_args.potential_energy, 100)
