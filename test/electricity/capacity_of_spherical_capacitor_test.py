from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import capacity_of_spherical_capacitor as capacity_law

# Description
## For a spherical capacitor with a first radius of 1e-2 meter and a second radius of 3e-2 meter,
## the capacitance is 8.34e-12 farad with a dielectric constant of 5.
## https://matematika-club.ru/ehlektroemkost-kondensatora-kalkulyator-onlajn


@fixture(name="test_args")
def test_args_fixture():
    relative_permittivity = 5
    first_radius = Quantity(1e-2 * units.meter)
    second_radius = Quantity(3e-2 * units.meter)

    Args = namedtuple("Args", ["relative_permittivity", "first_radius", "second_radius"])
    return Args(relative_permittivity=relative_permittivity,
        first_radius=first_radius,
        second_radius=second_radius)


def test_basic_capacity(test_args):
    result = capacity_law.calculate_capacity(test_args.relative_permittivity,
        test_args.first_radius, test_args.second_radius)
    assert_equal(result, 8.34e-12 * units.farad)


def test_swap_radius(test_args):
    result = capacity_law.calculate_capacity(test_args.relative_permittivity,
        test_args.second_radius, test_args.first_radius)
    assert_equal(result, 8.34e-12 * units.farad)


def test_bad_relative_permittivity(test_args):
    relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_capacity(relative_permittivity, test_args.first_radius,
            test_args.second_radius)


def test_bad_radius(test_args):
    first_radius = Quantity(1 * units.coulomb)
    second_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_capacity(test_args.relative_permittivity, first_radius,
            second_radius)
    with raises(TypeError):
        capacity_law.calculate_capacity(test_args.relative_permittivity, 100, 100)
