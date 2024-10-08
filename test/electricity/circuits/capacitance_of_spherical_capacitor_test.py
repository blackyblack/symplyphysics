from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits import capacitance_of_spherical_capacitor as capacity_law

# Description
## For a spherical capacitor with a first radius of 1e-2 meter and a second radius of 3e-2 meter,
## the capacitance is 8.34e-12 farad with a dielectric constant of 5.
## https://matematika-club.ru/ehlektroemkost-kondensatora-kalkulyator-onlajn

Args = namedtuple("Args", ["absolute_permittivity", "first_radius", "second_radius"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    _relative_permittivity = 5
    absolute_permittivity = Quantity(units.vacuum_permittivity * _relative_permittivity)
    first_radius = Quantity(1e-2 * units.meter)
    second_radius = Quantity(3e-2 * units.meter)
    return Args(absolute_permittivity=absolute_permittivity,
        first_radius=first_radius,
        second_radius=second_radius)


def test_basic_capacity(test_args: Args) -> None:
    result = capacity_law.calculate_capacity(test_args.absolute_permittivity,
        test_args.first_radius, test_args.second_radius)
    assert_equal(result, 8.34e-12 * units.farad)


def test_swap_radius(test_args: Args) -> None:
    result = capacity_law.calculate_capacity(test_args.absolute_permittivity,
        test_args.second_radius, test_args.first_radius)
    assert_equal(result, 8.34e-12 * units.farad)


def test_bad_relative_permittivity(test_args: Args) -> None:
    absolute_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_capacity(absolute_permittivity, test_args.first_radius,
            test_args.second_radius)


def test_bad_radius(test_args: Args) -> None:
    first_radius = Quantity(1 * units.coulomb)
    second_radius = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        capacity_law.calculate_capacity(test_args.absolute_permittivity, first_radius,
            second_radius)
    with raises(TypeError):
        capacity_law.calculate_capacity(test_args.absolute_permittivity, 100, 100)
