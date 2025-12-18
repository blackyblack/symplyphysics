from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.special_relativity.relativistic_dynamics.obsolete_concepts import relativistic_mass_via_rest_mass_and_speed as law

# From https://rechneronline.de/spectrum/relativistic-mass-growth.php

Args = namedtuple("Args", ["m", "v"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(6 * units.kilogram)
    v = Quantity(20_000_000 * (units.meter / units.second))
    return Args(m=m, v=v)


def test_basic_mass(test_args: Args) -> None:
    result = law.calculate_relativistic_mass(test_args.m, test_args.v)
    assert_equal(result, 6.01339 * units.kilogram)


def test_basic_zero_velocity(test_args: Args) -> None:
    velocity = Quantity(0 * units.meter / units.second)
    result = law.calculate_relativistic_mass(test_args.m, velocity)
    assert_equal(result, 6 * units.kilogram)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relativistic_mass(mb, test_args.v)
    with raises(TypeError):
        law.calculate_relativistic_mass(100, test_args.v)


def test_bad_velocity(test_args: Args) -> None:
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_relativistic_mass(test_args.m, mv)
    with raises(TypeError):
        law.calculate_relativistic_mass(test_args.m, 100)
