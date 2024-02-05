from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    convert_to,
)
from symplyphysics.laws.relativistic import relativistic_mass

# From https://rechneronline.de/spectrum/relativistic-mass-growth.php


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(6 * units.kilogram)
    v = Quantity(20_000_000 * (units.meter / units.second))
    Args = namedtuple("Args", ["m", "v"])
    return Args(m=m, v=v)


def test_basic_mass(test_args):
    result = relativistic_mass.calculate_relativistic_mass(test_args.m, test_args.v)
    result_mass = convert_to(result, units.kilogram).evalf(4)
    assert_approx(result_mass, 6.01339)


def test_basic_zero_velocity(test_args):
    velocity = Quantity(0 * units.meter / units.second)
    result = relativistic_mass.calculate_relativistic_mass(test_args.m, velocity)
    result_mass = convert_to(result, units.kilogram).evalf(4)
    assert_approx(result_mass, 6)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_mass.calculate_relativistic_mass(mb, test_args.v)
    with raises(TypeError):
        relativistic_mass.calculate_relativistic_mass(100, test_args.v)


def test_bad_velocity(test_args):
    mv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        relativistic_mass.calculate_relativistic_mass(test_args.m, mv)
    with raises(TypeError):
        relativistic_mass.calculate_relativistic_mass(test_args.m, 100)
