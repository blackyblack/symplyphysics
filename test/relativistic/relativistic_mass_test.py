from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    convert_to,
)
from symplyphysics.laws.relativistic import relativistic_mass


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(6 * units.kilogram)
    v = Quantity(2_000_000 * (units.meter / units.second))
    Args = namedtuple("Args", ["m", "v"])
    return Args(m=m, v=v)


def test_basic_mass(test_args):
    result = relativistic_mass.calculate_relativistic_mass(
        test_args.m, test_args.v)
    result_mass = convert_to(result, units.kilogram).evalf(4)
    assert result_mass == approx(6.04, 0.01)


def test_basic_zero_velocity(test_args):
    velocity = Quantity(0 * units.meter / units.second)
    result = relativistic_mass.calculate_relativistic_mass(
        test_args.m, velocity)
    result_mass = convert_to(result, units.kilogram).evalf(4)
    assert result_mass == approx(6, 0.01)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    # negative_mass = Quantity(-1 * units.kilogram)
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
