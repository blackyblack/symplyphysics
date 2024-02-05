from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import mass_fraction as mass_fraction_law


@fixture(name="test_args")
def test_args_fixture():
    m_i = Quantity(1 * units.kilogram)
    m = Quantity(3 * units.kilogram)
    Args = namedtuple("Args", ["m_i", "m"])
    return Args(m_i=m_i, m=m)


def test_basic_mass_fraction(test_args):
    result = mass_fraction_law.calculate_mass_fraction(test_args.m_i, test_args.m)
    assert_approx(result, 0.3333)


def test_bad_mass_fraction(test_args):
    with raises(AttributeError):
        mass_fraction_law.calculate_mass_fraction(test_args.m, test_args.m_i)
    with raises(AttributeError):
        mass_fraction_law.calculate_mass_fraction((-1) * test_args.m_i, test_args.m)


def test_bad_mass_of_component(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_fraction_law.calculate_mass_fraction(mb, test_args.m)
    with raises(TypeError):
        mass_fraction_law.calculate_mass_fraction(100, test_args.m)


def test_bad_mass_of_mixture(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_fraction_law.calculate_mass_fraction(test_args.m_i, mb)
    with raises(TypeError):
        mass_fraction_law.calculate_mass_fraction(test_args.m_i, 100)
