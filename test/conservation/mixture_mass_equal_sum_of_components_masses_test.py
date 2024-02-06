from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.conservation import mixture_mass_equal_sum_of_components_masses as sum_of_masses_law


@fixture(name="test_args")
def test_args_fixture():
    m1 = Quantity(10 * units.kilograms)
    m2 = Quantity(20 * units.kilograms)
    Args = namedtuple("Args", ["m1", "m2"])
    return Args(m1=m1, m2=m2)


def test_basic_superposition(test_args):
    result = sum_of_masses_law.calculate_mass_of_mixture([test_args.m1, test_args.m2])
    assert_equal(result, 30 * units.kilograms)


def test_three_masses_array(test_args):
    m3 = Quantity(5 * units.kilograms)
    result = sum_of_masses_law.calculate_mass_of_mixture([test_args.m1, test_args.m2, m3])
    assert_equal(result, 35 * units.kilograms)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        sum_of_masses_law.calculate_mass_of_mixture([mb, test_args.m2])
    with raises(TypeError):
        sum_of_masses_law.calculate_mass_of_mixture([100, test_args.m2])
    with raises(errors.UnitsError):
        sum_of_masses_law.calculate_mass_of_mixture([test_args.m1, mb])
    with raises(TypeError):
        sum_of_masses_law.calculate_mass_of_mixture([test_args.m1, 100])
    with raises(errors.UnitsError):
        sum_of_masses_law.calculate_mass_of_mixture([mb, mb])
    with raises(TypeError):
        sum_of_masses_law.calculate_mass_of_mixture([100, 100])
    with raises(TypeError):
        sum_of_masses_law.calculate_mass_of_mixture(test_args.m1)
