from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.relativistic import energy_is_mass as fundamental_formula

# 6kg of some object contains 539253107242090584 Joules of inner energy, which might be released with help of nucleal reaction.
# Calculated with "https://www.omnicalculator.com/physics/emc2"


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(6 * units.kilogram)
    Args = namedtuple("Args", ["m"])
    return Args(m=m)


def test_basic_energy(test_args):
    result = fundamental_formula.calculate_rest_energy(test_args.m)
    assert_equal(result, 539253107242090584 * units.joule)


def test_bad_mass():
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        fundamental_formula.calculate_rest_energy(mb)
    with raises(TypeError):
        fundamental_formula.calculate_rest_energy(100)
