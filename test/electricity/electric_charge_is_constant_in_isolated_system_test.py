from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import electric_charge_is_constant_in_isolated_system as charge_law


@fixture(name="test_args")
def test_args_fixture():
    Qs = Quantity(1 * units.coulomb)
    Args = namedtuple("Args", ["Qs"])
    return Args(Qs=Qs)


def test_basic_charge_conservation(test_args):
    result = charge_law.calculate_charge_after(test_args.Qs)
    assert_equal(result, 1 * units.coulomb)


def test_bad_charge():
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        charge_law.calculate_charge_after(Qb)
    with raises(TypeError):
        charge_law.calculate_charge_after(100)
