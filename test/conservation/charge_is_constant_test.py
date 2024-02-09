from collections import namedtuple
from pytest import raises, fixture
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import charge_is_constant

# An isolated system contained -1 C of charge in total. At any time in the future
# it will contain the same amount of charge.

Args = namedtuple("Args", "charge_before")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    charge_before = Quantity(-1 * units.coulomb)
    return Args(charge_before=charge_before)


def test_basic_law(test_args: Args) -> None:
    result = charge_is_constant.calculate_charge_after(test_args.charge_before)
    assert_equal(result, -1 * units.coulomb)


def test_bad_args() -> None:
    bad_charge_before = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        charge_is_constant.calculate_charge_after(bad_charge_before)
    with raises(TypeError):
        charge_is_constant.calculate_charge_after(100)
