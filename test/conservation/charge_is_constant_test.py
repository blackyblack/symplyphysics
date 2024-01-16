from collections import namedtuple
from pytest import raises, approx, fixture
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.conservation import charge_is_constant

# An isolated system contained -1 C of charge in total. At any time in the future 
# it will contain the same amount of charge.


@fixture(name="test_args")
def test_args_fixture():
    charge_before = Quantity(-1 * units.coulomb)
    Args = namedtuple("Args", "charge_before")
    return Args(charge_before=charge_before)


def test_basic_law(test_args):
    result = charge_is_constant.calculate_charge_after(test_args.charge_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge)
    result_charge = convert_to(result, units.coulomb).evalf(3)
    assert result_charge == approx(-1, 1e-3)


def test_bad_args(test_args):
    bad_charge_before = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        charge_is_constant.calculate_charge_after(bad_charge_before)
    with raises(TypeError):
        charge_is_constant.calculate_charge_after(100)
