from collections import namedtuple
from pytest import raises, fixture
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import charge_is_quantized


@fixture(name="test_args")
def test_args_fixture():
    n = 30
    Args = namedtuple("Args", "n")
    return Args(n=n)


def test_basic_law(test_args):
    result = charge_is_quantized.calculate_charge(test_args.n)
    assert_equal(result, 4.81e-18 * units.coulomb)


def test_bad_args():
    nb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        charge_is_quantized.calculate_charge(nb)
