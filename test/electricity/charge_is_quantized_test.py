from collections import namedtuple
from pytest import raises, approx, fixture
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import charge_is_quantized


@fixture(name="test_args")
def test_args_fixture():
    n = 30
    Args = namedtuple("Args", "n")
    return Args(n=n)


def test_basic_law(test_args):
    result = charge_is_quantized.calculate_charge(test_args.n)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge)
    result_charge = convert_to(result, units.coulomb).evalf(3)
    assert result_charge == approx(4.81e-18, 1e-18)


def test_bad_args():
    nb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        charge_is_quantized.calculate_charge(nb)
