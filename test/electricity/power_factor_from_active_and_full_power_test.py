from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI, dimensionless)
from symplyphysics.laws.electricity import power_factor_from_active_and_full_power as power_factor_law

# Description
## Assert we have a device wich consumes 10 Watt of power and makes 3 Watt of work. Power factor of this consumer should be 3/10.

@fixture(name="test_args")
def test_args_fixture():
    P = Quantity(3 * units.watt)
    S = Quantity(10 * units.watt)
    Args = namedtuple("Args", ["P", "S"])
    return Args(P=P, S=S)


def test_basic_power_factor(test_args):
    result = power_factor_law.calculate_power_factor(test_args.P, test_args.S)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_current = convert_to(result, dimensionless).evalf(6)
    assert result_current == approx(0.3, 0.001)


def test_bad_power(test_args):
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        power_factor_law.calculate_power_factor(pb, test_args.S)
        power_factor_law.calculate_power_factor(test_args.P, pb)
    with raises(TypeError):
        power_factor_law.calculate_power_factor(100, test_args.S)
        power_factor_law.calculate_power_factor(test_args.P, 100)
