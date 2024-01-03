from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (errors, units, convert_to, Quantity, SI, dimensionless)
from symplyphysics.laws.thermodynamics import efficiency_factor as efficiency_law

# Description
## Suppose we have a heat engine that receives 2 Joules of heat from the heater, and gives 1.5 Joules to the refrigerator. Then its efficiency should be (2 - 1.5) / 2 = 0.25


@fixture(name="test_args")
def test_args_fixture():
    Q_h = Quantity(2 * units.energy)
    Q_r = Quantity(1.5 * units.energy)
    Args = namedtuple("Args", ["Q_h", "Q_r"])
    return Args(Q_h=Q_h, Q_r=Q_r)


def test_basic_efficiency(test_args):
    result = efficiency_law.calculate_efficiency_factor(test_args.Q_h, test_args.Q_r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_efficiency = convert_to(result, dimensionless).evalf(6)
    assert result_efficiency == approx(0.25, 0.001)


def test_bad_efficiency(test_args):
    Q_rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        efficiency_law.calculate_efficiency_factor(test_args.Q_h, Q_rb)
        efficiency_law.calculate_efficiency_factor(Q_rb, test_args.Q_r)
    with raises(TypeError):
        efficiency_law.calculate_efficiency_factor(test_args.Q_h, 1.5)
        efficiency_law.calculate_efficiency_factor(2, test_args.Q_h)
