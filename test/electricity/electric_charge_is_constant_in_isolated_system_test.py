from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.electricity import electric_charge_is_constant_in_isolated_system as charge_law

@fixture
def test_args():
    Q_before = Quantity(units.charge, 1 * units.coulomb)
    Args = namedtuple("Args", ["Q_before"])
    return Args(Q_before = Q_before)

def test_basic_charge_conservation(test_args):
    result = charge_law.calculate_charge_after(test_args.Q_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge)
    result_charge = convert_to(result, units.coulomb).subs(units.coulomb, 1).evalf(2)
    assert result_charge == approx(1.0, 0.001)

def test_bad_charge():
    Qb = Quantity(units.length)
    with raises(errors.UnitsError):
        charge_law.calculate_charge_after(Qb)
    with raises(TypeError):
        charge_law.calculate_charge_after(100)
