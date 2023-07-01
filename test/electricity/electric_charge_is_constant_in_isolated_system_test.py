from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import electric_charge_is_constant_in_isolated_system as charge_law


@fixture(name="test_args")
def test_args_fixture():
    Qs = Quantity(1 * units.coulomb)
    Args = namedtuple("Args", ["Qs"])
    return Args(Qs=Qs)


def test_basic_charge_conservation(test_args):
    result = charge_law.calculate_charge_after(test_args.Qs)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge)
    result_charge = convert_to(result, units.coulomb).subs(units.coulomb, 1).evalf(2)
    assert result_charge == approx(1.0, 0.001)


def test_bad_charge():
    Qb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        charge_law.calculate_charge_after(Qb)
    with raises(AttributeError):
        charge_law.calculate_charge_after(100)
