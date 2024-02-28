from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import compliance_is_inverse_stiffness as compliance_def

# Description
## A spring has stiffness k = 1.0e3 N/m. Then its compliance is 1.0e-3 m/N.

Args = namedtuple("Args", "k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(1.0e3 * units.newton / units.meter)
    return Args(k=k)


def test_definition(test_args: Args) -> None:
    result = compliance_def.calculate_compliance(test_args.k)
    assert_equal(result, 1.0e-3 * units.meter / units.newton)


def test_bad_stiffness(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        compliance_def.calculate_compliance(kb)
    with raises(TypeError):
        compliance_def.calculate_compliance(100)
