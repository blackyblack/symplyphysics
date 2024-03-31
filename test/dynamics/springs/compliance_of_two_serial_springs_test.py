from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.springs import compliance_of_two_serial_springs as compliance_law

# Description
## Two springs are connected end-to-end, the first spring's compliance is 1.0e-2 m/N, the second's
## is 1.0e-3 m/N. The total compliance of the system of springs is 1.1e-2 m/N.

Args = namedtuple("Args", "c1 c2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c1 = Quantity(1.0e-2 * units.meter / units.newton)
    c2 = Quantity(1.0e-3 * units.meter / units.newton)
    return Args(c1=c1, c2=c2)


def test_law(test_args: Args) -> None:
    result = compliance_law.calculate_compliance(test_args.c1, test_args.c2)
    assert_equal(result, 1.1e-2 * units.meter / units.newton)


def test_bad_compliances(test_args: Args) -> None:
    cb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        compliance_law.calculate_compliance(cb, test_args.c2)
    with raises(TypeError):
        compliance_law.calculate_compliance(100, test_args.c2)
    with raises(errors.UnitsError):
        compliance_law.calculate_compliance(test_args.c1, cb)
    with raises(TypeError):
        compliance_law.calculate_compliance(test_args.c1, 100)
