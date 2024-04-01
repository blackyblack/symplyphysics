from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_via_amount_of_heat_and_work_done as first_law,)

# Description
## The amount of heat supplied to the system is 1 J and the work done by the system is 2 J.
## The change in internal energy of the system is -1 J.

Args = namedtuple("Args", "dq dw")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dq = Quantity(1 * units.joule)
    dw = Quantity(2 * units.joule)
    return Args(dq=dq, dw=dw)


def test_law(test_args: Args) -> None:
    result = first_law.calculate_internal_energy_change(test_args.dq, test_args.dw)
    assert_equal(result, -1 * units.joule)


def test_bad_energies(test_args: Args) -> None:
    eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        first_law.calculate_internal_energy_change(eb, test_args.dw)
    with raises(TypeError):
        first_law.calculate_internal_energy_change(100, test_args.dw)
    with raises(errors.UnitsError):
        first_law.calculate_internal_energy_change(test_args.dq, eb)
    with raises(TypeError):
        first_law.calculate_internal_energy_change(test_args.dq, 100)
