from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.springs import stiffness_of_two_serial_springs as spring_law

# Description
## Two springs are connected in parallel, the first spring's constant is 50 N/m and the second's
## is 100 N/m. The total spring constant of the system of springs is 150 N/m.

Args = namedtuple("Args", "k1 k2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k1 = Quantity(50.0 * units.newton / units.meter)
    k2 = Quantity(100.0 * units.newton / units.meter)
    return Args(k1=k1, k2=k2)


def test_law(test_args: Args) -> None:
    result = spring_law.calculate_total_stiffness(test_args.k1, test_args.k2)
    assert_equal(result, 150.0 * units.newton / units.meter)


def test_bad_spring_constants(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        spring_law.calculate_total_stiffness(kb, test_args.k2)
    with raises(TypeError):
        spring_law.calculate_total_stiffness(100, test_args.k2)
    with raises(errors.UnitsError):
        spring_law.calculate_total_stiffness(test_args.k1, kb)
    with raises(TypeError):
        spring_law.calculate_total_stiffness(test_args.k1, 100)
