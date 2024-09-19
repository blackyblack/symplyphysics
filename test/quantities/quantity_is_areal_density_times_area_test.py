from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import quantity_is_areal_density_times_area as law

Args = namedtuple("Args", "m q a")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(4.4 * units.kilogram / units.meter**2)
    q = Quantity(3.5 * units.coulomb / units.meter**2)
    a = Quantity(0.1 * units.meter**2)
    return Args(m=m, q=q, a=a)


def test_mass_law(test_args: Args) -> None:
    result = law.calculate_extensive_quantity(test_args.m, test_args.a)
    assert_equal(result, 0.44 * units.kilogram)


def test_charge_law(test_args: Args) -> None:
    result = law.calculate_extensive_quantity(test_args.q, test_args.a)
    assert_equal(result, 0.35 * units.coulomb)


def test_bad_area(test_args: Args) -> None:
    ab = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_extensive_quantity(test_args.m, ab)
    with raises(TypeError):
        law.calculate_extensive_quantity(test_args.m, 100)
