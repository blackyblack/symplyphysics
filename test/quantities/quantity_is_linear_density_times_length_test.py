from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import quantity_is_linear_density_times_length as law

Args = namedtuple("Args", "m q l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(4.4 * units.kilogram / units.meter)
    q = Quantity(3.5 * units.coulomb / units.meter)
    l = Quantity(0.1 * units.meter)
    return Args(m=m, q=q, l=l)


def test_mass_law(test_args: Args) -> None:
    result = law.calculate_extensive_quantity(test_args.m, test_args.l)
    assert_equal(result, 0.44 * units.kilogram)


def test_charge_law(test_args: Args) -> None:
    result = law.calculate_extensive_quantity(test_args.q, test_args.l)
    assert_equal(result, 0.35 * units.coulomb)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(units.second)
    with raises(errors.UnitsError):
        law.calculate_extensive_quantity(test_args.m, lb)
    with raises(TypeError):
        law.calculate_extensive_quantity(test_args.m, 100)
