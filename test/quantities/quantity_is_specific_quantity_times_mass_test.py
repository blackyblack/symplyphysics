from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import quantity_is_specific_quantity_times_mass as specific_qty_law

Args = namedtuple("Args", "c v m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    c = Quantity(1 * units.joule / units.kilogram)
    v = Quantity(1 * units.liter / units.kilogram)
    m = Quantity(1 * units.kilogram)
    return Args(c=c, v=v, m=m)


def test_specific_heat(test_args: Args) -> None:
    result = specific_qty_law.calculate_extensive_quantity(test_args.c, test_args.m)
    assert_equal(result, 1 * units.joule)


def test_specific_volume(test_args: Args) -> None:
    result = specific_qty_law.calculate_extensive_quantity(test_args.v, test_args.m)
    assert_equal(result, 1 * units.liter)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        specific_qty_law.calculate_extensive_quantity(test_args.v, mb)
    with raises(TypeError):
        specific_qty_law.calculate_extensive_quantity(test_args.v, 100)
