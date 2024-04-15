from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.quantities import quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law

Args = namedtuple("Args", "vm mm n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    vm = Quantity(1 * units.meter**3 / units.mole)
    mm = Quantity(40 * units.gram / units.mole)
    n = Quantity(1 * units.mole)
    return Args(vm=vm, mm=mm, n=n)


def test_law_volume(test_args: Args) -> None:
    result = molar_qty_law.calculate_extensive_quantity(test_args.vm, test_args.n)
    assert_equal(result, 1 * units.meter**3)


def test_law_mass(test_args: Args) -> None:
    result = molar_qty_law.calculate_extensive_quantity(test_args.mm, test_args.n)
    assert_equal(result, 40 * units.gram)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        molar_qty_law.calculate_extensive_quantity(test_args.vm, nb)
    with raises(TypeError):
        molar_qty_law.calculate_extensive_quantity(test_args.vm, 100)
