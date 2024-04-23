from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    isochoric_heat_capacity_via_isobaric_heat_capacity_for_ideal_gas as mayers_relation,)

# Description
## For some ideal gas the isobaric heat capacity is C_p = 5 J/K. The amount of gas substance
## is n = 0.24 mol. Its isochoric heat capacity amounts to C_v = 3 J/K.

Args = namedtuple("Args", "cp n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cp = Quantity(5 * units.joule / units.kelvin)
    n = Quantity(0.24 * units.mole)
    return Args(cp=cp, n=n)


def test_law(test_args: Args) -> None:
    result = mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, test_args.n)
    assert_equal(result, 3 * units.joule / units.kelvin, tolerance=2e-3)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(cb, test_args.n)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(100, test_args.n)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, nb)
    with raises(TypeError):
        mayers_relation.calculate_isochoric_heat_capacity(test_args.cp, 100)
