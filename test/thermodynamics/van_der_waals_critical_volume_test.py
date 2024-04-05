from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import van_der_waals_critical_volume as critical_law

# Description
## The critical volume of Argon in the van der Waals equation of state model is V_c = 96 cm**3/mol,
## for Argon b = 0.03201 L/mol.

Args = namedtuple("Args", "b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(0.03201 * units.liter / units.mole)
    return Args(b=b)


def test_law(test_args: Args) -> None:
    vc = critical_law.calculate_critical_molar_volume(test_args.b)
    assert_equal(vc, 96 * units.centimeter**3 / units.mole)


def test_bad_second_parameter() -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_law.calculate_critical_molar_volume(bb)
    with raises(TypeError):
        critical_law.calculate_critical_molar_volume(100)
