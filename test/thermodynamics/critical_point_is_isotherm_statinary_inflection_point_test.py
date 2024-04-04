from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    critical_point_is_isotherm_stationary_inflection_point as critical_point_law,
)

# Description
## The critical point for 1 mole of Argon (a = 1.355 bar*(L/mol)**2, b = 0.03201 L/mol)
## is V_c = 96 cm**3, p_c = 48.3 atm, T_c = 151 K.

Args = namedtuple("Args", "a b n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    a = Quantity(1.355 * units.bar * (units.liter / units.mole)**2)
    b = Quantity(0.03201 * units.liter / units.mole)
    n = Quantity(1 * units.mole)
    return Args(a=a, b=b, n=n)


def test_law(test_args: Args) -> None:
    vc, pc, tc = critical_point_law.calculate_critical_point(test_args.a, test_args.b, test_args.n)
    assert_equal(vc, 96 * units.centimeter**3)
    assert_equal(pc, 48.3 * units.atmosphere)
    assert_equal(tc, 151 * units.kelvin)


def test_bad_first_parameter(test_args: Args) -> None:
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_point_law.calculate_critical_point(ab, test_args.b, test_args.n)
    with raises(TypeError):
        critical_point_law.calculate_critical_point(100, test_args.b, test_args.n)


def test_bad_second_parameter(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_point_law.calculate_critical_point(test_args.a, bb, test_args.n)
    with raises(TypeError):
        critical_point_law.calculate_critical_point(test_args.a, 100, test_args.n)


def test_bad_amount_of_substance(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        critical_point_law.calculate_critical_point(test_args.a, test_args.b, nb)
    with raises(TypeError):
        critical_point_law.calculate_critical_point(test_args.a, test_args.b, 100)
