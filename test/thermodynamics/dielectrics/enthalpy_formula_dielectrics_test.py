from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.thermodynamics.dielectrics import enthalpy_formula as law

Args = namedtuple("Args", "u e d")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(20.0 * units.joule / units.meter**3)
    e = Quantity(10.0 * units.volt / units.meter)
    d = Quantity(3.0 * units.coulomb / units.meter**2)
    return Args(u=u, e=e, d=d)


def test_law(test_args: Args) -> None:
    result = law.calculate_enthalpy_density(test_args.u, test_args.e, test_args.d)
    assert_equal(result, -10 * units.joule / units.meter**3)


def test_bad_internal_energy_density(test_args: Args) -> None:
    ub = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density(ub, test_args.e, test_args.d)
    with raises(TypeError):
        law.calculate_enthalpy_density(100, test_args.e, test_args.d)


def test_bad_electric_field_strength(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density(test_args.u, eb, test_args.d)
    with raises(TypeError):
        law.calculate_enthalpy_density(test_args.u, 100, test_args.d)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_enthalpy_density(test_args.u, test_args.e, db)
    with raises(TypeError):
        law.calculate_enthalpy_density(test_args.u, test_args.e, 100)
