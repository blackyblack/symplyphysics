from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.thermodynamics.dielectrics import gibbs_energy_formula as law

Args = namedtuple("Args", "f e d")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = Quantity(10.0 * units.joule / units.meter**3)
    e = Quantity(9.0 * units.volt / units.meter)
    d = Quantity(2.0 * units.coulomb / units.meter**2)
    return Args(f=f, e=e, d=d)


def test_law(test_args: Args) -> None:
    result = law.calculate_gibbs_energy_density(test_args.f, test_args.e, test_args.d)
    assert_equal(result, -8 * units.joule / units.meter**3, relative_tolerance=4e-3)


def test_bad_free_energy_density(test_args: Args) -> None:
    fb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density(fb, test_args.e, test_args.d)
    with raises(TypeError):
        law.calculate_gibbs_energy_density(100, test_args.e, test_args.d)


def test_bad_electric_field_strength(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density(test_args.f, eb, test_args.d)
    with raises(TypeError):
        law.calculate_gibbs_energy_density(test_args.f, 100, test_args.d)


def test_bad_electric_displacement(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_gibbs_energy_density(test_args.f, test_args.e, db)
    with raises(TypeError):
        law.calculate_gibbs_energy_density(test_args.f, test_args.e, 100)
