from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors, quantities
from symplyphysics.laws.electricity import electric_flux_through_closed_surface_in_linear_material as law

Args = namedtuple("Args", "q e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q = Quantity(1e-5 * units.coulomb)
    e = Quantity(81 * quantities.vacuum_permittivity)
    return Args(q=q, e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_total_electric_flux(test_args.q, test_args.e)
    assert_equal(result, 14e3 * units.volt * units.meter, relative_tolerance=5e-3)


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_total_electric_flux(qb, test_args.e)
    with raises(TypeError):
        law.calculate_total_electric_flux(100, test_args.e)


def test_bad_absolute_permittivity(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_total_electric_flux(test_args.q, eb)
    with raises(TypeError):
        law.calculate_total_electric_flux(test_args.q, 100)
