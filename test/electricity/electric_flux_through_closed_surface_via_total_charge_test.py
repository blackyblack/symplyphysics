from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import electric_flux_through_closed_surface_via_total_charge as law

Args = namedtuple("Args", "q")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q = Quantity(-5e-3 * units.coulomb)
    return Args(q=q)


def test_law(test_args: Args) -> None:
    result = law.calculate_total_electric_flux(test_args.q)
    assert_equal(result, -5.65e8 * units.volt * units.meter)


def test_bad_charge() -> None:
    qb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_total_electric_flux(qb)
    with raises(TypeError):
        law.calculate_total_electric_flux(100)
