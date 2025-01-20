from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, units, Quantity, errors
from symplyphysics.laws.electricity import (
    electric_displacement_is_permittivity_times_electric_field as law,)

Args = namedtuple("Args", "eps e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    eps = Quantity(81 * units.vacuum_permittivity)
    e = Quantity(1e3 * units.volt / units.meter)
    return Args(eps=eps, e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_electric_displacement(test_args.eps, test_args.e)
    assert_equal(result, 7.2e-7 * units.coulomb / units.meter**2, tolerance=4e-3)


def test_bad_permittivity(test_args: Args) -> None:
    epsb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_electric_displacement(epsb, test_args.e)
    with raises(TypeError):
        law.calculate_electric_displacement(100, test_args.e)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_electric_displacement(test_args.eps, eb)
    with raises(TypeError):
        law.calculate_electric_displacement(test_args.eps, 100)
