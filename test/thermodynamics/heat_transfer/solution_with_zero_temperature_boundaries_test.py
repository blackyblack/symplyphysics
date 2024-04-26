from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.heat_transfer import solution_with_zero_temperature_boundaries as solution

# Description
## Thermal diffusivity of aluminium is chi = 97 mm**2/s. The length of the rod is 0.4 m.
## At t = 1 s and x = 0.3 m. The value of the solution of the heat equation is T = 67 K
## for n = 3 and B = 100 K.

Args = namedtuple("Args", "b k n l x t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    b = Quantity(100 * units.kelvin)
    k = Quantity(97 * units.millimeter**2 / units.second)
    n = 3
    l = Quantity(0.4 * units.meter)
    x = Quantity(0.3 * units.meter)
    t = Quantity(1 * units.second)
    return Args(b=b, k=k, n=n, l=l, x=x, t=t)


def test_law(test_args: Args) -> None:
    result = solution.calculate_temperature(test_args.b, test_args.k, test_args.n, test_args.l, test_args.x, test_args.t)
    assert_equal(result, 67 * units.kelvin)


def test_bad_temperature_coefficient(test_args: Args) -> None:
    bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution.calculate_temperature(bb, test_args.k, test_args.n, test_args.l, test_args.x, test_args.t)
    with raises(TypeError):
        solution.calculate_temperature(100, test_args.k, test_args.n, test_args.l, test_args.x, test_args.t)


def test_bad_diffusivity(test_args: Args) -> None:
    kb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution.calculate_temperature(test_args.b, kb, test_args.n, test_args.l, test_args.x, test_args.t)
    with raises(TypeError):
        solution.calculate_temperature(test_args.b, 100, test_args.n, test_args.l, test_args.x, test_args.t)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution.calculate_temperature(test_args.b, test_args.k, nb, test_args.l, test_args.x, test_args.t)


def test_bad_position(test_args: Args) -> None:
    xb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, xb, test_args.x, test_args.t)
    with raises(TypeError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, 100, test_args.x, test_args.t)
    with raises(errors.UnitsError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, test_args.l, xb, test_args.t)
    with raises(TypeError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, test_args.l, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, test_args.l, test_args.x, tb)
    with raises(TypeError):
        solution.calculate_temperature(test_args.b, test_args.k, test_args.n, test_args.l, test_args.x, 100)
