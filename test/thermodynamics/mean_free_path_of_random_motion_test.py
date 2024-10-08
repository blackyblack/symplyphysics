from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import mean_free_path_of_random_motion as mean_free_path_law

# Description
## The estimate of the mean free path of N = 10.000 particles of diameter d = 3 µm traveling randomly
## in a volume V = 1 litre is lambda = 2.5 km.

Args = namedtuple("Args", "d n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    d = Quantity(3 * prefixes.micro * units.meter)
    n_ = 10_000
    v = Quantity(1 * units.liter)
    n = Quantity(n_ / v)
    return Args(d=d, n=n)


def test_law(test_args: Args) -> None:
    result = mean_free_path_law.calculate_mean_free_path(test_args.d, test_args.n)
    assert_equal(result, 2.5 * units.kilometer)


def test_bad_diameter(test_args: Args) -> None:
    db = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mean_free_path_law.calculate_mean_free_path(db, test_args.n)
    with raises(TypeError):
        mean_free_path_law.calculate_mean_free_path(100, test_args.n)


def test_bad_number_density(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        mean_free_path_law.calculate_mean_free_path(test_args.d, nb)
    with raises(TypeError):
        mean_free_path_law.calculate_mean_free_path(test_args.d, 100)
