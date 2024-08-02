from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    bulk_modulus_via_young_modulus_and_poisson_ratio as law,
)

# Averaged values for stainless steel are being used here.
# Taken from [source](https://www.azom.com/properties.aspx?ArticleID=965)

Args = namedtuple("Args", "e n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(196 * prefixes.giga * units.pascal)
    n = 0.27
    return Args(e=e, n=n)


def test_law(test_args: Args) -> None:
    result = law.calculate_bulk_modulus(test_args.e, test_args.n)
    assert_equal(result, 142 * prefixes.giga * units.pascal)


def test_bad_pressure(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_bulk_modulus(eb, test_args.n)
    with raises(TypeError):
        law.calculate_bulk_modulus(100, test_args.n)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_bulk_modulus(test_args.e, nb)
