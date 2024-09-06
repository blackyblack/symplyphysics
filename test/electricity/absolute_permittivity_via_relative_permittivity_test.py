from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import absolute_permittivity_via_relative_permittivity as law

Args = namedtuple("Args", "e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = 5.5
    return Args(e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_absolute_permittivity(test_args.e)
    assert_equal(result, 4.87e-11 * units.farad / units.meter)


def test_bad_number() -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_absolute_permittivity(eb)
