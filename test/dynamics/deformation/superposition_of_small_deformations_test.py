from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    superposition_of_small_deformations as law,
)

Args = namedtuple("Args", "e1 e2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e1 = 0.01
    e2 = -0.005
    return Args(e1=e1, e2=e2)


def test_law(test_args: Args) -> None:
    result = law.calculate_total_strain(test_args.e1, test_args.e2)
    assert_equal(result, 0.005)


def test_bad_strain(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_total_strain(eb, test_args.e2)
    with raises(errors.UnitsError):
        law.calculate_total_strain(test_args.e1, eb)
