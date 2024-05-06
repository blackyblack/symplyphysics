from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import (
    statistical_weight_of_macrostate as weight_law,
)

Args = namedtuple("Args", "n1 n2 n3")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n1 = 3
    n2 = 10
    n3 = 8
    return Args(n1=n1, n2=n2, n3=n3)


def test_law_two_states(test_args: Args) -> None:
    result = weight_law.calculate_statistical_weight((test_args.n1, test_args.n2))
    assert_equal(result, 286)


def test_law_three_states(test_args: Args) -> None:
    result = weight_law.calculate_statistical_weight((test_args.n1, test_args.n2, test_args.n3))
    assert_equal(result, 58198140)


def test_bad_numbers(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        weight_law.calculate_statistical_weight((nb, test_args.n2))
    with raises(errors.UnitsError):
        weight_law.calculate_statistical_weight((test_args.n1, nb))
    with raises(errors.UnitsError):
        weight_law.calculate_statistical_weight(nb)
    with raises(TypeError):
        weight_law.calculate_statistical_weight(100)
