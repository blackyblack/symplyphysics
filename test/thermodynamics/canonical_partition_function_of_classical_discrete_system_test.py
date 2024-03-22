from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    canonical_partition_function_of_classical_discrete_system as partition_function_law,
)

# Description
## A canonical ensemble is described as having Boltzmann factors of 0.1, 0.66, 0.3 and 0.9.
## The partition function of the ensemble amounts to 1.96.

Args = namedtuple("Args", "f1 f2 f3 f4")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f1 = 0.1
    f2 = 0.66
    f3 = 0.3
    f4 = 0.9
    return Args(f1=f1, f2=f2, f3=f3, f4=f4)


def test_law_four_factors(test_args: Args) -> None:
    result = partition_function_law.calculate_partition_function([test_args.f1, test_args.f2, test_args.f3, test_args.f4])
    assert_equal(result, 1.96)


def test_law_three_factors(test_args: Args) -> None:
    result = partition_function_law.calculate_partition_function([test_args.f1, test_args.f2, test_args.f3])
    assert_equal(result, 1.06)


def test_bad_boltzmann_factors(test_args: Args) -> None:
    fb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        partition_function_law.calculate_partition_function([fb, test_args.f1])
    with raises(errors.UnitsError):
        partition_function_law.calculate_partition_function([test_args.f1, fb])
    with raises(errors.UnitsError):
        partition_function_law.calculate_partition_function(fb)
    with raises(TypeError):
        partition_function_law.calculate_partition_function(test_args.f1)
