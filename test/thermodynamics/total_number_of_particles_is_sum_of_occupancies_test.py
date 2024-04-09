from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    total_number_of_particles_is_sum_of_occupancies as normalization_law,
)

# Description
## The average numbers of particles are 0.3, 0.7, 2.1, and 0.9. The total number of particles is 4.

Args = namedtuple("Args", "n1 n2 n3 n4")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n1 = 0.3
    n2 = 0.7
    n3 = 2.1
    n4 = 0.9
    return Args(n1=n1, n2=n2, n3=n3, n4=n4)


def test_law_four_levels(test_args: Args) -> None:
    result = normalization_law.calculate_total_particle_count([test_args.n1, test_args.n2, test_args.n3, test_args.n4])
    assert_equal(result, 4)


def test_law_two_levels(test_args: Args) -> None:
    result = normalization_law.calculate_total_particle_count((test_args.n1, test_args.n2))
    assert_equal(result, 1)


def test_bad_numbers(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        normalization_law.calculate_total_particle_count([nb, test_args.n1])
    with raises(errors.UnitsError):
        normalization_law.calculate_total_particle_count([test_args.n1, test_args.n2, nb])
    with raises(errors.UnitsError):
        normalization_law.calculate_total_particle_count(nb)
    with raises(TypeError):
        normalization_law.calculate_total_particle_count(100)
    with raises(ValueError):
        normalization_law.calculate_total_particle_count((test_args.n1, test_args.n3))
