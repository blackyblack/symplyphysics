from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.fermi_dirac_statistics import normalization

# Description
## The average numbers of particles distributed according to the Fermi-Dirac discribution on possible
## energy levels are 0.3, 0.7, 2.1, and 0.9. The total number of particles is 4.

Args = namedtuple("Args", "n1 n2 n3 n4")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n1 = 0.3
    n2 = 0.7
    n3 = 2.1
    n4 = 0.9
    return Args(n1=n1, n2=n2, n3=n3, n4=n4)


def test_law_four_levels(test_args: Args) -> None:
    result = normalization.calculate_total_fermion_count([test_args.n1, test_args.n2, test_args.n3, test_args.n4])
    assert_equal(result, 4)


def test_law_two_levels(test_args: Args) -> None:
    result = normalization.calculate_total_fermion_count((test_args.n1, test_args.n2))
    assert_equal(result, 1)


def test_bad_numbers(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        normalization.calculate_total_fermion_count([nb, test_args.n1])
    with raises(errors.UnitsError):
        normalization.calculate_total_fermion_count([test_args.n1, test_args.n2, nb])
    with raises(errors.UnitsError):
        normalization.calculate_total_fermion_count(nb)
    with raises(TypeError):
        normalization.calculate_total_fermion_count(100)
    with raises(ValueError):
        normalization.calculate_total_fermion_count((test_args.n1, test_args.n3))
