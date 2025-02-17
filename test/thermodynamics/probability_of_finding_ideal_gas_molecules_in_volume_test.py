from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    probability_of_finding_ideal_gas_molecules_in_volume as probability_law,)

Args = namedtuple("Args", "v0 v n")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1 * units.liter)
    v = Quantity(0.9 * units.liter)
    n = 1000
    return Args(v0, v, n)


def test_law(test_args: Args) -> None:
    result = probability_law.calculate_probability(test_args.v0, test_args.v, test_args.n)
    assert_equal(result, 1.75e-46, relative_tolerance=2e-3)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        probability_law.calculate_probability(vb, test_args.v, test_args.n)
    with raises(TypeError):
        probability_law.calculate_probability(100, test_args.v, test_args.n)
    with raises(errors.UnitsError):
        probability_law.calculate_probability(test_args.v0, vb, test_args.n)
    with raises(TypeError):
        probability_law.calculate_probability(test_args.v0, 100, test_args.n)

    vb = Quantity(2 * test_args.v0)
    with raises(ValueError):
        probability_law.calculate_probability(test_args.v0, vb, test_args.n)


def test_bad_number(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        probability_law.calculate_probability(test_args.v0, test_args.v, nb)
