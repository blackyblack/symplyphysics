from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.thermodynamics import (
    entropy_of_independent_subsystems_is_sum_of_their_entropies as subadditivity_law,)

Args = namedtuple("Args", "s1 s2 s3")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s1 = Quantity(1 * units.joule / units.kelvin)
    s2 = Quantity(2 * units.joule / units.kelvin)
    s3 = Quantity(5 * units.joule / units.kelvin)
    return Args(s1=s1, s2=s2, s3=s3)


def test_law(test_args: Args) -> None:
    result = subadditivity_law.calculate_total_entropy((test_args.s1, test_args.s2))
    assert_equal(result, 3 * units.joule / units.kelvin)
    result = subadditivity_law.calculate_total_entropy((test_args.s1, test_args.s2, test_args.s3))
    assert_equal(result, 8 * units.joule / units.kelvin)


def test_bad_entropies(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        subadditivity_law.calculate_total_entropy(sb)
    with raises(errors.UnitsError):
        subadditivity_law.calculate_total_entropy((sb, test_args.s2))
    with raises(errors.UnitsError):
        subadditivity_law.calculate_total_entropy((test_args.s1, sb))
    with raises(TypeError):
        subadditivity_law.calculate_total_entropy(100)
    with raises(TypeError):
        subadditivity_law.calculate_total_entropy((100, test_args.s2))
    with raises(TypeError):
        subadditivity_law.calculate_total_entropy(test_args.s1)
