from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics import reduced_mass_of_two_body_system as law

Args = namedtuple("Args", "m1 m2")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m1 = Quantity(3.0 * units.kilogram)
    m2 = Quantity(2.0 * units.kilogram)
    return Args(m1=m1, m2=m2)


def test_law(test_args: Args) -> None:
    result = law.calculate_reduced_mass(test_args.m1, test_args.m2)
    assert_equal(result, 1.2 * units.kilogram)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_reduced_mass(mb, test_args.m2)
    with raises(TypeError):
        law.calculate_reduced_mass(100, test_args.m2)
    with raises(errors.UnitsError):
        law.calculate_reduced_mass(test_args.m1, mb)
    with raises(TypeError):
        law.calculate_reduced_mass(test_args.m1, 100)
