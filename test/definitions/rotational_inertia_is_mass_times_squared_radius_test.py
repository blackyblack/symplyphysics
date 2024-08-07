from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import rotational_inertia_is_mass_times_squared_radius as moment_of_inertia_def

# Description
## Assume particle with 5kgs of mass is about to spin around axle, and a distance to this axle is 3m.
## Moment of inertia of this system should be 45kg*m**2.

Args = namedtuple("Args", ["m", "R"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = Quantity(5 * units.kilogram)
    R = Quantity(3 * units.meter)
    return Args(m=m, R=R)


def test_basic_moment_of_inertia(test_args: Args) -> None:
    result = moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, test_args.R)
    assert_equal(result, 45 * units.kilogram * units.meter**2)


def test_inertia_with_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_of_inertia_def.calculate_moment_of_inertia(mb, test_args.R)
    with raises(TypeError):
        moment_of_inertia_def.calculate_moment_of_inertia(100, test_args.R)


def test_inertia_with_bad_radius(test_args: Args) -> None:
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, Rb)
    with raises(TypeError):
        moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, 100)
