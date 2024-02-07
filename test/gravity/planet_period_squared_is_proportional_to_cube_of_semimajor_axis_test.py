from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity import (
    planet_period_squared_is_proportional_to_cube_of_semimajor_axis as law_of_periods
)

# Description
## A planet orbits around a star of mass M = 3e30 kg. The semimajor axis of its orbit is a = 3e9 km.
## Its period of rotation around the star is about 73 years.


@fixture(name="test_args")
def test_args_fixture():
    M = Quantity(3e30 * units.kilogram)
    a = Quantity(3e9 * units.kilometer)
    Args = namedtuple("Args", "M a")
    return Args(M=M, a=a)


def test_law(test_args):
    result = law_of_periods.calculate_rotation_period(test_args.M, test_args.a)
    assert_equal(result, 73.1 * units.year, tolerance=1e-2)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law_of_periods.calculate_rotation_period(mb, test_args.a)
    with raises(TypeError):
        law_of_periods.calculate_rotation_period(100, test_args.a)


def test_bad_length(test_args):
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law_of_periods.calculate_rotation_period(test_args.M, ab)
    with raises(TypeError):
        law_of_periods.calculate_rotation_period(test_args.M, 100)
