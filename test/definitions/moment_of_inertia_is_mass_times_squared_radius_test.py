from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import moment_of_inertia_is_mass_times_squared_radius as moment_of_inertia_def

# Description
## Assume particle with 5kgs of mass is about to spin around axle, and a distance to this axle is 3m.
## Moment of inertia of this system should be 45kg*m**2.


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(5 * units.kilogram)
    R = Quantity(3 * units.meter)
    Args = namedtuple("Args", ["m", "R"])
    return Args(m=m, R=R)


def test_basic_moment_of_inertia(test_args):
    result = moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_ = convert_to(result, moment_of_inertia_def.definition_units_SI).evalf(2)
    assert result_ == approx(45.0, 0.01)


def test_inertia_with_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_of_inertia_def.calculate_moment_of_inertia(mb, test_args.R)
    with raises(TypeError):
        moment_of_inertia_def.calculate_moment_of_inertia(100, test_args.R)


def test_inertia_with_bad_radius(test_args):
    Rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, Rb)
    with raises(TypeError):
        moment_of_inertia_def.calculate_moment_of_inertia(test_args.m, 100)
