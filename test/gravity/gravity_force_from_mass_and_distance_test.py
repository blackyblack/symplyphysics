from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_law

# Description
## For example we are calculating gravity force for objects with masses 3000 and 5000 kilograms
## with 0.06 meters of distance betweem their centers of mass.
## According to matematika-club.ru gravity calculator, force should be equal to 0.27809583 Newtons.


@fixture(name="test_args")
def test_args_fixture():
    m1 = Quantity(3000 * units.kilogram)
    m2 = Quantity(5000 * units.kilogram)
    R = Quantity(0.06 * units.meter)
    Args = namedtuple("Args", ["m1", "m2", "R"])
    return Args(m1=m1, m2=m2, R=R)


def test_basic_force(test_args):
    result = gravity_law.calculate_force(test_args.m1, test_args.m2, test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).subs(units.newton, 1).evalf(7)
    assert result_force == approx(0.27809583, 0.0000001)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gravity_law.calculate_force(mb, test_args.m2, test_args.R)
    with raises(TypeError):
        gravity_law.calculate_force(100, test_args.m2, test_args.R)
    with raises(errors.UnitsError):
        gravity_law.calculate_force(test_args.m1, mb, test_args.R)
    with raises(TypeError):
        gravity_law.calculate_force(test_args.m1, 100, test_args.R)


def test_bad_distance(test_args):
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gravity_law.calculate_force(test_args.m1, test_args.m2, db)
    with raises(TypeError):
        gravity_law.calculate_force(test_args.m1, test_args.m2, 100)
