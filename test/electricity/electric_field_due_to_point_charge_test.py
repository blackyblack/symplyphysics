from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    prefixes,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as electric_field

# The electric field set up by a point charge of 1 nC at a distance of 3 m should be approximately 0.999 N/C


@fixture(name="test_args")
def test_args_fixture():
    q = Quantity(1 * prefixes.nano * units.coulomb)
    r = Quantity(3 * units.meter)
    Args = namedtuple("Args", "q r")
    return Args(q=q, r=r)


def test_basic_law(test_args):
    result = electric_field.calculate_electric_field(test_args.q, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force / units.charge)
    result_field = convert_to(result, units.newton / units.coulomb).evalf(3)
    assert_approx(result_field, 0.999)


def test_bad_charge(test_args):
    qb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(qb, test_args.r)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.r)


def test_bad_distance(test_args):
    rb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.q, rb)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.q, 100)
