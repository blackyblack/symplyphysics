from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as electric_field

# The electric field set up by a point charge of 1 nC at a distance of 3 m should be approximately 0.999 N/C

Args = namedtuple("Args", "q r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q = Quantity(1 * prefixes.nano * units.coulomb)
    r = Quantity(3 * units.meter)
    return Args(q=q, r=r)


def test_basic_law(test_args: Args) -> None:
    result = electric_field.calculate_electric_field(test_args.q, test_args.r)
    assert_equal(result, 0.999 * units.newton / units.coulomb)


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(qb, test_args.r)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.r)


def test_bad_distance(test_args: Args) -> None:
    rb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.q, rb)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.q, 100)
