from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import electric_field_value_is_force_over_test_charge as electric_field

# The electric field is 6 N/C if the force exerted on the test charge of 0.5 C is 3 N.

Args = namedtuple("Args", "q0 F")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q0 = Quantity(0.5 * units.coulomb)
    F = Quantity(3 * units.newton)
    return Args(q0=q0, F=F)


def test_basic_law(test_args: Args) -> None:
    result = electric_field.calculate_electric_field(test_args.F, test_args.q0)
    assert_equal(result, 6 * units.newton / units.coulomb)


def test_bad_force(test_args: Args) -> None:
    Fb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(Fb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.q0)


def test_bad_charge(test_args: Args) -> None:
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.F, q0b)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.F, 100)
