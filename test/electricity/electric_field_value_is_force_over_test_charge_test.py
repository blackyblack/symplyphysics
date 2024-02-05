from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import electric_field_value_is_force_over_test_charge as electric_field

# The electric field is 6 N/C if the force exerted on the test charge of 0.5 C is 3 N.


@fixture(name="test_args")
def test_args_fixture():
    q0 = Quantity(0.5 * units.coulomb)
    F = Quantity(3 * units.newton)
    Args = namedtuple("Args", "q0 F")
    return Args(q0=q0, F=F)


def test_basic_law(test_args):
    result = electric_field.calculate_electric_field(test_args.F, test_args.q0)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force / units.charge)
    result_field = convert_to(result, units.newton / units.coulomb).evalf(2)
    assert_approx(result_field, 6)


def test_bad_force(test_args):
    Fb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(Fb, test_args.q0)
    with raises(TypeError):
        electric_field.calculate_electric_field(100, test_args.q0)


def test_bad_charge(test_args):
    q0b = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        electric_field.calculate_electric_field(test_args.F, q0b)
    with raises(TypeError):
        electric_field.calculate_electric_field(test_args.F, 100)
