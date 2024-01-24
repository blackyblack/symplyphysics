from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import superposition_of_forces_is_sum as forces_law


@fixture(name="test_args")
def test_args_fixture():
    F1 = Quantity(10 * units.newton)
    F2 = Quantity(20 * units.newton)
    Args = namedtuple("Args", ["F1", "F2"])
    return Args(F1=F1, F2=F2)


def test_basic_superposition(test_args):
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2])
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        forces_law.definition_units_SI)
    result_force = convert_to(result, units.newton).evalf(3)
    assert result_force == approx(30, 0.001)


def test_three_forces_array(test_args):
    F3 = Quantity(-5 * units.newton)
    result = forces_law.calculate_resultant_force([test_args.F1, test_args.F2, F3])
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
        forces_law.definition_units_SI)
    result_force = convert_to(result, units.newton).evalf(3)
    assert result_force == approx(25, 0.01)


def test_bad_force(test_args):
    Fb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, test_args.F2])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, test_args.F2])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([test_args.F1, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([test_args.F1, 100])
    with raises(errors.UnitsError):
        forces_law.calculate_resultant_force([Fb, Fb])
    with raises(TypeError):
        forces_law.calculate_resultant_force([100, 100])
    with raises(TypeError):
        forces_law.calculate_resultant_force(test_args.F1)
