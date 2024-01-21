from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import rotational_inertia_is_additive as rotational_inertia_law

# Description
## The system consists of two particles, the rotational inertia of which is 0.1 kg*m^2 and
## 3.0 kg*m^2 respectively. The total rotational inertia of the system amounts to 3.1 kg*m^2.
## If we add a third particle with rotational inertia of 2.1 kg*m^2, the total rotational
## inertia will be 5.2 kg*m^2.


@fixture(name="test_args")
def test_args_fixture():
    I1 = Quantity(0.1 * units.kilogram * units.meter**2)
    I2 = Quantity(3.0 * units.kilogram * units.meter**2)
    I3 = Quantity(2.1 * units.kilogram * units.meter**2)
    Args = namedtuple("Args", "I1 I2 I3")
    return Args(I1=I1, I2=I2, I3=I3)


def test_law_two_particles(test_args):
    result = rotational_inertia_law.calculate_rotational_inertia([
        test_args.I1, test_args.I2
    ])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert result_value == approx(3.1, 1e-2)


def test_law_three_particles(test_args):
    result = rotational_inertia_law.calculate_rotational_inertia([
        test_args.I1, test_args.I2, test_args.I3
    ])
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert result_value == approx(5.2, 1e-2)


def test_bad_inertias(test_args):
    Ib = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia([Ib])
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia([test_args.I2, Ib])
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia(Ib)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.I1)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(100)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia([100])
