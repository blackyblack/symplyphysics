from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units,
    errors,
    convert_to,
    Quantity,
    SI,
    QuantityVector,
)
from symplyphysics.laws.electricity.vector import vector_of_electric_dipole_moment as dipole_moment

# Description
## The vector of the dipole moment of two point charges of magnitude Q = 2e-13 C, the displacement 
## vector between which is [1, -1, 0.5] m, amounts to [2e-13, -2e-13, 1e-13] C*m


@fixture(name="test_args")
def test_args_fixture():
    q = Quantity(2e-13 * units.coulomb)
    l = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(-1.0 * units.meter),
        Quantity(0.5 * units.meter),    
    ])
    Args = namedtuple("Args", "q l p")
    return Args(q=q, l=l, p=p)


def test_basic_law(test_args):
    result = dipole_moment.calculate_dipole_moment(test_args.q, test_args.l)
    assert len(result.components) == 3
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.charge * units.length)
    correct_values = [2e-13, -2e-13, 1e-13]
    for result_component, correct_value in zip(result.components, correct_values):
        result_value = convert_to(result_component, units.coulomb * units.meter).evalf(2)
        assert result_value == approx(correct_value, 1e-2)


def test_bad_charge(test_args):
    qb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_moment.calculate_dipole_moment(qb, test_args.l)
    with raises(TypeError):
        dipole_moment.calculate_dipole_moment(100, test_args.l)


def test_bad_displacement(test_args):
    lb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_moment.calculate_dipole_moment(test_args.q, lb)
    with raises(TypeError):
        dipole_moment.calculate_dipole_moment(test_args.q, 100)
