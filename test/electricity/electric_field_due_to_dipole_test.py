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
from symplyphysics.laws.electricity import electric_field_due_to_dipole as dipole_field

# Description
## A dipole of a magnitude of 0.5e-10 C*m produces an electric field of a magnitude of 0.899 N/C
## at a distance of 1 m.


@fixture(name="test_args")
def test_args_fixture():
    p = Quantity(0.5e-10 * units.coulomb * units.meter)
    z = Quantity(1 * units.meter)
    Args = namedtuple("Args", "p z")
    return Args(p=p, z=z)


def test_basic_law(test_args):
    result = dipole_field.calculate_electric_field(test_args.p, test_args.z)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force / units.charge)
    result_field = convert_to(result, units.newton / units.coulomb).evalf(3)
    assert_approx(result_field, 0.899)


def test_bad_dipole(test_args):
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        dipole_field.calculate_electric_field(pb, test_args.z)
    with raises(TypeError):
        dipole_field.calculate_electric_field(100, test_args.z)


def test_bad_distance(test_args):
    zb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_field.calculate_electric_field(test_args.p, zb)
    with raises(TypeError):
        dipole_field.calculate_electric_field(test_args.p, 100)
