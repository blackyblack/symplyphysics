from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import electric_field_due_to_dipole as dipole_field

# Description
## A dipole of a magnitude of 0.5e-10 C*m produces an electric field of a magnitude of 0.899 N/C
## at a distance of 1 m.

Args = namedtuple("Args", "p z")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(5e-11 * units.coulomb * units.meter)
    z = Quantity(1 * units.meter)
    return Args(p=p, z=z)


def test_basic_law(test_args: Args) -> None:
    result = dipole_field.calculate_electric_field(test_args.p, test_args.z)
    assert_equal(result, 0.899 * units.newton / units.coulomb)


def test_bad_dipole(test_args: Args) -> None:
    pb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        dipole_field.calculate_electric_field(pb, test_args.z)
    with raises(TypeError):
        dipole_field.calculate_electric_field(100, test_args.z)


def test_bad_distance(test_args: Args) -> None:
    zb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_field.calculate_electric_field(test_args.p, zb)
    with raises(TypeError):
        dipole_field.calculate_electric_field(test_args.p, 100)
