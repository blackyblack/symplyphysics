from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic.rotational_inertia_formulas import (
    solid_disk_about_central_axis as solid_disk
)

# Description
## A disk is rotating about its central axis. It is 2.0 kg heavy and it has a radius of 0.4 m.
## Its rotational inertia about that axis amounts to 0.160 kg*(m**2).


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(2.0 * units.kilogram)
    r = Quantity(0.4 * units.meter)
    Args = namedtuple("Args", "m r")
    return Args(m=m, r=r)


def test_basic_law(test_args):
    result = solid_disk.calculate_rotational_inertia(test_args.m, test_args.r)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.length**2)
    result_value = convert_to(result, units.kilogram * units.meter**2).evalf(3)
    assert result_value == approx(0.16, 1e-3)


def test_bad_mass(test_args):
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        solid_disk.calculate_rotational_inertia(mb, test_args.r)
    with raises(TypeError):
        solid_disk.calculate_rotational_inertia(100, test_args.r)


def test_bad_radius(test_args):
    rb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        solid_disk.calculate_rotational_inertia(test_args.m, rb)
    with raises(TypeError):
        solid_disk.calculate_rotational_inertia(test_args.m, 100)
