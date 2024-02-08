from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic.rotational_inertia import (
    rotational_inertia_cartesian_integral as rotational_inertia_law,)

# Description
## The rotational inertia of a unit cube about an axis passing through its center and the centers of
## the top and bottom faces is 1.0 kg*m^2. Its density is 6.0 kg/m^3.


@fixture(name="test_args")
def test_args_fixture():
    x = y = z = Quantity(0.5 * units.meter)
    rho = Quantity(6.0 * units.kilogram / units.meter**3)
    Args = namedtuple("Args", "x y z rho")
    return Args(x=x, y=y, z=z, rho=rho)


def test_basic_law(test_args):
    result = rotational_inertia_law.calculate_rotational_inertia(test_args.rho, test_args.x,
        test_args.y, test_args.z)
    assert_equal(result, 1 * units.kilogram * units.meter**2)


def test_bad_density(test_args):
    rhob = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia(rhob, test_args.x, test_args.y,
            test_args.z)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(100, test_args.x, test_args.y,
            test_args.z)


def test_bad_lengths(test_args):
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, lb, test_args.y,
            test_args.z)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, test_args.x, lb,
            test_args.z)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, test_args.x, test_args.y,
            lb)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, 100, test_args.y,
            test_args.z)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, test_args.x, 100,
            test_args.z)
    with raises(TypeError):
        rotational_inertia_law.calculate_rotational_inertia(test_args.rho, test_args.x, test_args.y,
            100)
