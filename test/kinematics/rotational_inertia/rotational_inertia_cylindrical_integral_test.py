from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematics.rotational_inertia import (
    rotational_inertia_cylindrical_integral as rotational_inertia_law,)

# Description
## A cylinder is rotating about its axis. Its density is 300 kg/m**3, and it has a radius of 0.3 m
## and a height of 0.1 m. Its rotational inertia amounts to 0.382 kg*(m**2).

Args = namedtuple("Args", "rho rmax zmax")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    rho = Quantity(300.0 * units.kilogram / units.meter**3)
    rmax = Quantity(0.3 * units.meter)
    zmax = Quantity(0.1 * units.meter)
    return Args(rho=rho, rmax=rmax, zmax=zmax)


def test_basic_law(test_args: Args) -> None:
    result = rotational_inertia_law.calculate_cylinder_rotational_inertia(
        test_args.rho, test_args.rmax, test_args.zmax)
    assert_equal(result, 0.382 * units.kilogram * units.meter**2)


def test_bad_density(test_args: Args) -> None:
    rhob = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(rhob, test_args.rmax,
            test_args.zmax)
    with raises(TypeError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(100, test_args.rmax,
            test_args.zmax)


def test_bad_lengths(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(test_args.rho, lb,
            test_args.zmax)
    with raises(errors.UnitsError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(test_args.rho, test_args.rmax,
            lb)
    with raises(TypeError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(test_args.rho, 100,
            test_args.zmax)
    with raises(TypeError):
        rotational_inertia_law.calculate_cylinder_rotational_inertia(test_args.rho, test_args.rmax,
            100)
