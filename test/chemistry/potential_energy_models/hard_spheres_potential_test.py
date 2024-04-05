from collections import namedtuple
from pytest import fixture, raises
from sympy import S
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)
from symplyphysics.laws.chemistry.potential_energy_models import hard_spheres_potential

# Description
## With in the model of hard spheres of diameter sigma = 10 nm, the value of the potential
## is zero when the distance between the two particles is 20 nm, and infinity at distance r0 = 5 nm.

Args = namedtuple("Args", "r0 r1 d")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    r0 = Quantity(20 * units.nanometer)
    r1 = Quantity(5 * units.nanometer)
    d = Quantity(10 * units.nanometer)
    return Args(r0=r0, r1=r1, d=d)


def test_law_zero_potential(test_args: Args) -> None:
    result = hard_spheres_potential.calculate_potential(test_args.r0, test_args.d)
    assert_equal(result, 0)


# FIXME: maybe allow infinity too in [here?](https://github.com/blackyblack/symplyphysics/blob/ef7a449974e977a2676e55339d38f4943e7d1ad3/symplyphysics/core/dimensions.py#L107)
# Otherwise this raises a UnitsError
## def test_law_infinite_potential(test_args: Args) -> None:
##     result = hard_spheres_potential.calculate_potential(test_args.r1, test_args.d)
##     assert_equal(result, S.Infinity * units.joule)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        hard_spheres_potential.calculate_potential(lb, test_args.d)
    with raises(TypeError):
        hard_spheres_potential.calculate_potential(100, test_args.d)
    with raises(errors.UnitsError):
        hard_spheres_potential.calculate_potential(test_args.r0, lb)
    with raises(TypeError):
        hard_spheres_potential.calculate_potential(test_args.r0, 100)
