from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    compressibility_factor_via_intermolecular_force_potential as compressibility_law,
)

# Description
## Assuming the model of hard spheres, for a gas with particles of radius r = 4 Å confined
## in a space of volume V = 1 m**3 with N = 3e25 number of particles, its compressibility factor
## is approximately Z = 1.004.

Args = namedtuple("Args", "n v r")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 3 * 10**25
    v = Quantity(1 * units.meter**3)
    r = Quantity(4 * units.angstrom)
    return Args(n=n, v=v, r=r)


def test_law(test_args: Args) -> None:
    result = compressibility_law.calculate_compressibility_factor(test_args.n, test_args.v, test_args.r)
    assert_equal(result, 1.004)


def test_bad_number_of_particles(test_args: Args) -> None:
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        compressibility_law.calculate_compressibility_factor(nb, test_args.v, test_args.r)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        compressibility_law.calculate_compressibility_factor(test_args.n, vb, test_args.r)
    with raises(TypeError):
        compressibility_law.calculate_compressibility_factor(test_args.n, 100, test_args.r)


def test_bad_radius(test_args: Args) -> None:
    rb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        compressibility_law.calculate_compressibility_factor(test_args.n, test_args.v, rb)
    with raises(TypeError):
        compressibility_law.calculate_compressibility_factor(test_args.n, test_args.v, 100)
