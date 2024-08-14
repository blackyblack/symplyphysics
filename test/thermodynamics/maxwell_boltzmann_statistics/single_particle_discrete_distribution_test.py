from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import single_particle_discrete_distribution as law

# Description
## A classical system of non-interacting particles in thermal equilibrium contains 1000 particles.
## The equilibrium temperature of the system is 293.15 K. Assuming the partition function is Z = 2.53,
## the average of particles found in a microstate with energy E_i = 3e-21 J is 188.

Args = namedtuple("Args", "n e t z")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    n = 1000
    e = Quantity(3e-21 * units.joule)
    t = Quantity(293.15 * units.kelvin)
    z = 2.53
    return Args(n=n, e=e, t=t, z=z)


def test_law(test_args: Args) -> None:
    result = law.calculate_particle_count_in_microstate(
        test_args.n, test_args.e, test_args.t, test_args.z)
    assert_equal(result, 188, tolerance=2e-3)


def test_bad_particle_count(test_args: Args) -> None:
    nb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_particle_count_in_microstate(nb, test_args.e, test_args.t,
            test_args.z)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_particle_count_in_microstate(test_args.n, eb, test_args.t,
            test_args.z)
    with raises(TypeError):
        law.calculate_particle_count_in_microstate(test_args.n, 100, test_args.t,
            test_args.z)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_particle_count_in_microstate(test_args.n, test_args.e, tb,
            test_args.z)
    with raises(TypeError):
        law.calculate_particle_count_in_microstate(test_args.n, test_args.e, 100,
            test_args.z)


def test_bad_partition_function(test_args: Args) -> None:
    zb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_particle_count_in_microstate(test_args.n, test_args.e,
            test_args.t, zb)
