from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.fermi_dirac_statistics import (
    single_particle_state_distribution as state_distribution,)

# Description
## The average number of fermions in a single-particle state of energy e = 4.01 eV is 0.24.
## The total chemical potential of the system is mu = 4 eV and the temperature of the system
## is 100 K.

Args = namedtuple("Args", "e mu t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(4.01 * units.electronvolt)
    mu = Quantity(4 * units.electronvolt)
    t = Quantity(100 * units.kelvin)
    return Args(e=e, mu=mu, t=t)


def test_law(test_args: Args) -> None:
    result = state_distribution.calculate_occupancy_of_state(test_args.e, test_args.mu, test_args.t)
    assert_equal(result, 0.24, tolerance=6e-3)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        state_distribution.calculate_occupancy_of_state(eb, test_args.mu, test_args.t)
    with raises(TypeError):
        state_distribution.calculate_occupancy_of_state(100, test_args.mu, test_args.t)
    with raises(errors.UnitsError):
        state_distribution.calculate_occupancy_of_state(test_args.e, eb, test_args.t)
    with raises(TypeError):
        state_distribution.calculate_occupancy_of_state(test_args.e, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        state_distribution.calculate_occupancy_of_state(test_args.e, test_args.mu, tb)
    with raises(TypeError):
        state_distribution.calculate_occupancy_of_state(test_args.e, test_args.mu, 100)
