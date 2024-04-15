from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.bose_einstein_statistics import (
    single_particle_state_distribution as distribution_law,
)

# Description
## Assuming the particles follow Bose-Einstein distribution, the occupancy of the single-particle state
## with energy e = 4.01 eV is 0.456. The total chemical potential of the system is mu = 4.00 eV, and the
## temperature of the system is 100 K.

Args = namedtuple("Args", "e mu t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(4.01 * units.electronvolt)
    mu = Quantity(4.00 * units.electronvolt)
    t = Quantity(100 * units.kelvin)
    return Args(e=e, mu=mu, t=t)


def test_law(test_args: Args) -> None:
    result = distribution_law.calculate_occupancy_of_state(test_args.e, test_args.mu, test_args.t)
    assert_equal(result, 0.456)


def test_bad_energy(test_args: Args) -> None:
    eb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distribution_law.calculate_occupancy_of_state(eb, test_args.mu, test_args.t)
    with raises(TypeError):
        distribution_law.calculate_occupancy_of_state(100, test_args.mu, test_args.t)
    with raises(errors.UnitsError):
        distribution_law.calculate_occupancy_of_state(test_args.e, eb, test_args.t)
    with raises(TypeError):
        distribution_law.calculate_occupancy_of_state(test_args.e, 100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        distribution_law.calculate_occupancy_of_state(test_args.e, test_args.mu, tb)
    with raises(TypeError):
        distribution_law.calculate_occupancy_of_state(test_args.e, test_args.mu, 100)
