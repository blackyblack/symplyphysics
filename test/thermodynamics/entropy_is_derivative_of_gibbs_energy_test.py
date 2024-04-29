from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import entropy_is_derivative_of_gibbs_energy as entropy_law

# Description
## As the system's temperature increases by 1 mK, its Gibbs energy decreases by 1 J. The entropy of the
## system is 1 kJ/K.

Args = namedtuple("Args", "dg dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dg = Quantity(-1 * units.joule)
    dt = Quantity(1 * prefixes.milli * units.kelvin)
    return Args(dg=dg, dt=dt)


def test_law(test_args: Args) -> None:
    result = entropy_law.calculate_entropy(test_args.dg, test_args.dt)
    assert_equal(result, 1 * prefixes.kilo * units.joule / units.kelvin)


def test_bad_energy(test_args: Args) -> None:
    gb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_entropy(gb, test_args.dt)    
    with raises(TypeError):
        entropy_law.calculate_entropy(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_entropy(test_args.dg, tb)
    with raises(TypeError):
        entropy_law.calculate_entropy(test_args.dg, 100)