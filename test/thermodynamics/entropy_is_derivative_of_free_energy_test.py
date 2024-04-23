from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import entropy_is_derivative_of_free_energy as entropy_law

# Description
## As the system's temperature increases by 1 mK, its free energy decreases by 1 J. The entropy of the
## system is 1 kJ/K.

Args = namedtuple("Args", "df dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    df = Quantity(-1 * units.joule)
    dt = Quantity(1 * prefixes.milli * units.kelvin)
    return Args(df=df, dt=dt)


def test_law(test_args: Args) -> None:
    result = entropy_law.calculate_entropy(test_args.df, test_args.dt)
    assert_equal(result, 1 * prefixes.kilo * units.joule / units.kelvin)


def test_bad_energy(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_entropy(fb, test_args.dt)    
    with raises(TypeError):
        entropy_law.calculate_entropy(100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_entropy(test_args.df, tb)
    with raises(TypeError):
        entropy_law.calculate_entropy(test_args.df, 100)
