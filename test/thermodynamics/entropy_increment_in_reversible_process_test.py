from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)

from symplyphysics.laws.thermodynamics import entropy_increment_in_reversible_process as entropy_law

# Description
## The entropy increment in a reversible prossess of heat exchange where the amount of heat transferred
## is δQ = 0.1 J and the common temperature of the systems is 200 K, is dS = 0.5 mJ/K

Args = namedtuple("Args", "dq t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dq = Quantity(0.1 * units.joule)
    t = Quantity(200 * units.kelvin)
    return Args(dq=dq, t=t)


def test_law(test_args: Args) -> None:
    result = entropy_law.calculate_infinitesimal_entropy_change(test_args.dq, test_args.t)
    assert_equal(result, 0.5 * prefixes.milli * units.joule / units.kelvin)


def test_bad_heat(test_args: Args) -> None:
    qb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_infinitesimal_entropy_change(qb, test_args.t)
    with raises(TypeError):
        entropy_law.calculate_infinitesimal_entropy_change(100, test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        entropy_law.calculate_infinitesimal_entropy_change(test_args.dq, tb)
    with raises(TypeError):
        entropy_law.calculate_infinitesimal_entropy_change(test_args.dq, 100)
