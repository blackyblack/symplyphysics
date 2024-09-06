from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    units,
    Quantity,
    errors,
)
from symplyphysics.laws.electricity import energy_via_constant_power_and_time as law

Args = namedtuple("Args", "p t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(4.5 * units.watt)
    t = Quantity(20 * units.second)
    return Args(p=p, t=t)


def test_law(test_args: Args) -> None:
    result = law.calculate_energy(test_args.p, test_args.t)
    assert_equal(result, 90 * units.joule)


def test_bad_power(test_args: Args) -> None:
    pb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_energy(pb, test_args.t)
    with raises(TypeError):
        law.calculate_energy(100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_energy(test_args.p, tb)
    with raises(TypeError):
        law.calculate_energy(test_args.p, 100)
