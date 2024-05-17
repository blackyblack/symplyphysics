from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.electricity.resonators import resonant_frequency_of_ring_resonator as frequency_law

# Description
## The length of the ring resonator is 850 micrometer. The effective dielectric permittivity is 4. The order of interference is 1.
## Then the resonant frequency is 176.35 gigahertz.

Args = namedtuple("Args", ["ring_length", "order_interference", "permittivity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ring_length = Quantity(850 * units.micrometer)
    order_interference = 1
    permittivity = 4

    return Args(ring_length=ring_length,
        order_interference=order_interference,
        permittivity=permittivity)


def test_basic_frequency(test_args: Args) -> None:
    result = frequency_law.calculate_frequency(test_args.ring_length, test_args.order_interference,
        test_args.permittivity)
    assert_equal(result, 176.35 * prefixes.giga * units.hertz)


def test_bad_ring_length(test_args: Args) -> None:
    ring_length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_frequency(ring_length, test_args.order_interference,
            test_args.permittivity)
    with raises(TypeError):
        frequency_law.calculate_frequency(100, test_args.order_interference, test_args.permittivity)


def test_bad_order_interference(test_args: Args) -> None:
    order_interference = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_frequency(test_args.ring_length, order_interference,
            test_args.permittivity)


def test_bad_permittivity(test_args: Args) -> None:
    permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        frequency_law.calculate_frequency(test_args.ring_length, test_args.order_interference,
            permittivity)
