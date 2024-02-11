from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import mass_flow_rate as mass_rate_definition

Args = namedtuple("Args", ["m0", "m1", "t"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m0 = Quantity(1 * units.kilograms)
    m1 = Quantity(20 * units.kilograms)
    t = Quantity(5 * units.second)
    return Args(m0=m0, m1=m1, t=t)


def test_basic_mass_flow_rate(test_args: Args) -> None:
    result = mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, test_args.t)
    assert_equal(result, 3.8 * units.kilogram / units.second)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(mb, test_args.m1, test_args.t)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, mb, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(100, test_args.m1, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, 100, test_args.t)


def test_bad_time(test_args: Args) -> None:
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, tb)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, 100)
