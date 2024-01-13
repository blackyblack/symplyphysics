from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import mass_flow_rate as mass_rate_definition


@fixture(name="test_args")
def test_args_fixture():
    m0 = Quantity(1 * units.kilograms)
    m1 = Quantity(20 * units.kilograms)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["m0", "m1", "t"])
    return Args(m0=m0, m1=m1, t=t)


def test_basic_mass_flow_rate(test_args):
    result = mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass / units.time)
    result_acceleration = convert_to(result, mass_rate_definition.definition_units_SI).evalf(2)
    assert result_acceleration == approx(3.8, 0.01)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(mb, test_args.m1, test_args.t)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, mb, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(100, test_args.m1, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, tb)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.m0, test_args.m1, 100)
