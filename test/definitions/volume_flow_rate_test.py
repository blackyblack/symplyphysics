from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import volume_flow_rate as mass_rate_definition


@fixture(name="test_args")
def test_args_fixture():
    v0 = Quantity(1 * units.meters ** 3)
    v1 = Quantity(20 * units.meters ** 3)
    t = Quantity(5 * units.second)
    Args = namedtuple("Args", ["v0", "v1", "t"])
    return Args(v0=v0, v1=v1, t=t)


def test_basic_volume_flow_rate(test_args):
    result = mass_rate_definition.calculate_mass_flow_rate(test_args.v0, test_args.v1, test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume / units.time)
    result_acceleration = convert_to(result, mass_rate_definition.definition_units_SI).evalf(2)
    assert result_acceleration == approx(3.8, 0.01)


def test_bad_volume(test_args):
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(vb, test_args.v1, test_args.t)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.v0, vb, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(100, test_args.v1, test_args.t)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.v0, 100, test_args.t)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.v0, test_args.v1, tb)
    with raises(TypeError):
        mass_rate_definition.calculate_mass_flow_rate(test_args.v0, test_args.v1, 100)
