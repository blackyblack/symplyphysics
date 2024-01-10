from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.hydro import body_weight_in_liquid as weight


@fixture(name="test_args")
def test_args_fixture():
    weight_air = Quantity(2.6 * units.newton)
    liquid_density = Quantity(1000 * units.kilogram / units.meter**3)
    body_density = Quantity(2600 * units.kilogram / units.meter**3)
    Args = namedtuple("Args", ["weight_air", "liquid_density", "body_density"])
    return Args(
        weight_air=weight_air,
        liquid_density=liquid_density,
        body_density=body_density
               )


def test_basic_weight(test_args):
    result = weight.calculate_weight(test_args.weight_air, test_args.liquid_density, test_args.body_density)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_weight = convert_to(result, units.newton).evalf(5)
    assert result_weight == approx(1.6, 0.001)


def test_bad_density(test_args):
    at = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        weight.calculate_weight(test_args.weight_air, at, test_args.body_density)
    with raises(errors.UnitsError):
        weight.calculate_weight(test_args.weight_air, test_args.liquid_density, at)
    with raises(TypeError):
        weight.calculate_weight(test_args.weight_air, 100, test_args.body_density)
    with raises(TypeError):
        weight.calculate_weight(test_args.weight_air, test_args.liquid_density, 100)


def test_bad_weight(test_args):
    bad_weight = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        weight.calculate_weight(bad_weight, test_args.liquid_density, test_args.body_density)
    with raises(TypeError):
        weight.calculate_weight(100, test_args.liquid_density, test_args.body_density)
