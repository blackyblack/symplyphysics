from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import errors, units, convert_to, Quantity, SI
from symplyphysics.laws.hydro import inner_friction_is_proportional_to_velocity_gradient as newtons_law


@fixture(name="test_args")
def test_args_fixture():
    dynamic_viscosity = Quantity(8.9e-4 * units.pascal * units.second)
    speed_before = Quantity(0 * units.meter / units.second)
    speed_after = Quantity(0.1 * units.meter / units.second)
    layer_separation = Quantity(5 * units.centimeter)
    layer_area = Quantity(30 * units.centimeter**2)
    Args = namedtuple("Args", "dynamic_viscosity speed_before speed_after layer_separation layer_area")
    return Args(
        dynamic_viscosity=dynamic_viscosity,
        speed_before=speed_before,
        speed_after=speed_after,
        layer_separation=layer_separation,
        layer_area=layer_area,
    )


def test_basic_law(test_args):
    result = newtons_law.calculate_inner_friction(
        test_args.dynamic_viscosity,
        test_args.speed_before,
        test_args.speed_after,
        test_args.layer_separation,
        test_args.layer_area,
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.force)
    result_force = convert_to(result, units.newton).evalf(6)
    assert result_force == approx(5.34e-6, 1e-6)


def test_bad_dynamic_viscosity(test_args):
    bad_dynamic_viscosity = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_inner_friction(
            bad_dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
            test_args.layer_area,
        )
    with raises(TypeError):
        newtons_law.calculate_inner_friction(
            100,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
            test_args.layer_area,
        )


def test_bad_fluid_speed(test_args):
    bad_fluid_speed = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            bad_fluid_speed,
            test_args.speed_after,
            test_args.layer_separation,
            test_args.layer_area,
        )
    with raises(TypeError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            100,
            test_args.speed_after,
            test_args.layer_separation,
            test_args.layer_area,
        )
    with raises(errors.UnitsError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            bad_fluid_speed,
            test_args.layer_separation,
            test_args.layer_area,
        )
    with raises(TypeError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            100,
            test_args.layer_separation,
            test_args.layer_area,
        )


def test_bad_layer_separation(test_args):
    bad_layer_separation = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            bad_layer_separation,
            test_args.layer_area,
        )
    with raises(TypeError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            100,
            test_args.layer_area,
        )


def test_bad_layer_area(test_args):
    bad_layer_area = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
            bad_layer_area,
        )
    with raises(TypeError):
        newtons_law.calculate_inner_friction(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
            100,
        )
