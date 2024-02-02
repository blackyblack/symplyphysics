from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import errors, units, convert_to, Quantity, SI
from symplyphysics.laws.hydro import shear_stress_is_proportional_to_velocity_gradient as newtons_law

# Description
## Water (dynamic viscosity of which is 8.9e-4 Pa*s) is flowing steadily on top of a flat solid plate.
## At the site of contact with the plate its speed is 0, and at the top it is 0.1 m/s. The flow is 5 cm thick.
## Assert that the shear stress in the water at its highest point of contact with the plate amounts to 1.78e-3 Pa.


@fixture(name="test_args")
def test_args_fixture():
    dynamic_viscosity = Quantity(8.9e-4 * units.pascal * units.second)
    speed_before = Quantity(0 * units.meter / units.second)
    speed_after = Quantity(0.1 * units.meter / units.second)
    layer_separation = Quantity(5 * units.centimeter)
    Args = namedtuple("Args", "dynamic_viscosity speed_before speed_after layer_separation")
    return Args(
        dynamic_viscosity=dynamic_viscosity,
        speed_before=speed_before,
        speed_after=speed_after,
        layer_separation=layer_separation,
    )


def test_basic_law(test_args):
    result = newtons_law.calculate_shear_stress(
        test_args.dynamic_viscosity,
        test_args.speed_before,
        test_args.speed_after,
        test_args.layer_separation,
    )
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_stress = convert_to(result, units.pascal).evalf(3)
    assert result_stress == approx(1.78e-3, 1e-3)


def test_bad_dynamic_viscosity(test_args):
    bad_dynamic_viscosity = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_shear_stress(
            bad_dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
        )
    with raises(TypeError):
        newtons_law.calculate_shear_stress(
            100,
            test_args.speed_before,
            test_args.speed_after,
            test_args.layer_separation,
        )


def test_bad_fluid_speed(test_args):
    bad_fluid_speed = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            bad_fluid_speed,
            test_args.speed_after,
            test_args.layer_separation,
        )
    with raises(TypeError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            100,
            test_args.speed_after,
            test_args.layer_separation,
        )
    with raises(errors.UnitsError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            bad_fluid_speed,
            test_args.layer_separation,
        )
    with raises(TypeError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            100,
            test_args.layer_separation,
        )


def test_bad_layer_separation(test_args):
    bad_layer_separation = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            bad_layer_separation,
        )
    with raises(TypeError):
        newtons_law.calculate_shear_stress(
            test_args.dynamic_viscosity,
            test_args.speed_before,
            test_args.speed_after,
            100,
        )
