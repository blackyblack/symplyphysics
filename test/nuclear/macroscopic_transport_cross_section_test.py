from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.nuclear import macroscopic_transport_cross_section as macro_tr_cs


@fixture(name="test_args")
def test_args_fixture():
    # carbon macroscopic scattering cross-section at 1 eV is 0.528 cm^-1
    macroscopic_scattering_cross_section = Quantity(0.528 / units.centimeter)
    # carbon (mass_number = 12) average scattering cosine angle is 0.0555
    average_scattering_angle_cosine = 0.0555
    Args = namedtuple("Args", ["scatter_cs", "u"])
    return Args(scatter_cs=macroscopic_scattering_cross_section, u=average_scattering_angle_cosine)


def test_basic_cross_section(test_args):
    result = macro_tr_cs.calculate_cross_section(test_args.scatter_cs, test_args.u)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, 1 / units.length)
    result_cross_section = convert_to(result, 1 / units.centimeter).evalf(2)
    assert result_cross_section == approx(0.498, 0.01)


def test_bad_scattering_cross_section(test_args):
    scatter_csb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        macro_tr_cs.calculate_cross_section(scatter_csb, test_args.u)
    with raises(TypeError):
        macro_tr_cs.calculate_cross_section(100, test_args.u)
