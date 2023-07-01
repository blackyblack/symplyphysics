from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.nuclear import infinite_multiplication_factor_from_macroscopic_fission_cross_section as multiplication_factor


@fixture(name="test_args")
def test_args_fixture():
    # Uranium-235 thermal neutrons per fission is 2.42
    neutrons_per_fission = 2.42
    # Uranium-235 macroscopic fission cross-section
    macro_fission_cross_section = Quantity(2.811 / units.centimeter)
    # Uranium-235 + Uranium-238 macroscopic absorption cross-section
    macro_abs_cross_section = Quantity((3.352 + 0.117) / units.centimeter)
    Args = namedtuple("Args", ["v", "Sf", "Sa"])
    return Args(v=neutrons_per_fission, Sf=macro_fission_cross_section, Sa=macro_abs_cross_section)


def test_basic_multiplication_factor(test_args):
    result = multiplication_factor.calculate_multiplication_factor(test_args.v, test_args.Sf,
        test_args.Sa)
    assert result == approx(1.96, 0.01)


def test_bad_macroscopic_cross_section(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        multiplication_factor.calculate_multiplication_factor(test_args.v, Sb, test_args.Sa)
    with raises(TypeError):
        multiplication_factor.calculate_multiplication_factor(test_args.v, 100, test_args.Sa)
    with raises(errors.UnitsError):
        multiplication_factor.calculate_multiplication_factor(test_args.v, test_args.Sf, Sb)
    with raises(TypeError):
        multiplication_factor.calculate_multiplication_factor(test_args.v, test_args.Sf, 100)
