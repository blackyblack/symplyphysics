from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    poisson_ratio_is_transverse_to_axial_strain_ratio as law)

Args = namedtuple("Args", "et ea")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    et = 0.001
    ea = -0.002
    return Args(et=et, ea=ea)


def test_law(test_args: Args) -> None:
    result = law.calculate_poisson_ratio(test_args.et, test_args.ea)
    assert_equal(result, 0.5)


def test_bad_strain(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_poisson_ratio(eb, test_args.ea)
    with raises(errors.UnitsError):
        law.calculate_poisson_ratio(test_args.et, eb)
