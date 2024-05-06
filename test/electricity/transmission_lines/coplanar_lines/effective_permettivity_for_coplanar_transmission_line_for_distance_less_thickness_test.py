from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.transmission_lines.coplanar_lines import effective_permettivity_for_coplanar_transmission_line_for_distance_less_thickness as permittivity_law

## The relative permittivity of the coplanar line dielectric is 4.
## Then the effective permittivity is 2.5.

Args = namedtuple("Args", [
    "relative_permittivity"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 4
    return Args(relative_permittivity=relative_permittivity)


def test_basic_effective_permittivity(test_args: Args) -> None:
    result = permittivity_law.calculate_effective_permittivity(test_args.relative_permittivity)
    assert_equal(result, 2.5)


def test_bad_permittivity() -> None:
    bad_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        permittivity_law.calculate_effective_permittivity(bad_permittivity)
