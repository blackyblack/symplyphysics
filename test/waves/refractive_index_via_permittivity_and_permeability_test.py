from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, Quantity, units, errors
from symplyphysics.laws.waves import refractive_index_via_permittivity_and_permeability as media_law

# Description.
## Refraction factor of air is 1.000292. Air dielectric permeability is 1.000576. Air magnetic permeability is 1.00000037.

Args = namedtuple("Args", ["dielectric_permeability", "magnetic_permeability"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    dielectric_permeability = 1.000576
    magnetic_permeability = 1.00000037
    return Args(dielectric_permeability=dielectric_permeability,
        magnetic_permeability=magnetic_permeability)


def test_basic_factor(test_args: Args) -> None:
    result = media_law.calculate_refraction_factor(test_args.dielectric_permeability,
        test_args.magnetic_permeability)
    assert_equal(result, 1.000292)


def test_bad_numbers(test_args: Args) -> None:
    nb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        media_law.calculate_refraction_factor(nb, test_args.magnetic_permeability)
    with raises(errors.UnitsError):
        media_law.calculate_refraction_factor(test_args.dielectric_permeability, nb)
