from collections import namedtuple
from pytest import approx, fixture
from symplyphysics.laws.waves import refraction_factor_from_media as media_law

# Description.
## Refraction factor of air is 1.000292. Air dielectric permeability is 1.000576. Air magnetic permeability is 1.00000037.


@fixture
def test_args():
    dielectric_permeability = 1.000576
    magnetic_permeability = 1.00000037
    Args = namedtuple("Args", ["dielectric_permeability", "magnetic_permeability"])
    return Args(dielectric_permeability=dielectric_permeability,
        magnetic_permeability=magnetic_permeability)


def test_basic_factor(test_args):
    result = media_law.calculate_refraction_factor(test_args.dielectric_permeability,
        test_args.magnetic_permeability).evalf(6)
    assert result == approx(1.000292, 0.00001)
