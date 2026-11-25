from sympy import limit
from sympy.physics.units.systems.si import dimsys_SI
from symplyphysics.core.dimensions.collect_expression import collect_expression_and_dimension
from symplyphysics.nuclear_physics.buckling import neutron_flux_for_uniform_sphere as sphere


def test_flux_dimension() -> None:
    _, dimension = collect_expression_and_dimension(sphere.law.rhs)
    assert dimsys_SI.equivalent_dims(dimension, sphere.neutron_flux.dimension)


def test_center_flux() -> None:
    assert limit(sphere.law.rhs, sphere.radial_distance, 0, dir="+") == sphere.dimension_factor


def test_boundary_flux() -> None:
    assert sphere.law.rhs.subs(sphere.radial_distance, sphere.radius) == 0
