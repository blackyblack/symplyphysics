from sympy import simplify
from sympy.vector import laplacian
from symplyphysics import (
    symbols, Function, Eq, pretty
)

# Description
## The quantity Bg^2 is called the geometrical buckling of the reactor and depends only on the geometry.
## This term is derived from the notion that the neutron flux distribution is somehow "buckled" in a homogeneous
## finite reactor. The neutron flux has more concave downward or "buckled" curvature (higher Bg^2) in a small
## reactor than in a large one. 

## Law: Bg^2 = ∇^2(Ф(x)) / Ф(x)
## Where:
## ∇^2 - Laplacian operator.
## Ф(x) (neutron flux density) - number of neutrons crossing through some arbitrary cross-sectional unit area in all
##   directions per unit time.
## Bg^2 - geometric buckling.

geometric_buckling = symbols('geometric_buckling')
distance_from_center = symbols('distance_from_center')
neutron_flux_function = symbols('neutron_flux_function', cls = Function)

flux_laplacian_function = laplacian(neutron_flux_function)

# see geometric buckling definition
law = simplify(Eq(geometric_buckling,
    -1 * flux_laplacian_function / neutron_flux_function(distance_from_center)))

def print():
    return pretty(law, use_unicode=False)
