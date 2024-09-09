"""
Electromotive force induced in rotating coil
============================================

Suppose a coil is being rotated around the axis that lies in the coil's cross section
(see `Figure <https://www.schoolphysics.co.uk/age16-19/Electricity%20and%20magnetism/Electromagnetic%20induction/text/Induced_emf_in_a_rotating_coil/index.html>`__)
in a magnetic field under the conditions described below. Then an electromotive will
be induced in the contour of the coil. It will be proportional to the number of turns
in the coil and the rate of change of magnetic flux through the coil.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

