"""
Nuclear Physics (Symbols)
=========================

Symbols related to electrodynamics.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

geometric_buckling = Symbol("B_g^2", 1 / units.area, display_latex="B_\\text{g}^2")
"""
**Geometric buckling** is a measure of neutron leakage in a nuclear reactor.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#>`__.
"""

multiplication_factor = Symbol("k", dimensionless)
"""
The **multiplication factor** denotes the rate of change of the neutron population in a
system, and is the ratio of the neutron population in the following generation to the
neutron population in the given generation.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Four_factor_formula#Multiplication>`__.
"""

infinite_multiplication_factor = Symbol("k_inf", dimensionless, display_latex="k_\\infty")
"""
The **infinite multiplication factor** is the :attr:`~multiplication_factor` when the
medium is infinite so that neutrons cannot leak out of the system.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Four_factor_formula#Multiplication>`__.
"""

effective_multiplication_factor = Symbol("k_eff", dimensionless, display_latex="k_\\text{eff}")
"""
The **effective multiplication factor** is most often defined as the ratio of the rate
of neutron production to the rate of neutron loss in a nuclear system.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Nuclear_chain_reaction#Effective_neutron_multiplication_factor>`__.
"""

neutron_diffusion_area = Symbol("L^2", units.area)
"""
**Diffusion area** is a quantity that appears when calculating the average distance
between the neutron's birth point as a thermal neutron and its absorption.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-length/>`__.
"""

neutron_flux = Symbol("Phi", 1 / (units.area * units.time), display_latex="\\Phi")
"""
**Neutron flux** is a scalar quantity defined as the total distance traveled by all free
neutrons per unit time and volume.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Neutron_flux#>`__.
"""

material_buckling = Symbol("B_m^2", 1 / units.area, display_latex="B_\\text{m}^2")
"""
**Material buckling** is a measure of the difference between neutron production and
neutron absorption in a nuclear reactor.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Geometric_and_material_buckling#>`__.
"""

fast_non_leakage_probability = Symbol("P_FNL", dimensionless, display_latex="P_\\text{FNL}")
"""
**Fast non-leakage probability** is defined as the ratio of the number of fast neutrons
that do not leak from the reactor to the number of fast neutrons produced by all
fissions.
"""

thermal_non_leakage_probability = Symbol("P_TNL", dimensionless, display_latex="P_\\text{TNL}")
"""
**Thermal non-leakage probability** is defined as the ratio of the number of thermal
neutrons that do not leak from the reactor to the number of thermal neutrons produced by
all fissions.
"""

fast_fission_factor = Symbol("epsilon", dimensionless, display_latex="\\varepsilon")
"""
**Fast fission factor** is the ratio of total number of fission neutrons to the number
of fission neutrons from just thermal fissions.

**Links:**

#. `Wikipedia, see row 4 of the table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

resonance_escape_probability = Symbol("p", dimensionless)
"""
**Resonance escape probability** is the probability that a neutron will slow down from
fission energy to thermal energies without being captured by a nuclear resonance.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Resonance_escape_probability#>`__.
"""

fast_absorption_fission_probability = Symbol("P_FAF", dimensionless, display_latex="P_\\text{FAF}")
"""
**Fast absorption fission probability** is the probability that a fast neutron
absorption in fuel causes fission.

**Links:**

#. `Wikipedia, see list under the table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

thermal_absorption_fission_probability = Symbol("P_TAF", dimensionless, display_latex="P_\\text{TAF}")
"""
**Thermal absorption fission probability** is the probability that a thermal neutron
absorption in fuel causes fission.

**Links:**

#. `Wikipedia, see list under the table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

thermal_utilization_factor = Symbol("f", dimensionless)
"""
**Thermal utilization factor** is defined as the ratio of the number of neutrons
absorbed by the fuel isotope to the number of neutrons absorbed anywhere.

**Links:**

#. `Wikipedia, second row of the table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

fast_utilization = Symbol("u_f", dimensionless, display_latex="u_\\text{f}")
"""
**Fast utilization** is the probability that a fast neutron is absorbed in fuel.

**Links:**

#. `Wikipedia, see list under the table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

neutron_fermi_age = Symbol("tau", units.area, display_latex="\\tau")
"""
**Fermi age** is a measure of how far a neutron travels during moderation (e.g. in a
graphite moderator), similar to the :symbols:`neutron_diffusion_area`.
"""

thermal_fission_factor = Symbol("eta", dimensionless, display_latex="\\eta")
"""
**Thermal fission factor** is defined as the ratio of the number of neutrons produced
from fission to the absorption in fuel isotope.

**Links:**

#. `Wikipedia, first line in table <https://en.wikipedia.org/wiki/Six_factor_formula#>`__.
"""

half_life = Symbol("t_1/2", units.time, display_latex="t_{1/2}")
"""
**Half-life** is the time required for a quantity to reduce to half of its initial
value.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Half-life#>`__.
"""

macroscopic_cross_section = Symbol("Sigma", 1 / units.length, display_latex="\\Sigma")
"""
The **macroscopic cross section** is an interaction rate per unit distance traveled by
the neutron obtained by multiplying the microscopic cross section (expressed as an area)
by atomic density (expressed as reciprocal volume).

**Links:**

#. `ScienceDirect <https://www.sciencedirect.com/topics/engineering/macroscopic-cross-section>`__.
"""

neutron_diffusion_coefficient = Symbol("D", units.length)
"""
**Diffusion coefficient** is the proportionality constant between the neutron current
density and the neutron current, as per the first Fick's law of diffusion.

**Links:**

#. `NuclearPower <https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/>`__.
"""

migration_area = Symbol("M^2", units.area)
"""
**Migration area**, or the square of the migration length, is the measure of the
distance a neutron travels while slowing down as a fast neutron and diffusing as a
thermal neutron.
"""

reproduction_factor = Symbol("eta", dimensionless, display_latex="\\eta")
"""
**Reproduction factor** is defined as the ratio of the number of neutrons produced from
thermal fissions to the thermal absorption in fuel isotope.

**Links:**

#. `Wikipedia, first row in table <https://en.wikipedia.org/wiki/Four_factor_formula>`__.
"""
