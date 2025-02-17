"""
Chemistry (Symbols)
===================

Symbols of physical quantities related to chemistry.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import Symbol

mass_fraction = Symbol("w", dimensionless)
"""
**Mass fraction** of a substance within a mixture is the ratio of the mass of the substance to the total mass
of the mixture.
"""

amount_of_substance = Symbol("n", units.amount_of_substance)
"""
**Amount of substance** in a given sample of matter is defined as a ratio between the number of elementary
entities and the Avogadro constant.
"""

density_of_states = Symbol("D", 1 / units.volume)
r"""
The **density of states** of a system describes describes the number of allowed modes or states
per unit energy range.
"""

band_gap = Symbol("E_g", units.energy, display_latex="E_\\text{g}")
"""
A **band gap**, or **energy gap**, is an energy range in a solid where no electronic states exist. 
"""

work_function = Symbol("W", units.energy)
"""
**Work function** is the minimum thermodynamic work (i.e., energy) needed to remove an electron from
a solid to a point in the vacuum immediately outside the solid surface.
"""

drift_velocity = Symbol("u", units.velocity)
"""
**Drift velocity** is the average velocity attained by charged particles, such as electrons, in a material
due to an electric field.
"""

molar_concentration = Symbol("c", units.amount_of_substance / units.volume)
"""
**Molar concentration**, or **molarity**, is a quantity most commonly defined as amount of substance of
solute per unit volume of solution, or per unit volume available to the species.
"""

ionization_coefficient = Symbol("alpha", 1 / units.length, display_latex="\\alpha")
"""
**Ionization coefficient** can be defined as the mean number of ionization processes
over the distance covered in the direction of the electric field.

**Links:**

#. `ETH Research Collection <https://www.research-collection.ethz.ch/bitstream/20.500.11850/186582/1/PostPrint.pdf>`__.
"""

cross_section = Symbol("sigma", units.area, display_latex="\\sigma")
"""
**Cross section** is a measure of the probability that a specific process will take place
in a collision of two particles.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Cross_section_(physics)>`__.
"""

atomic_number = Symbol("Z", dimensionless, integer=True, positive=True)
"""
The **atomic number** or **nuclear charge number** of a chemical element is the charge
number of its atomic nucleus.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Atomic_number>`__.
"""

valence = Symbol("v", dimensionless, nonnegative=True)
"""
**Valence** of an atom is a measure of its combining capacity with other atoms when it
forms chemical compounds or molecules, and is generally understood to be the number of
chemical bonds that each atom of a given chemical element typically forms.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Valence_(chemistry)>`__.
"""

equilibrium_constant = Symbol("K", dimensionless)
"""
The **equilibrium constant** expresses the relationship between products and reactants
of a reaction at equilibrium with respect to a specific unit.

**Links:**

#. `Chemistry LibreTexts <https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Equilibria/Chemical_Equilibria/The_Equilibrium_Constant>`__.
"""

electrochemical_equivalent = Symbol("Z", units.mass / units.charge)
"""
The **electrochemical equivalent** of a chemical element is the mass of that element
transported by a specific quantity of electricity, usually charge.
"""

mobility = Symbol("mu", units.area / units.voltage / units.time, display_latex="\\mu")
"""
Electrical **mobility** is the ability of charged particles to move through a medium in
response to an electric field that is pulling them.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electrical_mobility>`__.
"""

diffusion_coefficient = Symbol("D", units.area / units.time)
"""
**Diffusion coefficient**, also referred to as (mass) **diffusivity**, is the
proportionality constant between the molar flux due to molecular diffusion and the
negative value of the gradient in the concentration of the species. In simpler terms,
it is the amount of a particular substance that diffuses across a unit area in unit time
under the influence of a gradient of one unit.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_diffusivity>`__.
#. `ScienceDirect <https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/diffusion-coefficient>`__.
"""

mass_number = Symbol("A", dimensionless, integer=True)
"""
**Mass number**, also called **atomic mass number** or **nucleon number**, is the total
number of protons and neutrons in an atomic nucleus.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_number>`__.
"""
