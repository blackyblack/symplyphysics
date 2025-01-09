"""
Quantum Mechanics (Symbols)
===========================

Symbols related to quantum mechanics.
"""

from sympy import sqrt
from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

wave_function = SymbolNew("psi", 1 / sqrt(units.length), display_latex="\\psi")
"""
The **wave function** is a complex-valued quantity describing the quantum state of an isolated quantum system.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Wave_function>`__.
"""

probability_density = SymbolNew("rho", 1 / units.length, display_latex="\\rho")
"""
The **probability density**, or **probability amplitude**, is a real number that provide a relationship between
the quantum state vector of a system and the results of observations of that system. It is the likelihood of
finding an object or system in a particular state at a particular position and time.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Probability_amplitude>`__.
#. `AK Lectures <https://aklectures.com/lecture/schrodinger%27s-equations/probability-density-of-particles>`__.
"""
