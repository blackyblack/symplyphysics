"""
Relativistic mechanics (Symbols)
================================

Symbols related to relativistic mechanics.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

lorentz_factor = SymbolNew("gamma", dimensionless, display_latex="\\gamma")
"""
**Lorentz factor** is a quantity expressing how much the measurements of time, length, and other physical
properties change for an object while it moves.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Lorentz_factor#>`__.
"""

spacetime_interval = SymbolNew("s", units.length)
"""
The **spacetime interval** between two events in four-dimensional spacetime is the analog of the
distance between two points in three-dimensional space. It has a property of being invariant under
the change of reference frame.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Spacetime#Spacetime_interval>`__.
"""

proper_time = SymbolNew("tau", units.time, display_latex="\\tau")
"""
**Proper time** along a timelike world line is defined as the time as measured by a clock following
that line. The change in proper time, called **proper time interval**, between two events on a world
line has the property of being Lorentz invariant and is independent of coordinates.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Proper_time#>`__.
"""

rest_mass = SymbolNew("m_0", units.mass)
"""
**Rest mass**, also called **invariant mass**, **intrinsic mass**, or **proper mass**, is the portion of
the total mass of an object or a system of objects that is the same in all references frames related by
Lorentz transformations.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Invariant_mass#>`__.
"""

proper_length = SymbolNew("l_0", units.length)
"""
**Rest length**, or **proper length**, is the length of an object in the object's rest frame.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Proper_length#>`__.
"""

relativistic_mass = SymbolNew("m", units.mass)
"""
**Relativistic mass** of a moving object is its mass measured in an external inertial frame of reference.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mass_in_special_relativity#>`__.
"""
