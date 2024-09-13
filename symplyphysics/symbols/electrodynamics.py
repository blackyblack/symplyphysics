"""
Electrodynamics
===============

Symbols related to electrodynamics.
"""

from symplyphysics import units, SymbolNew, dimensionless

admittance = SymbolNew("Y", units.conductance)
"""
**Admittance** is a measure of how easily a circuit or device will allow a current to flow, defined as the reciprocal
of impedance.
"""

conductance = SymbolNew("G", units.conductance)
"""
**Conductance** is the ability of charge to flow in a certain path. It is the reciprocal of electrical resistance.
"""

susceptance = SymbolNew("B", units.conductance)
"""
**Susceptance** is the imaginary part of the electrical admittance.
"""

electrical_impedance = SymbolNew("Z", units.impedance)
"""
**Electrical impedance** is the opposition to current presented by the combined effect of
resistance and reactance in a circuit.
"""

electromotive_force = SymbolNew("E", units.voltage, display_latex="\\mathcal{E}")
"""
**Electromotive force**, also **electromotance**, abbreviated **emf**, an energy transfer to an electric circuit per unit
of electric charge, measured in volts.
"""

magnetic_flux = SymbolNew("Phi_B", units.magnetic_flux, display_latex="\\Phi_B")
"""
**Magnetic flux** through a surface is the surface integral of the normal component of the magnetic field :math:`\mathbf{B}`
over that surface.
"""

absolute_permittivity = SymbolNew("epsilon", units.capacitance / units.length, display_latex="\\varepsilon")
"""
**Absolute permittivity**, or often sometimes **permittivity**, is a measure of the electric polarizability of a dielectric material.
"""

relative_permittivity = SymbolNew("epsilon_r", dimensionless, display_latex="\\varepsilon_r")
"""
**Relative permittivity** is the permittivity of a medium relative to that of free space.
Also see :attr:`~symplyphysics.quantities.vacuum_permittivity`.
"""

capacitance = SymbolNew("C", units.capacitance)
"""
**Capacitance** is the capacity of a material object or device to store electric charge.
"""

charge = SymbolNew("q", units.charge)
"""
**Electric charge** is the physical property of matter that causes it to experience a force when placed in an electromagnetic field.
It can be positive or negative. Like charges repel each other and unlike charges attract each other.
"""

voltage = SymbolNew("V", units.voltage)
"""
**Voltage** is the difference in electric potential between two points.
"""
