"""
Electrodynamics (Symbols)
=========================

Symbols related to electrodynamics.
"""

from sympy.physics import units
from symplyphysics.core.dimensions import dimensionless
from symplyphysics.core.symbols.symbols import SymbolNew

admittance = SymbolNew("Y", units.conductance)
"""
**Admittance** is a measure of how easily a circuit or device will allow a current to flow, defined as the reciprocal
of impedance.
"""

electrical_conductance = SymbolNew("G", units.conductance)
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

magnetic_flux = SymbolNew("Phi_B", units.magnetic_flux, display_latex="\\Phi_\\mathbf{B}")
r"""
**Magnetic flux** through a surface is the surface integral of the normal component of the magnetic field :math:`\mathbf{B}`
over that surface.
"""

absolute_permittivity = SymbolNew("epsilon",
    units.capacitance / units.length,
    display_latex="\\varepsilon")
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

current = SymbolNew("I", units.current)
"""
**Current** is a flow of charged particles moving through an electrical conductor or space.
"""

electrical_resistance = SymbolNew("R", units.impedance)
"""
**Resistance** is the measure of the degree to which a conductor opposes an electric current through that conductor.
It is the real part of the complex-valued impedance.
"""

electric_dipole_moment = SymbolNew("p", units.charge * units.length)
"""
**Electric dipole moment** is a measure of the separation of positive and negative electrical charges within a system:
that is, a measure of the system's overall polarity.
"""

electric_field_strength = SymbolNew("E", units.voltage / units.length)
"""
**Electric field strength** refers to the magnitude of the electric field.
"""

surface_charge_density = SymbolNew("sigma", units.charge / units.area, display_latex="\\sigma")
"""
**Surface charge density** is charge per unit surface area.
"""

electric_flux = SymbolNew("Phi_E", units.voltage * units.length, display_latex="\\Phi_\\mathbf{E}")
r"""
**Electric flux** through a surface is the surface integral of the normal component of the electric field :math:`\mathbf{E}`
over that surface.
"""

magnetic_flux_density = SymbolNew("B", units.magnetic_flux_density)
"""
**Magnetic flux density**, also called **magnetic induction**, is a physical quantity that predicts the force on a charged
particle in the Lorentz force law.
"""

electric_potential = SymbolNew("U_E", units.voltage, display_latex="U_\\mathbf{E}")
"""
**Electric potential** is defined as the amount of work or energy needed per unit of electric charge to move the charge from
a reference point to a specific point in an electric field.
"""

power_factor = SymbolNew("pf", dimensionless, display_latex="\\mathrm{pf}")
"""
**Power factor** of an AC power system is defined as the ratio of the real power absorbed by the load to the apparent power
flowing in the circuit.
"""

electrical_resistivity = SymbolNew("rho", units.impedance * units.length, display_latex="\\rho")
"""
**Electrical resistivity** is a fundamental specific property of a material that measures its electrical resistance or how
strongly it resists electric current.
"""

inductance = SymbolNew("L", units.inductance)
"""
**Inductance** is the tendency of an electrical conductor to oppose a change in the electric current flowing through it.
"""

electric_time_constant = SymbolNew("tau", units.time, display_latex="\\tau")
"""
**Time constant** is the parameter characterizing the response to a step input of a first-order, linear time-invariant
system. It is related to the speed of the response.
"""

electrical_reactance = SymbolNew("X", units.impedance)
"""
**Reactance** is the opposition presented to alternating current by inductance and capacitance.
"""

current_density = SymbolNew("j", units.current / units.area)
"""
**Current density** is the amount of charge per unit time that flows through a unit area of a chosen cross section.
"""
