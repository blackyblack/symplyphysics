r"""
**Van der Waals Model**
=======================

**Van der Waals equation of state** is a model which extends the ideal gas law to include
the non-zero size of gas molecules and the interactions between them.

.. _vdw_critical_parameters_def:

**Critical parameters** of the van der Waals equation of state are such values of volume, pressure,
and temperature at which the isotherm has an inflection point whose tangent at that point is zero,
i.e. the first and second derivatives of pressure with respect to volume at constant temperature are
zero:

.. math::
    \left( \frac{\partial p}{\partial V} \right)_T = \left( \frac{\partial^2 p}{\partial V^2} \right)_T = 0

.. _vdw_reduced_units_def:

**Reduced units** are used in the :ref:`dimensionless <Dimensionless van der Waals equation>`
van der Waals equation of state. A reduced quantity :math:`X_r` is defined as the ratio of quantity
:math:`X` to the corresponding :ref:`critical parameter <vdw_critical_parameters_def>` :math:`X_\text{c}`.

**Links:**

#. `Wikipedia — Van der Waals equation <https://en.wikipedia.org/wiki/Van_der_Waals_equation>`__
"""
