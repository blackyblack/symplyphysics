#!/usr/bin/env python3
"""
Light enters a prism whose cross section is an equilateral triangle. The prism is made of fused
silica. Plot the dispersion of light after it passes the prism.
"""

from sympy import symbols as sym_symbols, Eq, solve, pi, plot, sin, cos, tan, Expr, Rational, sqrt
from symplyphysics import symbols, clone_as_symbol
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.laws.optics import refraction_angle_from_environments as snells_law

wavelength_in_um = sym_symbols("lambda_um", positive=True)

air_refractive_index = clone_as_symbol(
    symbols.relative_refractive_index,
    subscript="0",
    positive=True,
)

glass_refractive_index = clone_as_symbol(
    symbols.relative_refractive_index,
    subscript="1",
    positive=True,
)

## Part 1. Solve the problem analytically.

# TODO: create a figure to refer to.

lb = sym_symbols("LB", positive=True)

# Length of the triangle sides
a = sym_symbols("a", positive=True)

incoming_incidence_angle = clone_as_symbol(symbols.angle, positive=True)

refraction_angle_expr = solve(snells_law.law, snells_law.refraction_angle)[1]

incoming_refraction_angle_expr = refraction_angle_expr.subs({
    snells_law.incidence_refractive_index: air_refractive_index,
    snells_law.incidence_angle: incoming_incidence_angle,
    snells_law.resulting_refractive_index: glass_refractive_index,
})

lbc_expr = pi / 2 - incoming_refraction_angle_expr

lcb_expr = pi - pi / 3 - lbc_expr

outgoing_incidence_angle_expr = pi / 2 - lcb_expr

outgoing_refraction_angle_expr = refraction_angle_expr.subs({
    snells_law.incidence_refractive_index: glass_refractive_index,
    snells_law.incidence_angle: outgoing_incidence_angle_expr,
    snells_law.resulting_refractive_index: air_refractive_index,
})

lcd_expr = pi / 2 + outgoing_refraction_angle_expr

cdn_expr = lcd_expr - 2 * pi / 3

cdx_expr = pi - cdn_expr

lc = sym_symbols("LC", positive=True)

# Law of sines
sines_eqn = Eq(lb / sin(lcb_expr), lc / sin(lbc_expr))

lc_expr = solve(sines_eqn, lc)[0]

cn_expr = a - lc_expr

ch_expr = sin(pi / 3) * cn_expr

hn_expr = cos(pi / 3) * cn_expr

oh_expr = a / 2 - hn_expr

hcd_expr = pi / 2 - cdn_expr

hd_expr = ch_expr * tan(hcd_expr)

od_expr = (oh_expr + hd_expr).simplify()

## Part 2. Program the dependency of refractive index of fused silica on wavelength.
# Source: `https://doi.org/10.1364/JOSA.55.001205`

a1, a2, a3, a4, a5, a6 = sym_symbols("a_1:7", positive=True)

a_values = {
    a1: 0.6961663,
    a2: 0.0684043,
    a3: 0.4079426,
    a4: 0.1162414,
    a5: 0.8974794,
    a6: 0.9896161e1,
}

glass_refractive_index_eqn = Eq(
    glass_refractive_index**2 - 1,
    a1 * wavelength_in_um**2 / (wavelength_in_um**2 - a2**2) + a3 * wavelength_in_um**2 /
    (wavelength_in_um**2 - a4**2) + a5 * wavelength_in_um**2 / (wavelength_in_um**2 - a6**2),
)

# Choose the positive branch
(glass_refractive_index_expr,) = (sol
    for sol in solve(glass_refractive_index_eqn, glass_refractive_index)
    if not sol.could_extract_minus_sign())

glass_refractive_index_subs = glass_refractive_index_expr.subs(a_values)

## Part 3. Program the dependency of refractive index of air on wavelength.
# Source: `https://refractiveindex.info/?shelf=other&book=air&page=Ciddor`

b1, b2, b3, b4 = sym_symbols("b_1:5", positive=True)

b_values = {
    b1: 0.5792105e-1,
    b2: 0.2380185e+3,
    b3: 0.1679170e-2,
    b4: 0.5736200e+2,
}

air_refractive_index_expr = (1 + b1 / (b2 - wavelength_in_um**(-2)) + b3 /
    (b4 - wavelength_in_um**(-2)))

air_refractive_index_subs = air_refractive_index_expr.subs(b_values)

## Part 4. Plot the passage of light after the prism.

# Taken from `https://en.wikipedia.org/wiki/Visible_spectrum#Spectral_colors`
end_wavelengths_in_nm = [380e-3, 450e-3, 485e-3, 500e-3, 565e-3, 590e-3, 625e-3, 750e-3]
mid_wavelengths_in_nm = [
    (l1 + l2) / 2 for l1, l2 in zip(end_wavelengths_in_nm, end_wavelengths_in_nm[1:])
]
colors = "purple blue cyan green yellow orange red".split()
assert len(mid_wavelengths_in_nm) == len(colors)

x = sym_symbols("x", positive=True)

values = {
    a: 1,
    lb: 1 / 3,
    incoming_incidence_angle: pi / 6,
}

baseplot = plot(
    title="Refraction of light by triangular prism",
    xlabel="",
    ylabel="",
    xlim=(0, 1),
    ylim=(0, 1),
    size=(8, 8),
    legend=True,
    show=False,
)

triangle_edge = two_point_function(Point2D(0, sqrt(3) / 2), Point2D(Rational(1, 2), 0), x)
subplot = plot(triangle_edge, (x, 0, Rational(1, 2)), line_color="black", label="", show=False)
baseplot.extend(subplot)

for wavelength_in_um_, color in zip(mid_wavelengths_in_nm, colors):
    air_refractive_index_ = air_refractive_index_subs.subs(wavelength_in_um,
        wavelength_in_um_).evalf(5)

    glass_refractive_index_ = glass_refractive_index_subs.subs(wavelength_in_um,
        wavelength_in_um_).evalf(5)

    values_ = values | {
        air_refractive_index: air_refractive_index_,
        glass_refractive_index: glass_refractive_index_,
    }

    oh_ = oh_expr.subs(values_).n()
    ch_ = ch_expr.subs(values_).n()
    od_ = od_expr.subs(values_).n()

    ray = two_point_function(Point2D(oh_, ch_), Point2D(od_, 0), x)
    subplot = plot(
        ray,
        (x, oh_, 1),
        line_color=color,
        label=rf"$\lambda \approx {round(wavelength_in_um_, 3)} \, \mu\text{{m}}$",
        show=False,
    )
    baseplot.extend(subplot)

baseplot.show()
