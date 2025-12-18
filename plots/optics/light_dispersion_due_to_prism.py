#!/usr/bin/env python3
"""
Light enters a prism whose cross section is an equilateral triangle. The prism is made of fused
silica. Plot the dispersion of light after it passes the prism.
"""

from sympy import symbols as sym_symbols, Eq, solve, pi, plot, sin, cos, tan, Rational, sqrt
from symplyphysics import symbols, clone_as_symbol
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.laws.optics.geometrical_optics.refraction import refraction_angle_from_environments as snells_law

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

# NOTE Refer to the figure at `../../img/light_dispersion_due_to_prism.svg`

# Distance from the top vertex to the point of intersection of the light ray with the first side.
lb = sym_symbols("LB", positive=True)

# Length of the triangle sides
a = sym_symbols("a", positive=True)

# Angle `E1BA` in the figure.
incoming_incidence_angle = clone_as_symbol(symbols.angle, positive=True)

refraction_angle_expr = solve(snells_law.law, snells_law.refraction_angle)[1]

# Angle `E2BC` in the figure.
incoming_refraction_angle_expr = refraction_angle_expr.subs({
    snells_law.incidence_refractive_index: air_refractive_index,
    snells_law.incidence_angle: incoming_incidence_angle,
    snells_law.resulting_refractive_index: glass_refractive_index,
})

angle_lbe2_expr = pi / 2

angle_lbc_expr = angle_lbe2_expr - incoming_refraction_angle_expr

# The triangle `MLN` is equilateral
angle_blc_expr = pi / 3

# Per the angle sum property in the triangle `BLC`
angle_lcb_expr = pi - angle_blc_expr - angle_lbc_expr

angle_lcf2_expr = pi / 2

# Angle `BCF2` in the figure.
outgoing_incidence_angle_expr = angle_lcf2_expr - angle_lcb_expr

# Angle `F1CD` in the figure.
outgoing_refraction_angle_expr = refraction_angle_expr.subs({
    snells_law.incidence_refractive_index: glass_refractive_index,
    snells_law.incidence_angle: outgoing_incidence_angle_expr,
    snells_law.resulting_refractive_index: air_refractive_index,
})

angle_lcf1_expr = pi / 2

angle_lcd_expr = angle_lcf1_expr + outgoing_refraction_angle_expr

# The triangle `MLN` is equilateral.
angle_cnh_expr = pi / 3

# Angle outer to angle `CNH`
angle_cnd_expr = pi - angle_cnh_expr

# External angle of triangle: `angle(LCD) = angle(CND) + angle(CDN)`
angle_cdn_expr = angle_lcd_expr - angle_cnd_expr

# Angle `CDx` in the figure, i.e. the angle outer to angle `CDN`.
angle_cdx_expr = pi - angle_cdn_expr

lc = sym_symbols("LC", positive=True)

# Law of sines in the triangle `LBC` (see figure)
sines_eqn = Eq(lb / sin(angle_lcb_expr), lc / sin(angle_lbc_expr))

lc_expr = solve(sines_eqn, lc)[0]

cn_expr = a - lc_expr

# In the right triangle `CNH`
ch_expr = sin(angle_cnh_expr) * cn_expr

# In the right triangle `CNH`
hn_expr = cos(angle_cnh_expr) * cn_expr

# Holds for equilateral triangles
on_expr = a / 2

# `ON = OH + HN => OH = ON - HN`
oh_expr = on_expr - hn_expr

angle_chd_expr = pi / 2

# Per the angle sum property in the triangle `CHD`
angle_hcd_expr = pi - angle_chd_expr - angle_cdn_expr

# In the right triangle `CHD`
hd_expr = ch_expr * tan(angle_hcd_expr)

od_expr = (oh_expr + hd_expr).simplify()

## Part 2. Program the dependency of refractive index of fused silica on wavelength.
# Source: `https://doi.org/10.1364/JOSA.55.001205`

wavelength_in_um = sym_symbols("lambda_um", positive=True)

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
