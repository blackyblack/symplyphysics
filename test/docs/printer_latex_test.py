from collections import namedtuple
from pytest import fixture
from sympy import (I, Integral, evaluate, exp, log, pi, sin, sqrt, Matrix, ImmutableMatrix, Symbol
    as SymSymbol)
from symplyphysics import (Quantity, Symbol, clone_as_function, clone_as_symbol, units,
    IndexedSymbol, IndexedSum, IndexedProduct, global_index, Function)
from symplyphysics.docs.printer_latex import latex_str
from symplyphysics.core.operations import symbolic

from symplyphysics.core.experimental.vectors import (VectorSymbol, VectorNorm, VectorDot,
    VectorCross, VectorMixedProduct, VectorFunction, vector_diff)
from symplyphysics.core.experimental.coordinate_systems import (CoordinateScalar, CoordinateVector,
    CartesianCoordinateSystem, CylindricalCoordinateSystem, QuantityCoordinateVector)

Args = namedtuple(
    "Args",
    [
    "mass",
    "temperature",
    "boltzmann_constant",
    "energy",
    "time",
    "force",
    "speed",
    "speed_of_light",
    "charge",
    "distance",
    "vacuum_permittivity",
    "electric_dipole_moment",
    "intensity",
    ],
)


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = Symbol("m")
    temperature = Symbol("T")
    boltzmann_constant = Quantity(display_latex="k_\\text{B}")
    energy = Symbol("E")
    time = Symbol("t")
    force = Symbol("F")
    speed = Symbol("v")
    speed_of_light = Quantity(display_symbol="c")
    charge = Symbol("q")
    distance = Symbol("d")
    vacuum_permittivity = Quantity(display_latex="\\varepsilon_0")
    electric_dipole_moment = Symbol("p")
    intensity = Symbol("I")
    return Args(
        mass=mass,
        temperature=temperature,
        boltzmann_constant=boltzmann_constant,
        energy=energy,
        time=time,
        force=force,
        speed=speed,
        speed_of_light=speed_of_light,
        charge=charge,
        distance=distance,
        vacuum_permittivity=vacuum_permittivity,
        electric_dipole_moment=electric_dipole_moment,
        intensity=intensity,
    )


def test_add() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a + b
    assert latex_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a + b + c
    assert latex_str(expr) == "a + b + c"


def test_mul() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * b
    assert latex_str(expr) == "a b"


def test_div() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a / b
    assert latex_str(expr) == "\\frac{a}{b}"


def test_div_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = a / b + c / d
    assert latex_str(expr) == "\\frac{a}{b} + \\frac{c}{d}"


def test_div_no_braces_common_frac() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert latex_str(expr) == "\\frac{a + b}{c + d}"


def test_sub() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a - b
    assert latex_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a - b + c
    assert latex_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * (b + a)
    assert latex_str(expr) == "a \\left(b + a\\right)"
    with evaluate(False):
        expr = (a + b) * a
    assert latex_str(expr) == "\\left(a + b\\right) a"


def test_mul_with_double_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = (a + b) * (c + d)
    assert latex_str(expr) == "\\left(a + b\\right) \\left(c + d\\right)"


def test_mul_with_sqrt_and_quantity(test_args: Args) -> None:
    with evaluate(False):
        expr = sqrt(2 * pi /
            (test_args.mass * test_args.boltzmann_constant * test_args.temperature))
    assert latex_str(expr) == "\\sqrt{\\frac{2 \\pi}{m k_\\text{B} T}}"


def test_mul_minus_one(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * test_args.energy / (test_args.boltzmann_constant * test_args.temperature))
    assert latex_str(expr) == "\\exp{\\left(- \\frac{E}{k_\\text{B} T} \\right)}"


def test_mul_integral(test_args: Args) -> None:
    time = test_args.time
    time_before = clone_as_symbol(time, subscript="0")
    time_after = clone_as_symbol(time, subscript="1")
    force = clone_as_function(test_args.force)
    with evaluate(False):
        expr = Integral(force(time), (time, time_before, time_after))
    assert latex_str(expr) == "\\int\\limits_{t_{0}}^{t_{1}} F{\\left(t \\right)}\\, dt"


def test_mul_sqrt(test_args: Args) -> None:
    with evaluate(False):
        expr = 1 / sqrt(1 - test_args.speed**2 / test_args.speed_of_light**2)
    assert latex_str(expr) == "\\frac{1}{\\sqrt{1 - \\frac{v^{2}}{c^{2}}}}"


def test_mul_fraction_and_sin(test_args: Args) -> None:
    natural_angular_frequency = Symbol("w_0", display_latex="\\omega_0")
    time = test_args.time
    with evaluate(False):
        expr = (test_args.force / (2 * test_args.mass * natural_angular_frequency)) * time * sin(
            natural_angular_frequency * time)
    assert latex_str(expr) == "\\frac{F}{2 m \\omega_{0}} t \\sin{\\left(\\omega_{0} t \\right)}"


def test_fraction_mul_inner_fraction(test_args: Args) -> None:
    first_charge = clone_as_symbol(test_args.charge, subscript="1")
    second_charge = clone_as_symbol(test_args.charge, subscript="2")
    with evaluate(False):
        expr = first_charge * second_charge / (4 * pi *
            test_args.vacuum_permittivity) / test_args.distance**2
    assert latex_str(expr) == "\\frac{q_{1} q_{2}}{4 \\pi \\varepsilon_0} \\frac{1}{d^{2}}"


def test_fraction_mul_fraction(test_args: Args) -> None:
    with evaluate(False):
        expr = 1 / (2 * pi * test_args.vacuum_permittivity) * (test_args.electric_dipole_moment /
            test_args.distance**3)
    assert latex_str(expr) == "\\frac{1}{2 \\pi \\varepsilon_0} \\frac{p}{d^{3}}"


def test_log10(test_args: Args) -> None:
    reference_intensity = Quantity(1e-12 * units.watt / units.meter**2, display_symbol="I_0")
    with evaluate(False):
        expr = log(test_args.intensity / reference_intensity, 10)
    assert latex_str(expr) == "\\log_{10} \\left( \\frac{I}{I_0} \\right)"


def test_minus_imaginary(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * (I / test_args.boltzmann_constant))
    assert latex_str(expr) == "\\exp{\\left(- \\frac{i}{k_\\text{B}} \\right)}"


def test_log_squared(test_args: Args) -> None:
    with evaluate(False):
        expr = log(test_args.boltzmann_constant)**2
    assert latex_str(expr) == "\\log \\left( k_\\text{B} \\right)^{2}"


def test_mutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = Matrix([[a, b], [c, d]])
    assert latex_str(expr) == "\\begin{pmatrix} a & b \\\\ c & d \\end{pmatrix}"

    # 4-row vector
    with evaluate(False):
        expr = Matrix([a, b, c, d])
    assert latex_str(expr) == "\\begin{pmatrix} a \\\\ b \\\\ c \\\\ d \\end{pmatrix}"

    # 4-col vector
    with evaluate(False):
        expr = Matrix([a, b, c, d]).T
    assert latex_str(expr) == "\\begin{pmatrix} a & b & c & d \\end{pmatrix}"


def test_immutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = ImmutableMatrix([[a, b], [c, d]])
    assert latex_str(expr) == "\\begin{pmatrix} a & b \\\\ c & d \\end{pmatrix}"

    # 4-row vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d])
    assert latex_str(expr) == "\\begin{pmatrix} a \\\\ b \\\\ c \\\\ d \\end{pmatrix}"

    # 4-col vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d]).T
    assert latex_str(expr) == "\\begin{pmatrix} a & b & c & d \\end{pmatrix}"


def test_indexed_symbol() -> None:
    f = IndexedSymbol("F")
    expr = f[global_index]
    assert latex_str(expr) == "{F}_{i}"


def test_indexed_sum() -> None:
    f = IndexedSymbol("F")
    expr = IndexedSum(f[global_index], global_index)
    assert latex_str(expr) == "\\sum_i {F}_{i}"

    with evaluate(False):
        expr = IndexedSum(f[global_index]**2, global_index)
    assert latex_str(expr) == "\\sum_i {F}_{i}^{2}"


def test_indexed_product() -> None:
    f = IndexedSymbol("F")

    expr = IndexedProduct(f[global_index], global_index)
    assert latex_str(expr) == "\\prod_i {F}_{i}"

    with evaluate(False):
        expr = IndexedProduct(f[global_index]**2, global_index)
    assert latex_str(expr) == "\\prod_i {F}_{i}^{2}"


def test_average() -> None:
    a = Symbol("a")
    expr = symbolic.Average(a)
    assert latex_str(expr) == "\\langle a \\rangle"


def test_finite_difference() -> None:
    a = Symbol("a")

    expr = symbolic.FiniteDifference(a)
    assert latex_str(expr) == "\\Delta a"

    expr = symbolic.FiniteDifference(a, wrap_latex=True)
    assert latex_str(expr) == "\\Delta \\left( a \\right)"


def test_exact_differential() -> None:
    a = Symbol("a")

    expr = symbolic.ExactDifferential(a)
    assert latex_str(expr) == "d a"

    expr = symbolic.ExactDifferential(a, wrap_latex=True)
    assert latex_str(expr) == "d \\left( a \\right)"


def test_quantity() -> None:
    expr = Quantity(0)
    assert latex_str(expr) == "0"

    expr = Quantity(4 * units.meter)
    assert latex_str(expr) == "4.0 m"

    expr = units.acceleration_due_to_gravity
    assert latex_str(expr) == "\\text{g}"


def test_vector_symbol() -> None:
    expr = VectorSymbol("v")
    assert latex_str(expr) == "{\\vec v}"

    expr = VectorSymbol("v", units.velocity)
    assert latex_str(expr) == "{\\vec v}"

    expr = VectorSymbol("w", display_latex="\\omega")
    assert latex_str(expr) == "\\omega"

    expr = VectorSymbol("v", units.velocity, display_latex="\\omega")
    assert latex_str(expr) == "\\omega"


def test_vector_combination(test_args: Args) -> None:
    v = VectorSymbol("v", units.velocity)
    q = VectorSymbol("q", display_latex="\\hat{q}")

    with evaluate(False):
        expr = test_args.speed_of_light * q + test_args.mass * v
    assert latex_str(expr) == "c \\hat{q} + m {\\vec v}"


def test_vector_norm() -> None:
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    expr = VectorNorm(v)
    assert latex_str(expr) == "\\left \\Vert {\\vec v} \\right \\Vert"

    with evaluate(False):
        expr = VectorNorm(v + 2 * w)
    assert latex_str(expr) == "\\left \\Vert {\\vec v} + 2 {\\vec w} \\right \\Vert"


def test_vector_dot() -> None:
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    with evaluate(False):
        expr = VectorDot(v, w)

    assert latex_str(expr) == "\\left( {\\vec v}, {\\vec w} \\right)"


def test_vector_cross() -> None:
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    with evaluate(False):
        expr = VectorCross(v, w)
    assert latex_str(expr) == "\\left[ {\\vec v}, {\\vec w} \\right]"


def test_vector_mixed_product() -> None:
    u = VectorSymbol("u")
    v = VectorSymbol("v")
    w = VectorSymbol("w")

    with evaluate(False):
        expr = VectorMixedProduct(u, v, w)
    assert latex_str(expr) == "\\left( {\\vec u}, {\\vec v}, {\\vec w} \\right)"


def test_applied_vector_function() -> None:
    x = Symbol("x")
    y = Symbol("y")
    f = VectorFunction("F")

    expr = f(x)
    assert latex_str(expr) == "{\\vec F} \\left( x \\right)"

    expr = f(x, y)
    assert latex_str(expr) == "{\\vec F} \\left( x, y \\right)"

    f = VectorFunction("F", dimension=units.force)

    expr = f(x)
    assert latex_str(expr) == "{\\vec F} \\left( x \\right)"

    expr = f(x, y)
    assert latex_str(expr) == "{\\vec F} \\left( x, y \\right)"


def test_vector_function() -> None:
    x = Symbol("x")
    y = Symbol("y")

    expr = VectorFunction("F", [x])
    assert latex_str(expr) == "{\\vec F} \\left( x \\right)"

    expr = VectorFunction("F", [x, y])
    assert latex_str(expr) == "{\\vec F} \\left( x, y \\right)"

    expr = VectorFunction("F", [x], dimension=units.force)
    assert latex_str(expr) == "{\\vec F} \\left( x \\right)"

    expr = VectorFunction("F", [x, y], dimension=units.force)
    assert latex_str(expr) == "{\\vec F} \\left( x, y \\right)"

    expr = VectorFunction("w", [x], display_latex="\\omega")
    assert latex_str(expr) == "\\omega \\left( x \\right)"

    expr = VectorFunction("w", [x, y], display_latex="\\omega")
    assert latex_str(expr) == "\\omega \\left( x, y \\right)"


def test_vector_derivative() -> None:
    x = Symbol("x")
    y = Symbol("y")
    f = VectorFunction("F", [x, y])

    expr = vector_diff(f(x, y), x)
    assert latex_str(expr) == "\\frac{\\partial}{\\partial x} {\\vec F} \\left( x, y \\right)"

    f = VectorFunction("F", [x, y], dimension=units.force)
    expr = vector_diff(f(x, y), x)
    assert latex_str(expr) == "\\frac{\\partial}{\\partial x} {\\vec F} \\left( x, y \\right)"


def test_cartesian_coordinate_scalar() -> None:
    f = Function("F")

    sys = CartesianCoordinateSystem()
    x, y, z = sys.base_scalars

    with evaluate(False):
        scalar = x + y + z

    expr = CoordinateScalar(scalar, sys)
    assert latex_str(expr) == latex_str(scalar)

    with evaluate(False):
        scalar = f(x, y, z)

    expr = CoordinateScalar(scalar, sys)
    assert latex_str(expr) == latex_str(scalar)


def test_non_cartesian_coordinate_scalar() -> None:
    f = Function("F")

    sys = CylindricalCoordinateSystem()
    rho, phi, z = sys.base_scalars
    p = SymSymbol("P")

    with evaluate(False):
        scalar = rho + z

    expr = CoordinateScalar(scalar, sys, p)
    assert latex_str(expr) == latex_str(scalar)

    with evaluate(False):
        scalar = f(rho, phi, z)

    expr = CoordinateScalar(scalar, sys, p)
    assert latex_str(expr) == latex_str(scalar)


def test_cartesian_coordinate_vector() -> None:
    sys = CartesianCoordinateSystem()
    x, y, z = sys.base_scalars

    with evaluate(False):
        components = ImmutableMatrix([x, y, z])

    expr = CoordinateVector(components, sys)
    assert latex_str(expr) == latex_str(components)

    components = ImmutableMatrix([Quantity(4 * units.meter), Quantity(-3 * units.meter), 0])
    expr = QuantityCoordinateVector(components, sys)
    assert latex_str(expr) == latex_str(components)


def test_non_cartesian_coordinate_vector() -> None:
    sys = CylindricalCoordinateSystem()
    p = SymSymbol("P")
    rho, _, z = sys.base_scalars

    with evaluate(False):
        components = ImmutableMatrix([rho, 0, z])

    expr = CoordinateVector(components, sys, p)
    assert latex_str(expr) == latex_str(components)

    components = ImmutableMatrix([Quantity(4 * units.meter), Quantity(-3 * units.meter), 0])
    expr = QuantityCoordinateVector(components, sys, p)
    assert latex_str(expr) == latex_str(components)
