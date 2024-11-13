from sympy import Integral, evaluate, exp, pi, sin, sqrt
from symplyphysics import SymbolNew, quantities, symbols, clone_as_function, clone_as_symbol
from symplyphysics.core.operations.average import Average
from symplyphysics.docs.printer_latex import latex_str


def test_add() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a + b
    assert latex_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a + b + c
    assert latex_str(expr) == "a + b + c"


def test_mul() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * b
    assert latex_str(expr) == "a b"


def test_div() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a / b
    assert latex_str(expr) == "\\frac{a}{b}"


def test_div_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = a / b + c / d
    assert latex_str(expr) == "\\frac{a}{b} + \\frac{c}{d}"


def test_div_no_braces_common_frac() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert latex_str(expr) == "\\frac{a + b}{c + d}"


def test_sub() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a - b
    assert latex_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a - b + c
    assert latex_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * (b + a)
    assert latex_str(expr) == "a \\left(b + a\\right)"
    with evaluate(False):
        expr = (a + b) * a
    assert latex_str(expr) == "\\left(a + b\\right) a"


def test_mul_with_double_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) * (c + d)
    assert latex_str(expr) == "\\left(a + b\\right) \\left(c + d\\right)"


def test_mul_with_sqrt_and_quantity() -> None:
    mass = symbols.mass
    temperature = symbols.temperature
    with evaluate(False):
        expr = sqrt(2 * pi / (mass * quantities.boltzmann_constant * temperature))
    assert latex_str(expr) == "\\sqrt{\\frac{2 \\pi}{m k_\\text{B} T}}"


def test_mul_minus_one() -> None:
    energy_of_state = symbols.energy
    temperature = symbols.temperature
    with evaluate(False):
        expr = exp(-1 * energy_of_state / (quantities.boltzmann_constant * temperature))
    assert latex_str(expr) == "\\exp{\\left(- \\frac{E}{k_\\text{B} T} \\right)}"


def test_mul_integral() -> None:
    time = symbols.time
    time_before = SymbolNew("t_0", display_latex="t_0")
    time_after = SymbolNew("t_1", display_latex="t_1")
    force = clone_as_function(symbols.force, display_symbol="F(t)")
    with evaluate(False):
        expr = Integral(force(time), (time, time_before, time_after))
    assert latex_str(expr) == "\\int\\limits_{t_{0}}^{t_{1}} F{\\left(t \\right)}\\, dt"


def test_mul_sqrt() -> None:
    speed = symbols.speed
    with evaluate(False):
        expr = 1 / sqrt(1 - speed**2 / quantities.speed_of_light**2)
    assert latex_str(expr) == "\\frac{1}{\\sqrt{1 - \\frac{v^{2}}{c^{2}}}}"


def test_mul_fraction_and_sin() -> None:
    mass = symbols.mass
    driving_force_amplitude = symbols.force
    natural_angular_frequency = SymbolNew("w_0", display_latex="\\omega_0")
    time = symbols.time
    with evaluate(False):
        expr = (driving_force_amplitude /
            (2 * mass * natural_angular_frequency)) * time * sin(natural_angular_frequency * time)
    assert latex_str(expr) == "\\frac{F}{2 m \\omega_{0}} t \\sin{\\left(\\omega_{0} t \\right)}"


def test_fraction_mul_inner_fraction() -> None:
    first_charge = clone_as_symbol(symbols.charge, display_symbol="q_1", display_latex="q_1")
    second_charge = clone_as_symbol(symbols.charge, display_symbol="q_2", display_latex="q_2")
    distance = symbols.distance
    with evaluate(False):
        expr = first_charge * second_charge / (4 * pi *
            quantities.vacuum_permittivity) / distance**2
    assert latex_str(expr) == "\\frac{q_{1} q_{2}}{4 \\pi \\varepsilon_0} \\frac{1}{d^{2}}"


def test_fraction_mul_fraction() -> None:
    electric_dipole_moment = symbols.electric_dipole_moment
    distance = symbols.distance
    with evaluate(False):
        expr = 1 / (2 * pi * quantities.vacuum_permittivity) * (electric_dipole_moment /
            distance**3)
    assert latex_str(expr) == "\\frac{1}{2 \\pi \\varepsilon_0} \\frac{p}{d^{3}}"


def test_average_operation() -> None:
    v = clone_as_symbol(symbols.speed, display_symbol="v", display_latex="v")
    t = clone_as_symbol(symbols.time, display_symbol="t", display_latex="t")
    avg = Average(v * t)
    assert latex_str(avg) == "\\langle v t \\rangle"
