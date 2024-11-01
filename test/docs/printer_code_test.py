from sympy import Integral, evaluate, exp, pi, sin, sqrt
from symplyphysics import SymbolNew, clone_as_symbol, quantities, symbols, clone_as_function
from symplyphysics.docs.printer_code import code_str


def test_add() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a + b
    assert code_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a + b + c
    assert code_str(expr) == "a + b + c"


def test_mul() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * b
    assert code_str(expr) == "a * b"


def test_div() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a / b
    assert code_str(expr) == "a / b"


def test_div_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = a / b + c / d
    assert code_str(expr) == "a / b + c / d"


def test_div_braces_common_frac() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert code_str(expr) == "(a + b) / (c + d)"


def test_sub() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a - b
    assert code_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a - b + c
    assert code_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * (b + a)
    assert code_str(expr) == "a * (b + a)"
    with evaluate(False):
        expr = (a + b) * a
    assert code_str(expr) == "(a + b) * a"


def test_mul_with_double_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) * (c + d)
    assert code_str(expr) == "(a + b) * (c + d)"


def test_mul_with_sqrt_and_quantity() -> None:
    mass = symbols.mass
    temperature = symbols.temperature
    with evaluate(False):
        expr = sqrt(2 * pi / (mass * quantities.boltzmann_constant * temperature))
    assert code_str(expr) == "sqrt(2 * pi / (m * k_B * T))"


def test_mul_minus_one() -> None:
    energy_of_state = symbols.energy
    temperature = symbols.temperature
    with evaluate(False):
        expr = exp(-1 * energy_of_state / (quantities.boltzmann_constant * temperature))
    assert code_str(expr) == "exp(-E / (k_B * T))"


def test_mul_integral() -> None:
    time = symbols.time
    time_before = SymbolNew("t_0", display_latex="t_0")
    time_after = SymbolNew("t_1", display_latex="t_1")
    force = clone_as_function(symbols.force, display_symbol="F(t)")
    with evaluate(False):
        expr = Integral(force(time), (time, time_before, time_after))
    assert code_str(expr) == "Integral(F(t), (t, t_0, t_1))"


def test_mul_sqrt() -> None:
    speed = symbols.speed
    with evaluate(False):
        expr = 1 / sqrt(1 - speed**2 / quantities.speed_of_light**2)
    assert code_str(expr) == "1 / sqrt(1 - v^2 / c^2)"


def test_mul_fraction_and_sin() -> None:
    mass = symbols.mass
    driving_force_amplitude = symbols.force
    natural_angular_frequency = SymbolNew("w_0", display_latex="\\omega_0")
    time = symbols.time
    with evaluate(False):
        expr = (driving_force_amplitude /
            (2 * mass * natural_angular_frequency)) * time * sin(natural_angular_frequency * time)
    assert code_str(expr) == "F / (2 * m * w_0) * t * sin(w_0 * t)"


def test_fraction_mul_inner_fraction() -> None:
    first_charge = clone_as_symbol(symbols.charge, display_symbol="q_1", display_latex="q_1")
    second_charge = clone_as_symbol(symbols.charge, display_symbol="q_2", display_latex="q_2")
    distance = symbols.distance
    with evaluate(False):
        expr = first_charge * second_charge / (4 * pi *
            quantities.vacuum_permittivity) / distance**2
    assert code_str(expr) == "q_1 * q_2 / (4 * pi * epsilon_0) / d^2"


def test_fraction_mul_inner_fraction_with_braces() -> None:
    first_charge = clone_as_symbol(symbols.charge, display_symbol="q_1", display_latex="q_1")
    second_charge = clone_as_symbol(symbols.charge, display_symbol="q_2", display_latex="q_2")
    distance = symbols.distance
    with evaluate(False):
        expr = 1 / (4 * pi *
            quantities.vacuum_permittivity) * (first_charge * second_charge / distance**2)
    assert code_str(expr) == "q_1 * q_2 / d^2 / (4 * pi * epsilon_0)"
