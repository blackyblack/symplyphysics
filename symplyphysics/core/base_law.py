from sympy import solve, Eq

from symplyphysics import Symbol, Quantity, validate_input, validate_output
from symplyphysics.core.symbols.fraction import Fraction
from symplyphysics.core.symbols.probability import Probability


class BaseLaw:
    def __init__(self, law: Eq) -> None:
        self.law = law
        self.symbols = law.atoms(Symbol)

        variables_keys = list(str(symbol) for symbol in self.symbols)
        for key, symbol in zip(variables_keys, self.symbols):
            setattr(self, key[:-1], symbol)     # key[:-1] for deleting index in end of string

    def calculate_symbol_value(
            self,
            variable: Symbol,
            dict_quantities: dict[Symbol, Quantity | Fraction | Probability]
        ) -> Quantity | Fraction | Probability:
        if variable in self.symbols:
            pass
        else:
            raise ValueError("The specified variable is not in the law")

        symbols_set = self.symbols.copy()
        symbols_set.remove(variable)

        variables_list = list(symbols_set)

        variables_for_calculate = list(str(symbol) for symbol in variables_list)
        dimensions = list(symbol.dimension for symbol in variables_list)

        dict_validate = dict(
            zip(variables_for_calculate, dimensions)
        )

        @validate_input(**dict_validate)
        @validate_output(variable.dimension)
        def __calculate_variable(*quantities_: Quantity | Fraction | Probability) -> Quantity | Fraction | Probability:
            solved = solve(self.law, variable, dict=True)[0][variable]
            result_expr = solved.subs(zip(variables_list, quantities_))
            result = Quantity(result_expr)
            return result

        return __calculate_variable(*dict_quantities.values())
