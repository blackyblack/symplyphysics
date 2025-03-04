from __future__ import annotations

from collections import defaultdict
from typing import Any, Optional, Sequence

from sympy import Add, Atom, Basic, Expr, Mul, S, sympify
from sympy.core.evalf import EvalfMixin
from sympy.core.parameters import global_parameters
from sympy.physics.units import Dimension
from sympy.physics.units.systems.si import dimsys_SI
from sympy.printing.printer import Printer

from symplyphysics.core.dimensions import collect_expression_and_dimension
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.symbols.id_generator import next_id
from symplyphysics.core.symbols.symbols import DimensionSymbol


class VectorExpr(Basic, EvalfMixin):  # type: ignore[misc]
    """
    Base class for all vector expressions.
    """

    def doit(self, **_hints: Any) -> VectorExpr:
        return self

    def __add__(self, other: VectorExpr) -> VectorExpr:
        return VectorAdd(self, other).doit()

    def __sub__(self, other: VectorExpr) -> VectorExpr:
        return VectorAdd(self, -other).doit()

    def __mul__(self, other: Any) -> VectorExpr:
        return VectorScale(self, other).doit()

    def __neg__(self) -> VectorExpr:
        return VectorScale(self, -1).doit()

    def __pos__(self) -> VectorExpr:
        return self

    def __truediv__(self, other: Any) -> VectorExpr:
        # NOTE: probably need to check if `other` is not `0` to return a special "NaN" vector

        return VectorScale(self, 1 / other).doit()

    @property
    def is_zero(self) -> bool:
        # check with `isinstance` in case the user instanciates their own zero vector.
        return isinstance(self, _VectorZero)


class _VectorZero(VectorExpr):
    """
    Class expressing the notion of a zero vector. This class isn't intended to be instantiated
    except for the constant `ZERO`.
    """

    def _sympystr(self, _p: Printer) -> str:
        return "0"


ZERO = _VectorZero()


# NOTE: Instead of marking the `norm` of the `VectorSymbol` on the symbol itself, another
#       possibility is to create a separate class for `UnitVector`s. This way, when we add support
#       for vector components, we wouldn't need to check for the fact that the norm calculated
#       using the supplied components equals the norm given at the instantiation of the symbol.
#       But perhaps this simply gives us additional information about the vector and there's no
#       need to worry.
class VectorSymbol(DimensionSymbol, VectorExpr, Atom):  # type: ignore[misc]
    """
    Class representing a symbolic vector.

    Note that the norm of the vector cannot be `0`, instead, use the pre-defined constant `ZERO`
    for the zero vector.
    """

    _norm: Optional[Expr]
    _id: int

    is_symbol = True

    def __new__(
        cls,
        display_symbol: Optional[str] = None,
        dimension: Optional[Dimension] = None,
        *,
        norm: Optional[Any] = None,  # pylint: disable=redefined-outer-name
        display_latex: Optional[str] = None,
    ) -> VectorSymbol:
        return VectorExpr.__new__(cls)

    def __init__(
        self,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(1),
        *,
        norm: Optional[Any] = None,  # pylint: disable=redefined-outer-name
        display_latex: Optional[str] = None,
    ) -> None:
        self._id = next_id("VEC")

        if not display_symbol:
            display_symbol = f"VEC{self._id}"
            if not display_latex:
                display_latex = f"\\mathbf{{v}}_{{{self._id}}}"
        elif not display_latex:
            display_latex = f"\\mathbf{{{display_symbol}}}"

        DimensionSymbol.__init__(
            self,
            display_name=display_symbol,
            dimension=dimension,
            display_latex=display_latex,
        )

        if norm is not None:
            norm = sympify(norm)

            if not isinstance(norm, Expr):
                raise TypeError(f"Norm {norm} must be an Expr, got {type(norm).__name__}.")

            try:
                is_non_negative = bool(norm >= 0)
            except TypeError:
                # The case of more complex expressions that we can't compare.
                pass
            else:
                if not is_non_negative:
                    raise ValueError(f"Norm must be non-negative, got {norm}.")

            # TODO: move it into __new__?
            if norm == 0:
                raise ValueError("Use the constant ZERO for a zero vector.")

            norm_dim = collect_expression_and_dimension(norm)[1]
            if not dimsys_SI.equivalent_dims(norm_dim, self.dimension):
                raise UnitsError(f"The norm must be {self.dimension}, got {norm_dim}.")

        self._norm = norm

    @property
    def norm(self) -> Optional[Expr]:
        return self._norm

    def _hashable_content(self) -> tuple[Any, ...]:
        return (self._id,)


class VectorNorm(Expr):  # type: ignore[misc]

    @property
    def argument(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    def __init__(self, vector: VectorExpr) -> None:
        self._args = (vector,)

    def doit(self, **_hints: Any) -> Expr:
        vector = self.argument

        if vector.is_zero:
            return S.Zero

        if isinstance(vector, VectorSymbol) and vector.norm is not None:
            return vector.norm

        if isinstance(vector, VectorScale):
            return VectorNorm(vector.args[0]) * abs(vector.args[1])

        return self

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"norm({p.doprint(self.argument)})"

    def _eval_evalf(self, prec: int) -> VectorNorm:
        return VectorNorm(self.argument.evalf(n=prec)).doit()


def norm(vector: VectorExpr) -> Expr:
    return VectorNorm(vector).doit()


class VectorScale(VectorExpr):

    @property
    def vector(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def scale(self) -> Expr:
        return self.args[1]

    def __new__(cls, vector: VectorExpr, scale: Any, **kwargs: Any) -> VectorScale:
        return super().__new__(cls)  # type: ignore[no-any-return]

    def __init__(self, vector: VectorExpr, scale: Any, **kwargs: Any) -> None:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            scale = sympify(scale)
            if not isinstance(scale, Expr):
                raise TypeError(f"Scale {scale} must be an Expr, got {type(scale).__name__}.")

        self._args = (vector, scale)

    def doit(self, **_hints: Any) -> VectorExpr:
        vector = self.vector
        scale = self.scale

        while isinstance(vector, VectorScale):
            scale *= vector.scale
            vector = vector.vector

        if scale == 0:
            return ZERO

        if vector.is_zero or scale == 1:
            return vector

        if isinstance(vector, VectorAdd):
            addends = [VectorScale(addend, scale) for addend in vector.args]
            return VectorAdd(*addends)

        return VectorScale(vector, scale, evaluate=False)

    def _sympystr(self, p: Printer) -> str:
        vector, value = self.args

        if isinstance(value, (Add, Mul)):
            return f"{p.doprint(vector)}*({p.doprint(value)})"

        return f"{p.doprint(vector)}*{p.doprint(value)}"

    def _eval_evalf(self, prec: int) -> VectorScale:
        return VectorScale(self.vector.evalf(n=prec), self.scale.evalf(n=prec)).doit()


class Scale(Expr):  # type: ignore[misc]

    def __init__(self, scale: Any) -> None:
        self._args = (scale,)

    def __mul__(self, other: Any) -> VectorScale:
        if isinstance(other, VectorExpr):
            return VectorScale(other, self.args[0]).doit()

        raise TypeError(
            f"Scale can only be multiplied with a VectorExpr, got {type(other).__name__}.")

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass


class VectorAdd(VectorExpr):

    identity = ZERO

    def __init__(self, *vectors: VectorExpr) -> None:
        # TODO: add dimension check for the arguments

        # TODO: if `vectors` is empty, return `ZERO` (in __new__?)

        for vector in vectors:
            if not isinstance(vector, VectorExpr):
                raise TypeError(f"All addends must be VectorExpr, got {type(vector).__name__}")

        self._args = vectors

    def doit(self, **_hints: Any) -> VectorExpr:

        def flatten_additions(addends: Sequence[VectorExpr]) -> list[VectorExpr]:
            addends = list(addends)
            i = 0

            while i < len(addends):
                addend = addends[i]

                if addend.is_zero:
                    addends.pop(i)
                    continue

                if isinstance(addend, VectorAdd):
                    addends.pop(i)
                    addends.extend(addend.args)
                    continue

                i += 1

            return addends

        def collect_scales(addends: Sequence[VectorExpr]) -> dict[VectorExpr, Expr]:
            mapping: defaultdict[VectorExpr, Expr] = defaultdict(lambda: S.Zero)

            for addend in addends:
                if isinstance(addend, VectorScale):
                    mapping[addend.args[0]] += addend.args[1]
                else:
                    mapping[addend] += S.One

            return mapping

        def filter_scales(mapping: dict[VectorExpr, Expr]) -> dict[VectorExpr, Expr]:
            excluded_keys = []

            for vector, scale in mapping.items():
                if scale == 0:
                    excluded_keys.append(vector)

            for key in excluded_keys:
                del mapping[key]

            return mapping

        mapping = filter_scales(collect_scales(flatten_additions(self.args)))

        scaled_addends = [vector * scale for vector, scale in mapping.items()]

        match len(scaled_addends):
            case 0:
                return ZERO
            case 1:
                return scaled_addends[0]
            case _:
                return VectorAdd(*scaled_addends)

    # TODO: order the addends in __init__ so that we could return a consistent `_hashable_contents`
    # tuple and get rid of custom __eq__
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, VectorAdd):
            return False

        return set(self.args) == set(other.args)

    def __hash__(self) -> int:
        return hash(self.args)

    def _sympystr(self, p: Printer) -> str:
        return " + ".join(map(p.doprint, self.args))

    def _eval_evalf(self, prec: int) -> VectorAdd:
        return VectorAdd(*(addend.evalf(n=prec) for addend in self.args)).doit()


__all__ = [
    "ZERO",
    "Scale",
    "VectorAdd",
    "VectorExpr",
    "VectorNorm",
    "VectorScale",
    "VectorSymbol",
    "norm",
]
