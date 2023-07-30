# Description

Physics laws implemented as code. Based on [SymPy](https://github.com/sympy/sympy) Python library.

# Sample generated plots

![Carnot Cycle](./img/carnot_cycle.png)
![Floating Body](./img/floating_body.png)

# How to install

```sh
pip install .
```

Install with **matplotlib** for plotting support:

```sh
pip install .[plots]
```

# How to install for development (local installation)

```sh
pip install -e .[dev,plots]
```

# How to run

Install **symplyphysics** and run:

```sh
cd examples/force_from_acceleration
python3 main.py
```

# How to test

Install with **pytest**:

```sh
pip install .[dev]
```

Run tests:

```sh
pytest
```

> **_NOTE:_**  for Windows users **Python/Scripts** folder should be added to the PATH environment variable
