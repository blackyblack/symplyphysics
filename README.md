# Description

Physics laws implemented as code. Based on [SymPy](https://github.com/sympy/sympy) Python library.

# How to install

```sh
python3 setup.py install --user
```

Install **matplotlib** for plotting support:

```sh
pip install matplotlib
```

# How to run

With **symplyphysics** installed:

```sh
cd examples/force_from_acceleration
python3 main.py
```

Without **symplyphysics** installed:

```sh
PYTHONPATH=. python3 ./examples/force_from_acceleration/main.py
```

