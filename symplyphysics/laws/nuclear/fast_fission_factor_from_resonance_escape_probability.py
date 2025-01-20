# Description
## The fast fission factor is the ratio of the fast neutrons produced by fissions at all energies to
## the number of fast neutrons produced in thermal fission.

## Law: ε ≈ 1 + ((1 - p) / p) * uf * vf * Pfaf / (f * vt * Ptaf * Ptnl)
## Where:
## p - resonance escape probability.
##   See [resonance escape probability](./resonance_escape_probability_from_resonance_absorption_integral.py) implementation.
## uf (fast utilization) - probability that a fast neutron is absorbed in fuel.
## vf - average number of fast neutrons produced per fission.
## vt - average number of thermal neutrons produced per fission.
## f - thermal utilisation factor.
##   See [thermal utilisation factor](./thermal_utilisation_factor_from_macroscopic_absorption_cross_sections.py) implementation.
## Pfaf -  probability that a fast neutron absorption in fuel causes fission.
## Ptaf -  probability that a thermal neutron absorption in fuel causes fission.
## Ptnl - thermal non-leakage factor.
##   See [thermal non_leakage probability](./thermal_non_leakage_probability_from_diffusion_length.py) implementation.
## ε - fast fission factor

## Unfortunately proper fast fission factor formula is too complex for actual use. Approximate formula is not
## very useful since there is no actual data on Pfaf, Ptaf and uf values. Therefore this module is empty,
## only present as a documentation for the fast fission factor.

# Links:
## Wikipedia, fourth row in table <https://en.wikipedia.org/wiki/Six_factor_formula>
