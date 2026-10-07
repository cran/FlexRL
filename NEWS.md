# FlexRL 1.0.0 (release candidate)

- Functions and arguments renamed, data configuration functions added.
- New `survival_model()`: the change over time of a structured PIV can be
  modelled with an exponential (default), Weibull, Gompertz, piecewise-constant
  or user-defined survival function, passed to `StEM()` through the 
  `survival_model` argument. `alpha` chains are named after the model parameters. 
  `SurvivalUnstable()` and `loglikSurvival()` are obsolete. `cond_hazard_params` 
  are now given on the parameter scale of the model (log-hazards for the
  exponential model).
- New `RL_diagnostics()`: provides record linkage diagnostics over linkage 
  decision rules (thresholds from 0.50 to 0.99); internally computes 
  FDP_measures and discrepancy_measures useful for downstream inference;
  `print(x, threshold = )` and `plot(x, type = )` display the diagnostics,
- New `link_with_<pkg>()` record linkage wrappers created for several packages:
  multilink, fastLink, BRL, reclin2, diyar, fedmatch, FlexRL. They return a 
  named list (`idxA`, `idxB`, `LinkScore`) that can be passed directly 
  to `RL_diagnostics()`.
- Useless dependencies dropped; R >= 4.1.0.
- Bugs corrected, warnings added.

# FlexRL 0.1.0 (development)

- Package creation. Initial CRAN submission.
