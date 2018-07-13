# Changelog
Torsten project changelog

## [Unrelesed]

## [0.84] - 2018-02-24
### Added
- Piecewise linear interpolation function.
- Univariate integral functions.

### Changed
- Update with Stan version 2.17.1.
- Minor revisions to User Manual.
- Bugfixes.

## [0.83] - 2017-08-02
### Added
- Work with TorstenHeaders
- Each chain has a different initial estimate

### Changed
- User manual
- Fix misspecification in ODE system for TwoCpt example.
- Other bugfixes


## [0.82] - 2017-01-29
### Added
- Allow parameter arguments to be passed as 1D or 2D arrays
- More unit tests
- Unit tests check automatic differentiation against finite differentiation.

### Changed
- Split the parameter argument into three arguments: pMatrix
  (parameters for the ODEs -- note: for linOdeModel, pMatrix
  is replaced by the constant rate matrix K), biovar
  (parameters for the biovariability), and tlag (parameters
  for the lag time).
- bugfixes

## [0.81] - 2016-09-27
### Added
- linCptModel (linear compartmental model) function

## [0.80a] - 2016-09-21
### Added
- check_finite statements in pred_1 and pred_2 to reject metropolis proposal if initial conditions are not finite
