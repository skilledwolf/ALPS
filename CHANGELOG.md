# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Fixed

- Release HDF5 scalar string buffers after reads, including failed conversions, and owned observables when clearing an `ObservableSet`.
- Correct scalar/vector Monte Carlo result negation and keep uncertainties nonnegative under negative scaling or reciprocal arithmetic. Preserve empty histogram ranges, avoid invalid access during empty vector/valarray conversions, and create private temporary files on Unix independently of the caller's umask.

### Changed

- Remove unused deprecated numerical containers, patched Boost accumulator headers and 295 unused Numeric Bindings headers. Use active ALPS matrix/vector APIs and upstream Boost equivalents; direct users of removed Numeric Bindings headers must provide their own installation. Retain the complete include closure needed by ALPS, IETL and the solvers.
