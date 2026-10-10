# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Fixed

- Release HDF5 scalar string buffers after reads, including failed conversions, and owned observables when clearing an `ObservableSet`.
- Correct scalar/vector Monte Carlo result negation and keep uncertainties nonnegative under negative scaling or reciprocal arithmetic. Preserve empty histogram ranges, avoid invalid access during empty vector/valarray conversions, and create private temporary files on Unix independently of the caller's umask.

### Changed

- Keep Python tests and fixtures alongside pyalps and run standalone numerical tutorial tests independently. Retire the optional Boost.Python differential build audit after retaining application checks and focused regressions; see [the coverage assessment](tests/python-migration.md) for historical compatibility limits.
