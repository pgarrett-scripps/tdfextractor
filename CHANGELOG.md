# Changelog

All notable changes to this project will be documented in this file.

## [0.4.0]
### Changed (BREAKING)
- `write_ms2_file`, `write_mgf_file`, and `write_mzml_file` now take a single
  dataclass argument (`Ms2Args`, `MgfArgs`, `MzmlArgs`) instead of 25-31
  individual keyword arguments. The dataclasses live in
  `tdfextractor.args` and are re-exported from the package root.
- `generate_header` in `ms2_extractor` now takes an `Ms2Args` instead of 23
  individual kwargs.
- `mz_encoding` and `intensity_encoding` are now typed
  `EncodingBitWidth = Literal[32, 64]`. Invalid values raise `ValueError`
  at `MzmlArgs` construction time instead of inside the writer.
- `mz_compression` / `intensity_compression` / `mobility_compression` are
  typed `CompressionName` (a `Literal` over the supported codec names).
- CLI defaults for `--mz-precision` / `--intensity-precision` are now `5` /
  `0` to match the dataclass defaults (previously `None`, downstream code
  treated `None` as 5/0 anyway, so behavior is unchanged for CLI users).
- mzML extractor: `--min-precursor-rt` / `--max-precursor-rt` now also bound
  the MS1 frames written to the mzML file. Previously these flags only
  filtered MS2 spectra/windows, leaving every MS1 frame in the file
  regardless of the RT window. The new behavior produces a coherent RT-
  bounded slice of the run for both DDA and DIA/PRM acquisitions.

### Added
- New `tdfextractor.args` module exposing `BaseExtractorArgs`, `Ms2Args`,
  `MgfArgs`, `MzmlArgs`, `EncodingBitWidth`, and `CompressionName`. Each
  args class has a `from_namespace(argparse.Namespace)` classmethod.
- `write_ms2_file` and `write_mgf_file` are now part of the public API
  (`from tdfextractor import ...`).
- New `tests/test_mgf_extractor.py` with full MGF format coverage.
- New `tests/test_args.py` with fast unit tests for the new dataclasses.
- New `tests/conftest.py` with session-scoped extraction fixtures so each
  acquisition type only runs through the writer once per test session.
- Added a `slow` pytest marker; run `pytest -m "not slow"` to skip the
  end-to-end extraction tests during iteration.
- `--workers` arg for parallel processing of multiple `.d` folders.

### Migration
```python
# before
write_mzml_file("/path/to/foo.d", output_file="foo.mzML", top_n_peaks=150)

# after
from tdfextractor import MzmlArgs, write_mzml_file
write_mzml_file(MzmlArgs(
    analysis_dir="/path/to/foo.d",
    output_file="foo.mzML",
    top_n_peaks=150,
))
```

## [0.3.0]
### Added
- more args
- mgf-ex & ms2-ex command shorthand
- fixed mgf pepmass

## [0.2.0]
### Added
- mgf file
- cli
- more args
- readme
- updated ms2 args
- rm linting action
- formatted with black

## [0.1.3]
### Added
- switched to serenipy as the backend for ms2 file creation
- removed string templates for ddams2spectra and ddapeakline
- get_contents now returns ms2_spectra, rather than strings
- added tests
- using context manager for timsdata

## [0.1.4]
### Added
- batch process msms spectra
- updated serenipy to 0.2.6

## [0.1.5]
### Added
- added iso width and mz
- converted intensity values to ints
- added min_intensity option
- added tqdm support
- improved logging/readability
- changed constants.MS2_VERSION to be extractor version

## [0.1.6]
### Added
- updated to tdfpy==0.1.6 

## [0.1.7]
### Added
- fixed requirements

## [0.1.7]
### Changed
- src based
- tdfpy==0.1.7
- Ms2 header
- Merged dataframes instead of keeping dicts
### Added
- workflows: pylint, pytest, pypi
- PRM workflow to ms2 extractor
- More I lines 
