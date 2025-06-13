# Changelog

All notable changes to this project will be documented in this file.

## [0.4.0]
### Added
- workers arg
- benchmark script

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
