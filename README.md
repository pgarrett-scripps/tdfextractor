# tdfextractor

A Python package to extract MS/MS spectra from Bruker TimsTOF .D folders and convert them to standard formats (MS2 and MGF).

## Installation

```bash
pip install tdfextractor
```

## Usage

tdfextractor provides two command-line tools for extracting spectra:

### MS2 Extraction
Extract MS2 format files (compatible with MS-GF+, Comet, etc.):

```bash
ms2-extractor /path/to/sample.d
ms2-extractor /path/to/sample.d --output custom_output.ms2 --min-intensity 100 --min-charge 2
```

### MGF Extraction  
Extract MGF format files (compatible with Mascot, MaxQuant, etc.):

```bash
mgf-extractor /path/to/sample.d
mgf-extractor /path/to/sample.d --casanovo  # Optimized for Casanovo de novo sequencing
```

## Command Line Arguments

### MS2 Extractor Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `analysis_dir` | str | - | Path to the .D analysis directory |
| `-o, --output` | str | `<analysis_dir_name>.ms2` | Output MS2 file path |
| `--remove-precursor` | flag | False | Remove precursor peaks from MS/MS spectra |
| `--precursor-peak-width` | float | 2.0 | Width around precursor m/z to remove (Da) |
| `--batch-size` | int | 100 | Batch size for processing spectra |
| `--top-n-spectra` | int | None | Keep only top N most intense peaks per spectrum |
| `--min-intensity` | float | 0.0 | Minimum intensity threshold for peaks |
| `--min-charge` | int | None | Minimum charge state filter |
| `--max-charge` | int | None | Maximum charge state filter |
| `--min-mz` | float | None | Minimum m/z filter |
| `--max-mz` | float | None | Maximum m/z filter |
| `--min-rt` | float | None | Minimum retention time filter (seconds) |
| `--max-rt` | float | None | Maximum retention time filter (seconds) |
| `--min-ccs` | float | None | Minimum CCS filter |
| `--max-ccs` | float | None | Maximum CCS filter |
| `-v, --verbose` | flag | False | Enable verbose logging |

### MGF Extractor Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `analysis_dir` | str | - | Path to the .D analysis directory |
| `-o, --output` | str | `<analysis_dir_name>.mgf` | Output MGF file path |
| `--remove-precursor` | flag | False | Remove precursor peaks from MS/MS spectra |
| `--precursor-peak-width` | float | 2.0 | Width around precursor m/z to remove (Da) |
| `--batch-size` | int | 100 | Batch size for processing spectra |
| `--top-n-spectra` | int | None | Keep only top N most intense peaks per spectrum |
| `--min-intensity` | float | 0.0 | Minimum intensity threshold for peaks |
| `--min-charge` | int | None | Minimum charge state filter |
| `--max-charge` | int | None | Maximum charge state filter |
| `--min-mz` | float | None | Minimum m/z filter |
| `--max-mz` | float | None | Maximum m/z filter |
| `--min-rt` | float | None | Minimum retention time filter (seconds) |
| `--max-rt` | float | None | Maximum retention time filter (seconds) |
| `--min-ccs` | float | None | Minimum CCS filter |
| `--max-ccs` | float | None | Maximum CCS filter |
| `-v, --verbose` | flag | False | Enable verbose logging |
| `--casanovo` | flag | False | Preset for Casanovo: enables precursor removal, top-150 peaks, min intensity 0.01, m/z 50-2500 |
