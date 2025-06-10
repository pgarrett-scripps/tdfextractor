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
ms2-extractor /path/to/directory_with_multiple_d_folders --output /path/to/output_directory
```

### MGF Extraction  
Extract MGF format files (compatible with Mascot, MaxQuant, etc.):

```bash
mgf-extractor /path/to/sample.d
mgf-extractor /path/to/sample.d --casanovo  # Optimized for Casanovo de novo sequencing
mgf-extractor /path/to/directory_with_multiple_d_folders --output /path/to/output_directory
```

## Output Options

Both extractors support flexible output options:

1. **No output specified**: Files are created within each .D folder with auto-generated names
2. **Specific file path**: Use `-o filename.ms2` or `-o filename.mgf` for single .D folder processing
3. **Output directory**: Use `-o /path/to/output_dir` for batch processing multiple .D folders
4. **Overwrite protection**: Use `--overwrite` to replace existing output files

### Batch Processing

When processing multiple .D folders, the extractors will:
- Automatically find all .D folders in the specified directory
- Create output files with names matching the .D folder names
- Skip existing files unless `--overwrite` is specified
- Create the output directory if it doesn't exist

## Command Line Arguments

### MS2 Extractor Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `analysis_dir` | str | - | Path to the .D analysis directory or directory containing .D folders |
| `-o, --output` | str | `<analysis_dir_name>.ms2` | Output MS2 file path or directory |
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
| `--overwrite` | flag | False | Overwrite existing output files |
| `-v, --verbose` | flag | False | Enable verbose logging |

### MGF Extractor Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `analysis_dir` | str | - | Path to the .D analysis directory or directory containing .D folders |
| `-o, --output` | str | `<analysis_dir_name>.mgf` | Output MGF file path or directory |
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
| `--overwrite` | flag | False | Overwrite existing output files |
| `-v, --verbose` | flag | False | Enable verbose logging |
| `--casanovo` | flag | False | Preset for Casanovo: enables precursor removal, top-150 peaks, min intensity 0.01, m/z 50-2500 |

## Features

- **Multiple format support**: Export to MS2 and MGF formats
- **Flexible output options**: Single files, batch processing, custom directories
- **Flexible filtering**: Filter by charge state, m/z range, retention time, CCS, and intensity
- **Batch processing**: Process multiple .D folders at once
- **Precursor removal**: Option to remove precursor peaks from spectra
- **Peak selection**: Keep only the most intense peaks per spectrum
- **DDA and PRM support**: Works with both Data-Dependent Acquisition and Parallel Reaction Monitoring data
- **Overwrite protection**: Prevents accidental file overwrites unless explicitly requested

## Requirements

- Python ≥ 3.8
- tdfpy ≥ 0.1.7
- serenipy ≥ 0.2.6
- tqdm
- pandas
