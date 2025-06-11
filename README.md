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

# shorthand
ms2-ex 
ms2-ex /path/to/sample.d --output custom_output.ms2 --min-intensity 100 --min-charge 2
ms2-ex /path/to/directory_with_multiple_d_folders --output /path/to/output_directory
```

### MGF Extraction  
Extract MGF format files

```bash
mgf-extractor /path/to/sample.d

#shorthand
mgf-ex
mgf-ex /path/to/sample.d --casanovo  # Optimized for Casanovo de novo sequencing
mgf-ex /path/to/directory_with_multiple_d_folders --output /path/to/output_directory
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

Both MS2 and MGF extractors share the same arguments, with only a few format-specific options:

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `analysis_dir` | str | - | Path to the .D analysis directory or directory containing .D folders |
| `-o, --output` | str | `<analysis_dir_name>.<ext>` | Output file path or directory |
| `--remove-precursor` | flag | False | Remove precursor peaks from MS/MS spectra |
| `--precursor-peak-width` | float | 2.0 | Width around precursor m/z to remove (Da) |
| `--batch-size` | int | 100 | Batch size for processing spectra |
| `--top-n-peaks` | int | None | Keep only top N most intense peaks per spectrum |
| `--min-spectra-intensity` | float | None | Minimum intensity threshold for MS/MS peaks (absolute or 0.0-1.0 for percentage) |
| `--max-spectra-intensity` | float | None | Maximum intensity threshold for MS/MS peaks (absolute or 0.0-1.0 for percentage) |
| `--min-spectra-mz` | float | None | Minimum m/z filter for MS/MS peaks |
| `--max-spectra-mz` | float | None | Maximum m/z filter for MS/MS peaks |
| `--min-precursor-intensity` | int | None | Minimum precursor intensity filter |
| `--max-precursor-intensity` | int | None | Maximum precursor intensity filter |
| `--min-precursor-charge` | int | None | Minimum precursor charge state filter |
| `--max-precursor-charge` | int | None | Maximum precursor charge state filter |
| `--min-precursor-mz` | float | None | Minimum precursor m/z filter |
| `--max-precursor-mz` | float | None | Maximum precursor m/z filter |
| `--min-precursor-rt` | float | None | Minimum precursor retention time filter (seconds) |
| `--max-precursor-rt` | float | None | Maximum precursor retention time filter (seconds) |
| `--min-precursor-ccs` | float | None | Minimum precursor CCS filter |
| `--max-precursor-ccs` | float | None | Maximum precursor CCS filter |
| `--min-precursor-neutral-mass` | float | None | Minimum precursor neutral mass filter |
| `--max-precursor-neutral-mass` | float | None | Maximum precursor neutral mass filter |
| `--mz-precision` | int | 5 | Number of decimal places for m/z values |
| `--intensity-precision` | int | 0 | Number of decimal places for intensity values |
| `--keep-empty-spectra` | flag | False | Write empty spectra to output file |
| `--overwrite` | flag | False | Overwrite existing output files |
| `-v, --verbose` | flag | False | Enable verbose logging |

### Format-Specific Arguments

**MS2 Extractor Only:**
- `--ip2`: Use IP2 preset settings (sets min charge to 1)

**MGF Extractor Only:**
- `--casanovo`: Use Casanovo preset settings (enables precursor removal, top-150 peaks, min intensity 0.01, m/z range 50-2500, min charge 1)
