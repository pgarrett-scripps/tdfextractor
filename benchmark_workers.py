"""
Benchmark script to compare MGF extraction performance with different worker counts.
"""

import time
import tempfile
import shutil
from pathlib import Path
import logging
import argparse

from src.tdfextractor.mgf_exctractor import process_single_d_folder
from concurrent.futures import ThreadPoolExecutor, as_completed

# Set up logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

class BenchmarkArgs:
    """Mock args object for benchmarking."""
    def __init__(self):
        self.remove_precursor = False
        self.precursor_peak_width = 2.0
        self.batch_size = 100
        self.top_n_peaks = None
        self.min_spectra_intensity = None
        self.max_spectra_intensity = None
        self.min_spectra_mz = None
        self.max_spectra_mz = None
        self.min_precursor_intensity = None
        self.max_precursor_intensity = None
        self.min_precursor_charge = None
        self.max_precursor_charge = None
        self.min_precursor_mz = None
        self.max_precursor_mz = None
        self.min_precursor_rt = None
        self.max_precursor_rt = None
        self.min_precursor_ccs = None
        self.max_precursor_ccs = None
        self.min_precursor_neutral_mass = None
        self.max_precursor_neutral_mass = None
        self.keep_empty_spectra = False
        self.overwrite = True

def benchmark_sequential(d_folders, output_dir, args):
    """Benchmark sequential processing."""
    start_time = time.time()
    
    successful_count = 0
    for d_folder in d_folders:
        try:
            success = process_single_d_folder(d_folder, args, output_dir, None)
            if success:
                successful_count += 1
        except Exception as e:
            logger.error(f"Error processing {d_folder}: {e}")
    
    end_time = time.time()
    return end_time - start_time, successful_count

def benchmark_parallel(d_folders, output_dir, args, num_workers):
    """Benchmark parallel processing with specified number of workers."""
    start_time = time.time()
    
    successful_count = 0
    failed_count = 0
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all jobs
        future_to_folder = {
            executor.submit(process_single_d_folder, d_folder, args, output_dir, None): d_folder
            for d_folder in d_folders
        }
        
        # Process completed jobs
        for future in as_completed(future_to_folder):
            d_folder = future_to_folder[future]
            try:
                success = future.result()
                if success:
                    successful_count += 1
                else:
                    failed_count += 1
            except Exception as e:
                logger.error(f"Unexpected error processing {d_folder}: {e}")
                failed_count += 1
    
    end_time = time.time()
    return end_time - start_time, successful_count

def run_benchmark():
    """Run the benchmark comparison."""
    # Find test data directory
    test_data_dir = Path("tests/data")
    if not test_data_dir.exists():
        logger.error(f"Test data directory not found: {test_data_dir}")
        return
    
    # Find .d folders
    d_folders = list(test_data_dir.glob("*.d"))
    if not d_folders:
        logger.error(f"No .d folders found in: {test_data_dir}")
        return
    
    logger.info(f"Found {len(d_folders)} .d folders: {[f.name for f in d_folders]}")
    
    # Validate .d folders
    valid_d_folders = []
    for d_folder in d_folders:
        if (d_folder / "analysis.tdf").exists() and (d_folder / "analysis.tdf_bin").exists():
            valid_d_folders.append(d_folder)
        else:
            logger.warning(f"Skipping invalid .d folder: {d_folder}")
    
    if not valid_d_folders:
        logger.error("No valid .d folders found")
        return
    
    logger.info(f"Using {len(valid_d_folders)} valid .d folders for benchmark")
    
    # Create mock args
    args = BenchmarkArgs()
    
    # Run benchmarks
    num_runs = 3  # Run multiple times for better average
    results = {}
    
    for workers in [1, 2]:
        times = []
        logger.info(f"\n{'='*50}")
        logger.info(f"Benchmarking with {workers} worker(s)")
        logger.info(f"{'='*50}")
        
        for run in range(num_runs):
            # Create temporary output directory
            with tempfile.TemporaryDirectory() as temp_dir:
                output_dir = Path(temp_dir)
                
                logger.info(f"Run {run + 1}/{num_runs} with {workers} worker(s)...")
                
                if workers == 1:
                    elapsed_time, successful_count = benchmark_sequential(valid_d_folders, output_dir, args)
                else:
                    elapsed_time, successful_count = benchmark_parallel(valid_d_folders, output_dir, args, workers)
                
                times.append(elapsed_time)
                logger.info(f"Run {run + 1} completed: {elapsed_time:.2f}s, {successful_count} files processed")
        
        avg_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)
        
        results[workers] = {
            'times': times,
            'avg': avg_time,
            'min': min_time,
            'max': max_time
        }
        
        logger.info(f"Average time with {workers} worker(s): {avg_time:.2f}s")
        logger.info(f"Min: {min_time:.2f}s, Max: {max_time:.2f}s")
    
    # Print comparison
    logger.info(f"\n{'='*60}")
    logger.info("BENCHMARK RESULTS SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"Files processed: {len(valid_d_folders)}")
    logger.info(f"Runs per configuration: {num_runs}")
    logger.info("")
    
    for workers in sorted(results.keys()):
        result = results[workers]
        logger.info(f"{workers} Worker(s):")
        logger.info(f"  Average: {result['avg']:.2f}s")
        logger.info(f"  Range: {result['min']:.2f}s - {result['max']:.2f}s")
        logger.info(f"  Individual runs: {[f'{t:.2f}s' for t in result['times']]}")
        logger.info("")
    
    # Calculate speedup
    if 1 in results and 2 in results:
        speedup = results[1]['avg'] / results[2]['avg']
        logger.info(f"Speedup with 2 workers: {speedup:.2f}x")
        
        if speedup > 1.1:
            logger.info("✅ Parallel processing shows significant improvement!")
        elif speedup > 1.0:
            logger.info("✅ Parallel processing shows modest improvement")
        else:
            logger.info("⚠️  Parallel processing is slower (likely due to overhead)")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark MGF extraction with different worker counts")
    args = parser.parse_args()
    
    run_benchmark()
