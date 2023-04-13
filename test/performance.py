import time
from pathlib import Path

from tdfextractor.ms2_extractor import get_ms2_content, write_ms2_file

d_folder = str(Path('data') / '200ngHeLaPASEF_1min.d')

start_time = time.time()
write_ms2_file(d_folder, include_spectra=True, batch_size=1000)
print(time.time()-start_time)

"""for i in range(1, 1000, 100):
    start_time = time.process_time()

    _ = list(get_ms2_content(d_folder, include_spectra=True, batch_size=i))

    print(time.process_time()-start_time, i)"""
