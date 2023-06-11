"""
String template for ms2 file header
"""

MS2_HEADER = 'H\tExtractor\tTimsTOF_extractor\n' \
                 'H\tExtractorVersion\t{version}\n' \
                 'H\tPublicationDate\t20-02-2020\n' \
                 'H\tDate Converted\t{date_of_creation}\n' \
                 'H\tComments\tTimsTOF_extractor written by Yu Gao, 2018\n' \
                 'H\tComments\tTimsTOF_extractor modified by Titus Jung, 2019\n' \
                 'H\tComments\tTimsTOF_extractor modified by Patrick Garrett, April 13th 2023\n' \
                 'H\tMinimumMsMsIntensity\t{minimum_intensity}\n' \
                 'H\tRemoveCharge1\t{remove_charge1}\n' \
                 'H\tRemoveEmptySpectra\t{remove_empty_spectra}\n' \
                 'H\tIncludeSpectra\t{include_spectra}\n' \
                 'H\tExtractorOptions\tMSn\n' \
                 'H\tAcquisitionMethod\t{method}\n' \
                 'H\tInstrumentType\tTIMSTOF\n' \
                 'H\tDataType\tCentroid\n' \
                 'H\tScanType\tMS2\n' \
                 'H\tResolution{resolution}\n' \
                 'H\tIsolationWindow\n' \
                 'H\tFirstScan\t{first_scan:d}\n' \
                 'H\tLastScan\t{last_scan:d}\n' \
                 'H\tMonoIsotopic PrecMz\tTrue\n'
