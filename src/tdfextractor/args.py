"""
Argument dataclasses for tdfextractor writers.

The :class:`BaseExtractorArgs` dataclass holds every option shared between the
MS2, MGF, and mzML writers. Format-specific subclasses add their own fields:

* :class:`Ms2Args` and :class:`MgfArgs` add text-output precision controls.
* :class:`MzmlArgs` adds mzML-specific compression and encoding controls.

Each subclass exposes a :meth:`from_namespace` classmethod that builds an
instance from an :class:`argparse.Namespace` produced by one of the parsers in
:mod:`tdfextractor.cli_args`. CLI flags that are not fields on the dataclass
(e.g. ``workers``, ``verbose``, ``overwrite``, preset flags) are silently
dropped — those are handled by the ``main()`` driver, not the writer.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields
from typing import Literal

EncodingBitWidth = Literal[32, 64]
"""Allowed bit widths for mzML binary array encoding."""

CompressionName = Literal[
    "none",
    "zlib",
    "zstd",
    "numpress-linear",
    "numpress-slof",
    "numpress-pic",
]
"""Allowed compression identifiers for mzML binary arrays."""


@dataclass
class BaseExtractorArgs:
    """Fields common to MS2, MGF, and mzML extractors."""

    # I/O
    analysis_dir: str
    output_file: str | None = None

    # Precursor removal
    remove_precursor: bool = False
    precursor_peak_width: float = 2.0

    # Processing
    batch_size: int = 100
    keep_empty_spectra: bool = False

    # Spectra/peak filters
    top_n_peaks: int | None = None
    min_spectra_intensity: float | None = None
    max_spectra_intensity: float | None = None
    min_spectra_mz: float | None = None
    max_spectra_mz: float | None = None

    # Precursor filters
    min_precursor_intensity: float | None = None
    max_precursor_intensity: float | None = None
    min_precursor_charge: int | None = None
    max_precursor_charge: int | None = None
    min_precursor_mz: float | None = None
    max_precursor_mz: float | None = None
    min_precursor_rt: float | None = None
    max_precursor_rt: float | None = None
    min_precursor_ccs: float | None = None
    max_precursor_ccs: float | None = None
    min_precursor_neutral_mass: float | None = None
    max_precursor_neutral_mass: float | None = None

    @classmethod
    def from_namespace(cls, ns: argparse.Namespace):
        """Build an args instance from an argparse.Namespace.

        Only attributes whose names match a field on ``cls`` are pulled across.
        ``ns.output`` is translated to ``output_file``. Unrelated CLI flags
        (``workers``, ``verbose``, ``overwrite``, preset flags) are silently
        dropped. Argparse defaults must line up with the dataclass defaults so
        that ``None`` from a missing CLI flag is the intended ``None``-as-no-
        filter sentinel and not an accidental override.
        """
        ns_dict = vars(ns)
        kwargs: dict = {}
        for f in fields(cls):
            if f.name == "output_file":
                if "output" in ns_dict:
                    kwargs["output_file"] = ns_dict["output"]
            elif f.name in ns_dict:
                kwargs[f.name] = ns_dict[f.name]
        return cls(**kwargs)


@dataclass
class Ms2Args(BaseExtractorArgs):
    """Arguments for :func:`tdfextractor.ms2_extractor.write_ms2_file`."""

    mz_precision: int = 5
    intensity_precision: int = 0


@dataclass
class MgfArgs(BaseExtractorArgs):
    """Arguments for :func:`tdfextractor.mgf_exctractor.write_mgf_file`."""

    mz_precision: int = 5
    intensity_precision: int = 0


@dataclass
class MzmlArgs(BaseExtractorArgs):
    """Arguments for :func:`tdfextractor.mzml_extractor.write_mzml_file`."""

    include_ms1: bool = True

    mz_compression: CompressionName = "zlib"
    intensity_compression: CompressionName = "zlib"
    mobility_compression: CompressionName = "zlib"

    mz_encoding: EncodingBitWidth = 64
    intensity_encoding: EncodingBitWidth = 32

    def __post_init__(self) -> None:
        for name in ("mz_encoding", "intensity_encoding"):
            value = getattr(self, name)
            if value not in (32, 64):
                raise ValueError(
                    f"{name} must be 32 or 64, got {value!r}"
                )

    @classmethod
    def from_namespace(cls, ns: argparse.Namespace) -> "MzmlArgs":
        obj = super().from_namespace(ns)
        # The CLI exposes --no-ms1 (store_true) but the dataclass uses
        # the positive include_ms1 form, so invert here.
        if hasattr(ns, "no_ms1"):
            obj.include_ms1 = not ns.no_ms1
        return obj  # type: ignore[return-value]
