"""
InsightFold – Variants v1
Thin public API for the notebook and future CLI.
"""

from .io import (
    UNIPROT_VARIATION_API,
    AFDB_PREDICTION_URL,
    ModelInfo,
    build_afdb_id,
    fetch_uniprot_entry,
    afdb_prediction_meta,
    fetch_model_info,
    crc64_iso,
)
from .normalise import (
    UniProtFeature,
    parse_uniprot_features,
)
from .map import tidy_from_features
from .viz import show_structure_with_pins, centre_on_position
from .telemetry import (
    build_prefilled_google_form_link,
    build_github_issue_link,
    send_webhook_event,
)

__all__ = [
    # IO
    "UNIPROT_VARIATION_API",
    "AFDB_PREDICTION_URL",
    "ModelInfo",
    "build_afdb_id",
    "fetch_uniprot_entry",
    "afdb_prediction_meta",
    "fetch_model_info",
    "crc64_iso",
    # Normalisation
    "UniProtFeature",
    "parse_uniprot_features",
    # Mapping
    "tidy_from_features",
    # Visuals
    "show_structure_with_pins",
    "centre_on_position",
    # Telemetry
    "build_prefilled_google_form_link",
    "build_github_issue_link",
    "send_webhook_event",
]

__version__ = "0.1.0"
