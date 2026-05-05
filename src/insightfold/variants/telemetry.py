from __future__ import annotations

import json
import urllib.parse
from typing import Dict, Optional

import requests


def build_prefilled_google_form_link(
    payload: Dict[str, str],
    *,
    base_url: str,
    field_map: Optional[Dict[str, str]] = None,
) -> str:
    """
    Build a prefilled Google Form URL.

    - base_url: the public Form URL ending with '/viewform'.
    - field_map: optional mapping from payload keys to Google 'entry.<id>' field names.
      If omitted, we simply append the payload as a query string for convenience.
    """
    if field_map:
        q = {}
        for k, v in payload.items():
            entry_key = field_map.get(k)
            if entry_key:
                q[entry_key] = v
        query = urllib.parse.urlencode(q, doseq=True, safe=",:;{}[]\"'")
    else:
        query = urllib.parse.urlencode(payload, doseq=True, safe=",:;{}[]\"'")
    return f"{base_url}?{query}"


def build_github_issue_link(
    *, new_issue_base_url: str, title: str, body: str
) -> str:
    """
    Build a link that opens a prefilled GitHub Issue or Discussion.

    - new_issue_base_url typically looks like:
      https://github.com/<org>/<repo>/issues/new
    """
    params = {"title": title, "body": body}
    return f"{new_issue_base_url}?{urllib.parse.urlencode(params)}"


def send_webhook_event(url: str, payload: Dict, *, timeout: int = 10) -> bool:
    """
    Optional: POST a telemetry event to a webhook you control.
    Returns True on 2xx, False otherwise. Never raises in normal use.
    """
    try:
        r = requests.post(url, json=payload, timeout=timeout)
        return 200 <= r.status_code < 300
    except Exception:
        return False
