"""Where the manuscript sources live, for the checks that compare text to data.

The repository carries no manuscript.  Where one is present, its directory is
named by the KIE_MANUSCRIPT environment variable or, failing that, by the first
line of analysis/manuscript_root.local, a file kept out of version control.
With neither, the checks that need the text are skipped and everything computed
from the data still runs.

Inside that directory the main text is the one file matching *_manuscript.tex
and the supplement body is si/si_body.tex, so no document name is written into
the code.
"""
from __future__ import annotations

import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOCAL = HERE / "manuscript_root.local"


def requested():
    """The manuscript directory the caller asked for, or None."""
    env = os.environ.get("KIE_MANUSCRIPT")
    if env:
        return Path(env).resolve()
    if LOCAL.exists():
        lines = [ln.strip() for ln in LOCAL.read_text().splitlines() if ln.strip()]
        if lines:
            p = Path(lines[0])
            return (p if p.is_absolute() else HERE / p).resolve()
    return None


def root() -> Path:
    """The requested directory, or a placeholder that does not exist."""
    return requested() or (HERE.parent / "manuscript")


def one(directory: Path, pattern: str) -> Path:
    """The single file matching pattern; a placeholder path if there is none."""
    hits = sorted(p for p in directory.glob(pattern) if p.is_file())
    if len(hits) > 1:
        raise RuntimeError(f"{len(hits)} files match {pattern} in {directory}; "
                           "expected one")
    return hits[0] if hits else directory / pattern.replace("*", "missing")
