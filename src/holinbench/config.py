"""Configuration loading and path resolution.

The whole pipeline is driven by ``config/config.yaml``. This module loads it into
a light ``Config`` wrapper that supports dotted access and resolves relative paths
against the repository root (the parent of the ``config/`` directory).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

log = logging.getLogger("holinbench")


class Config:
    """Dotted-access wrapper around the parsed YAML config."""

    def __init__(self, data: dict, config_path: Path):
        self._data = data
        self.config_path = Path(config_path).resolve()
        # Repo root = parent of the directory containing the config file
        # (config/config.yaml -> repo root). Falls back to config dir's parent.
        self.root = self.config_path.parent.parent

    # -- mapping-ish access ---------------------------------------------------
    def __getitem__(self, key: str) -> Any:
        return self._data[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self._data.get(key, default)

    def section(self, name: str) -> dict:
        """Return a top-level section as a dict (empty dict if absent)."""
        return self._data.get(name, {}) or {}

    def dotted(self, path: str, default: Any = None) -> Any:
        """Fetch a nested value with a dotted path, e.g. 'scoring.weights.hmm'."""
        node: Any = self._data
        for part in path.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    # -- paths ----------------------------------------------------------------
    def resolve(self, path_like: str | Path) -> Path:
        """Resolve a config path against the repo root unless absolute."""
        p = Path(path_like)
        return p if p.is_absolute() else (self.root / p)

    def path(self, key: str) -> Path:
        """Resolve a path from the ``paths`` section by key."""
        return self.resolve(self.section("paths")[key])

    @property
    def output_dir(self) -> Path:
        return self.resolve(self.section("paths").get("output_dir", "results"))

    # -- standard output sub-directories --------------------------------------
    def out(self, *parts: str, mkdir: bool = True) -> Path:
        p = self.output_dir.joinpath(*parts)
        if mkdir:
            p.parent.mkdir(parents=True, exist_ok=True)
        return p

    def out_dir(self, *parts: str) -> Path:
        p = self.output_dir.joinpath(*parts)
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def seed(self) -> int:
        return int(self.dotted("project.random_seed", 1729))

    @property
    def raw(self) -> dict:
        return self._data


def load_config(config_path: str | Path) -> Config:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return Config(data, path)
