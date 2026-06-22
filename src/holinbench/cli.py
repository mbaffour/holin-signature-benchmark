"""holinbench command-line interface.

Each subcommand maps to a pipeline stage and takes ``--config``. A lazy
``Pipeline`` cache lets individual commands recompute only the prerequisites they
need, while ``run-all`` chains every stage in a single in-memory pass.
"""
from __future__ import annotations

import logging

import click

from . import (alignment, benchmark, clustering, context, features, hmmer,
               motif, plots, report, scoring, synteny, validate)
from . import evidence, literature, seqmap
from .config import load_config
from .utils import setup_logging


# --------------------------------------------------------- lazy pipeline -------
class Pipeline:
    def __init__(self, cfg):
        self.cfg = cfg
        self._c = {}

    def _memo(self, key, fn):
        if key not in self._c:
            self._c[key] = fn()
        return self._c[key]

    @property
    def clean(self):
        return self._memo("clean", lambda: validate.run(self.cfg)["clean"])

    @property
    def features(self):
        return self._memo("features", lambda: features.run(self.cfg, self.clean))

    @property
    def clusters(self):
        return self._memo("clusters", lambda: clustering.run(self.cfg, self.clean))

    @property
    def alignment(self):
        return self._memo("alignment", lambda: alignment.run(
            self.cfg, self.clean, self.features, self.clusters))

    @property
    def built(self):
        return self._memo("built", lambda: hmmer.build(self.cfg, self.alignment))

    @property
    def hits(self):
        return self._memo("hits", lambda: hmmer.scan(
            self.cfg, self.built["models"], self.clean))

    @property
    def context(self):
        return self._memo("context", lambda: context.run(self.cfg, self.features))

    @property
    def ranking(self):
        return self._memo("ranking", lambda: scoring.run(
            self.cfg, self.features, self.context, self.hits,
            restrict_categories=["unknown", "weak"]))

    @property
    def benchmark(self):
        return self._memo("benchmark", lambda: benchmark.run(
            self.cfg, self.features, self.context, self.hits, self.clean))


def _cfg(config):
    return load_config(config)


_config_opt = click.option("--config", "-c", default="config/config.yaml",
                           show_default=True, help="Path to config.yaml")
_verbose_opt = click.option("--verbose", "-v", is_flag=True, help="Debug logging")


@click.group()
def cli():
    """Benchmark universal / topology / family / context-aware holin signatures."""


def _setup(verbose):
    setup_logging(logging.DEBUG if verbose else logging.INFO)


# --------------------------------------------------------- core stages ---------
@cli.command()
@_config_opt
@_verbose_opt
def validate_cmd(config, verbose):
    """Stage 1: validate & clean datasets."""
    _setup(verbose)
    validate.run(_cfg(config))


cli.add_command(validate_cmd, name="validate")


@cli.command(name="features")
@_config_opt
@_verbose_opt
def features_cmd(config, verbose):
    """Stage 6: compute architecture/topology features."""
    _setup(verbose)
    p = Pipeline(_cfg(config))
    p.features


@cli.command(name="cluster")
@_config_opt
@_verbose_opt
def cluster_cmd(config, verbose):
    """Stage 2: cluster gold positives at multiple identities."""
    _setup(verbose)
    Pipeline(_cfg(config)).clusters


@cli.command(name="align")
@_config_opt
@_verbose_opt
def align_cmd(config, verbose):
    """Stage 3: build MSAs + alignment quality metrics."""
    _setup(verbose)
    Pipeline(_cfg(config)).alignment


@cli.command(name="build-hmms")
@_config_opt
@_verbose_opt
def build_hmms_cmd(config, verbose):
    """Stage 4: build universal/topology/family HMMs."""
    _setup(verbose)
    Pipeline(_cfg(config)).built


@cli.command(name="scan")
@_config_opt
@_verbose_opt
def scan_cmd(config, verbose):
    """Stage 5: scan HMMs against all datasets."""
    _setup(verbose)
    Pipeline(_cfg(config)).hits


@cli.command(name="context")
@_config_opt
@_verbose_opt
def context_cmd(config, verbose):
    """Stage 7: genomic-context / lysis-cassette scoring."""
    _setup(verbose)
    Pipeline(_cfg(config)).context


@cli.command(name="score")
@_config_opt
@_verbose_opt
def score_cmd(config, verbose):
    """Stage 8: composite candidate ranking."""
    _setup(verbose)
    Pipeline(_cfg(config)).ranking


@cli.command(name="benchmark")
@_config_opt
@_verbose_opt
def benchmark_cmd(config, verbose):
    """Stage 9: benchmark models A-G (naive + leave-one-family-out)."""
    _setup(verbose)
    Pipeline(_cfg(config)).benchmark


@cli.command(name="motifs")
@_config_opt
@_verbose_opt
def motifs_cmd(config, verbose):
    """Stage 10: conservation/motif analysis + logos."""
    _setup(verbose)
    cfg = _cfg(config)
    Pipeline(cfg).alignment  # ensure alignments exist
    motif.run(cfg)


@cli.command(name="synteny")
@_config_opt
@_verbose_opt
def synteny_cmd(config, verbose):
    """Stage 11: lysis-cassette synteny maps."""
    _setup(verbose)
    p = Pipeline(_cfg(config))
    synteny.run(p.cfg, p.ranking)


@cli.command(name="report")
@_config_opt
@_verbose_opt
def report_cmd(config, verbose):
    """Stage 12: assemble the manuscript-style report."""
    _setup(verbose)
    report.run(_cfg(config))


@cli.command(name="run-all")
@_config_opt
@_verbose_opt
def run_all_cmd(config, verbose):
    """Run stages 1-12 end-to-end on the configured data."""
    _setup(verbose)
    cfg = _cfg(config)
    p = Pipeline(cfg)
    p.clean; p.features; p.clusters; p.alignment; p.built; p.hits; p.context
    ranking = p.ranking
    bench = p.benchmark
    motif.run(cfg)
    plots.run(cfg, feats=p.features, clusters=p.clusters, naive=bench["naive"],
              scores=bench["scores"], fp_fn=bench["fp_fn"], ranking=ranking)
    synteny.run(cfg, ranking)
    out = report.run(cfg)
    click.echo(f"Done. Report: {out}")


# --------------------------------------------------------- Stage 0 -------------
@cli.command(name="literature-search")
@_config_opt
@_verbose_opt
def literature_search_cmd(config, verbose):
    """Stage 0: search PubMed / Europe PMC for candidate holin papers."""
    _setup(verbose)
    literature.run_search(_cfg(config))


@cli.command(name="extract-evidence")
@_config_opt
@_verbose_opt
def extract_evidence_cmd(config, verbose):
    """Stage 0: extract candidate proteins + experimental-evidence sentences."""
    _setup(verbose)
    evidence.run_extract(_cfg(config))


@cli.command(name="map-sequences")
@_config_opt
@_verbose_opt
def map_sequences_cmd(config, verbose):
    """Stage 0: map candidates to accessions / sequences."""
    _setup(verbose)
    seqmap.run_map(_cfg(config))


@cli.command(name="prepare-review")
@_config_opt
@_verbose_opt
def prepare_review_cmd(config, verbose):
    """Stage 0: build the manual review template + curation summary."""
    _setup(verbose)
    evidence.run_prepare_review(_cfg(config))


@cli.command(name="export-curated")
@_config_opt
@_verbose_opt
def export_curated_cmd(config, verbose):
    """Stage 0: export gold_holins.csv from manually verified candidates only."""
    _setup(verbose)
    evidence.run_export_curated(_cfg(config))


def main():
    cli()


if __name__ == "__main__":
    main()
