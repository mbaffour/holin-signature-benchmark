"""Stage 0a — literature search (PubMed E-utilities + Europe PMC).

Retrieves candidate papers for a set of query groups and writes a normalized
``literature_search_results.tsv``. HTTP responses are cached on disk so a run is
reproducible and re-runnable offline. Network failures degrade gracefully: the
command logs a warning and continues with whatever it has.

IMPORTANT: finding a paper that says "holin" is NOT evidence. This stage only
collects papers; evidence.py does the conservative evidence extraction, and
nothing reaches the gold set without manual verification (see PROJECT_SPEC.md).
"""
from __future__ import annotations

import hashlib
import json
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from .config import Config
from .utils import log

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"


# --------------------------------------------------------- HTTP + cache --------
# Credential params that must NOT affect the cache key: the same scientific query
# should reuse cache regardless of who ran it / which API key was supplied.
_CACHE_IGNORED_PARAMS = ("email", "api_key", "tool")


def _cache_path(cache_dir: Path, url: str, params: dict) -> Path:
    # Exclude credentials so cache hits are determined by the query alone.
    key_params = {k: v for k, v in params.items() if k not in _CACHE_IGNORED_PARAMS}
    key = url + "?" + json.dumps(key_params, sort_keys=True)
    h = hashlib.sha1(key.encode()).hexdigest()[:16]
    return cache_dir / f"{h}.cache"


def _http_get(cfg: Config, url: str, params: dict, *, is_json: bool, delay: float):
    cache_dir = cfg.resolve(cfg.dotted("literature.cache_dir", "data/litcache"))
    cache_dir.mkdir(parents=True, exist_ok=True)
    cpath = _cache_path(cache_dir, url, params)
    if cpath.exists():
        text = cpath.read_text(encoding="utf-8")
        return json.loads(text) if is_json else text
    try:
        import requests
    except Exception:
        log.warning("`requests` not available; cannot perform live search.")
        return None
    try:
        resp = requests.get(url, params=params, timeout=30,
                            headers={"User-Agent": "holinbench/0.1 (literature mining)"})
        resp.raise_for_status()
        time.sleep(delay)
    except Exception as exc:
        log.warning("HTTP error for %s: %s", url, type(exc).__name__)
        return None
    text = resp.text
    cpath.write_text(text, encoding="utf-8")
    return resp.json() if is_json else text


# --------------------------------------------------------- dedup helpers -------
def _title_key(p: dict) -> str:
    """Fallback dedup key when no DOI/PMID/PMCID is available.

    A bare 60-char title prefix can collide distinct papers (e.g. shared boiler-
    plate prefixes). Normalize the WHOLE title (lowercase, alnum-only) and append
    the publication year so different papers with similar leading words don't
    merge.
    """
    title = str(p.get("title", "") or "")
    norm = "".join(ch for ch in title.lower() if ch.isalnum())
    if not norm:
        return ""
    year = str(p.get("year", "") or "")
    return f"title:{norm}|{year}"


# Fields whose presence makes a record "richer" for OA full-text retrieval.
_OA_MERGE_FIELDS = ("is_oa", "epmc_source", "epmc_id", "pmcid")


def _merge_paper_records(existing: dict, incoming: dict) -> None:
    """Merge an incoming duplicate record into the kept record in place.

    The first record (often PubMed) wins for stable bibliographic fields, but
    OA/full-text-enabling fields and any field the existing record is missing are
    carried over from whichever source has them, so OA full-text fetch is not lost
    when the OA-aware Europe PMC record arrives second.
    """
    # Always prefer a truthy OA flag and carry over EPMC fetch identifiers.
    for field in _OA_MERGE_FIELDS:
        inc = incoming.get(field)
        if field == "is_oa":
            existing["is_oa"] = bool(existing.get("is_oa")) or bool(inc)
        elif inc and not existing.get(field):
            existing[field] = inc
    # Fill any other field the existing record lacks (e.g. abstract, doi, pmid).
    for field, inc in incoming.items():
        if field in ("_queries",):
            continue
        if inc and not existing.get(field):
            existing[field] = inc


# --------------------------------------------------------- PubMed --------------
def _pubmed_search(cfg, query, retmax, email, api_key, delay) -> list[str]:
    params = {"db": "pubmed", "term": query, "retmax": retmax,
              "retmode": "json", "email": email}
    if api_key:
        params["api_key"] = api_key
    data = _http_get(cfg, f"{EUTILS}/esearch.fcgi", params, is_json=True, delay=delay)
    if not data:
        return []
    return data.get("esearchresult", {}).get("idlist", [])


def _pubmed_fetch(cfg, pmids, email, api_key, delay) -> list[dict]:
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
              "rettype": "abstract", "email": email}
    if api_key:
        params["api_key"] = api_key
    xml = _http_get(cfg, f"{EUTILS}/efetch.fcgi", params, is_json=False, delay=delay)
    if not xml:
        return []
    out = []
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return []
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID", default="")
        title = art.findtext(".//ArticleTitle", default="")
        abstract = " ".join(t.text or "" for t in art.findall(".//AbstractText"))
        year = art.findtext(".//JournalIssue/PubDate/Year", default="") or \
            art.findtext(".//PubDate/MedlineDate", default="")[:4]
        journal = art.findtext(".//Journal/Title", default="")
        authors = "; ".join(
            f"{a.findtext('LastName', '')} {a.findtext('Initials', '')}".strip()
            for a in art.findall(".//Author")[:6])
        doi = ""
        pmcid = ""
        for aid in art.findall(".//ArticleId"):
            t = aid.get("IdType")
            if t == "doi":
                doi = aid.text or ""
            elif t == "pmc":
                pmcid = aid.text or ""
        pubtypes = [pt.text or "" for pt in art.findall(".//PublicationType")]
        out.append({"source": "pubmed", "pmid": pmid, "pmcid": pmcid, "doi": doi,
                    "title": title, "abstract": abstract, "year": year,
                    "journal": journal, "authors": authors,
                    "is_review": any("review" in p.lower() for p in pubtypes)})
    return out


# --------------------------------------------------------- Europe PMC ----------
def _epmc_search(cfg, query, pagesize, delay) -> list[dict]:
    params = {"query": query, "format": "json", "pageSize": pagesize,
              "resultType": "core"}
    data = _http_get(cfg, f"{EPMC}/search", params, is_json=True, delay=delay)
    if not data:
        return []
    out = []
    for r in data.get("resultList", {}).get("result", []):
        out.append({
            "source": "europe_pmc", "pmid": r.get("pmid", ""),
            "pmcid": r.get("pmcid", ""), "doi": r.get("doi", ""),
            "title": r.get("title", ""), "abstract": r.get("abstractText", ""),
            "year": str(r.get("pubYear", "")), "journal": r.get("journalTitle", ""),
            "authors": r.get("authorString", ""),
            "is_review": (r.get("pubType", "") or "").lower().find("review") >= 0,
            "epmc_source": r.get("source", ""), "epmc_id": r.get("id", ""),
            "is_oa": r.get("isOpenAccess", "N") == "Y",
        })
    return out


def _epmc_fulltext(cfg, source, ext_id, delay) -> str:
    if not source or not ext_id:
        return ""
    text = _http_get(cfg, f"{EPMC}/{source}/{ext_id}/fullTextXML", {},
                     is_json=False, delay=delay)
    if not text:
        return ""
    try:
        root = ET.fromstring(text)
    except ET.ParseError:
        return ""
    # Concatenate all text content from the body.
    parts = []
    for elem in root.iter():
        if elem.text and elem.text.strip():
            parts.append(elem.text.strip())
    return " ".join(parts)


# --------------------------------------------------------- driver --------------
def _load_queries(cfg: Config) -> list[tuple[str, str]]:
    qfile = cfg.resolve(cfg.dotted("literature.search_queries_file",
                                   "data/example/search_queries.txt"))
    queries = []
    if qfile.exists():
        for line in qfile.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                grp, q = line.split("\t", 1)
            else:
                grp, q = "general", line
            queries.append((grp.strip(), q.strip()))
    return queries


def run_search(cfg: Config) -> pd.DataFrame:
    lit = cfg.section("literature")
    out_dir = cfg.resolve(lit.get("output_dir", "results/literature"))
    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["paper_id", "query_group", "query", "source", "pmid", "pmcid", "doi",
            "title", "authors", "year", "journal", "abstract", "is_review",
            "has_fulltext", "url"]

    if not lit.get("literature_search_enabled", True):
        log.warning("literature_search_enabled is false; skipping Stage 0 search.")
        empty = pd.DataFrame(columns=cols)
        empty.to_csv(out_dir / "literature_search_results.tsv", sep="\t", index=False)
        return empty

    email = lit.get("ncbi_email") or ""
    api_key = lit.get("ncbi_api_key")
    retmax = int(lit.get("max_papers_per_query", 20))
    include_reviews = lit.get("include_reviews", True)
    delay = float(lit.get("request_delay_seconds", 0.4))
    if api_key:
        delay = min(delay, 0.12)
    fulltext_dir = out_dir / "fulltext"
    fulltext_dir.mkdir(exist_ok=True)

    queries = _load_queries(cfg)
    # add user-supplied PMIDs/DOIs as pseudo-queries
    for pmid in lit.get("user_supplied_pmids", []) or []:
        queries.append(("user_pmid", f"{pmid}[uid]"))

    records: dict[str, dict] = {}
    for grp, q in queries:
        papers = []
        if lit.get("pubmed_enabled", True):
            pmids = _pubmed_search(cfg, q, retmax, email, api_key, delay)
            papers += _pubmed_fetch(cfg, pmids, email, api_key, delay)
        if lit.get("europe_pmc_enabled", True):
            papers += _epmc_search(cfg, q, retmax, delay)
        for p in papers:
            if p.get("is_review") and not include_reviews:
                continue
            key = (p.get("doi") or p.get("pmid") or p.get("pmcid")
                   or _title_key(p))
            if not key:
                continue
            if key in records:
                # Same paper from another source (e.g. PubMed + Europe PMC). Merge
                # so OA / full-text-enabling fields (is_oa, epmc_source, epmc_id)
                # survive even when the first-seen record (PubMed) lacked them.
                _merge_paper_records(records[key], p)
                records[key].setdefault("_queries", set()).add(f"{grp}:{q}")
                continue
            paper_id = "lit_" + hashlib.sha1(key.encode()).hexdigest()[:10]
            url = (f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/" if p.get("pmid")
                   else (f"https://doi.org/{p['doi']}" if p.get("doi") else ""))
            rec = {**p, "paper_id": paper_id, "query_group": grp, "query": q,
                   "has_fulltext": False, "url": url,
                   "_queries": {f"{grp}:{q}"}}
            records[key] = rec

    # Full-text fetch as a POST-PASS over the deduped+merged records, so OA fields
    # that arrived on a later merge (e.g. PubMed first, Europe PMC second) still
    # trigger the OA full-text download.
    if lit.get("pmc_fulltext_enabled", True):
        for rec in records.values():
            if rec.get("has_fulltext"):
                continue
            if rec.get("is_oa") and rec.get("epmc_source") and rec.get("epmc_id"):
                ft = _epmc_fulltext(cfg, rec.get("epmc_source"), rec.get("epmc_id"), delay)
                if ft:
                    (fulltext_dir / f"{rec['paper_id']}.txt").write_text(ft, encoding="utf-8")
                    rec["has_fulltext"] = True

    rows = []
    for rec in records.values():
        rows.append({c: rec.get(c, "") for c in cols})
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(out_dir / "literature_search_results.tsv", sep="\t", index=False)
    log.info("Stage 0: retrieved %d unique paper(s) across %d query group(s) "
             "(%d with full text).", len(df),
             len({g for g, _ in queries}), int(df["has_fulltext"].sum()) if not df.empty else 0)
    return df
