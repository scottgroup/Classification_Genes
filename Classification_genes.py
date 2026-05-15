from __future__ import annotations

"""
Standalone host-gene classifier for public use.

Minimum input:
    - A CSV file
    - The column name containing host gene Ensembl IDs

Typical usage:
    python3 Classification_genes.py input.csv --host-gene-column host_gene_id

What the script does:
    1. Reads the input CSV.
    2. Uses a MyGene-style annotation file next to this script by default.
    3. Reuses existing annotations from that file and downloads only missing
       gene annotations from MyGene.info.
    4. Downloads the GO ontology file (go-basic.obo) next to this script if
       it is not already present.
    5. Classifies each unique host gene into a broad functional class.
    6. Keeps the annotation file as annotation-only data.
    7. Merges classification results back into the input table.
    8. Writes the merged table back to the input CSV by default, or to a
       separate file when --output is provided.
"""

import argparse
import ast
import json
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd
from goatools.obo_parser import GODag


SCRIPT_DIR = Path(__file__).resolve().parent

GO_BASIC_URL = "https://current.geneontology.org/ontology/go-basic.obo"
MYGENE_GENE_URL = "https://mygene.info/v3/gene"

GO_BASIC_PATH = SCRIPT_DIR / "go-basic.obo"

MYGENE_FIELDS = "symbol,name,taxid,type_of_gene,ensembl,go.BP,go.MF,go.CC"
MYGENE_BATCH_SIZE = 1000

RIBOSOMAL_PROTEIN_MF = {"GO:0003735"}


STRUCTURAL_RIBOSOMAL_SYMBOL_RE = re.compile(
    r"^(?:"
    r"MRP[LS]\d+[A-Z]?"
    r"|RPLP[0-2]"
    r"|RPSA2?"
    r"|RPL36AL"
    r"|RPL\d+(?:[A-Z]|L\d*)?"
    r"|RPS\d+(?:[A-Z]|X|Y\d*|L)?"
    r")$",
    re.IGNORECASE,
)

STRUCTURAL_RIBOSOMAL_NAME_PATTERNS = (
    re.compile(
        r"^(?:40S|60S)\s+ribosomal protein\s+[LS]\d+[A-Z0-9.-]*"
        r"(?:\s+(?:X-linked|Y-linked(?:\s+\d+)?|like(?:\s+\d+)?))?$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^ribosomal protein\s+"
        r"(?:lateral stalk subunit\s+P[0-2]|SA(?:\s+\d+)?|[LS]\d+[A-Z0-9.-]*)"
        r"(?:\s+(?:X-linked|Y-linked(?:\s+\d+)?|like(?:\s+\d+)?))?$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^mitochondrial ribosomal protein\s+[LS]\d+[A-Z0-9.-]*$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?:large|small)\s+ribosomal subunit protein\s+[A-Z]?[LS]\d+[A-Z0-9.-]*$",
        re.IGNORECASE,
    ),
)

# These gene-name fallbacks are used only when no GO annotation is available.
NO_GO_GENE_NAME_KEYWORDS = {
    "Ribosomal protein": (
        STRUCTURAL_RIBOSOMAL_SYMBOL_RE,
    ),
    "Ribosome biogenesis": (
        re.compile(
            r"^(?:BMS1|DKC1|EMG1|FBL|GAR1|IMP[34]|MPHOSPH10|NAF1|NHP2|NOB1|"
            r"NOP(?:10|14|16|2|53|56|58)|NOL(?:6|8|9|10|11|12C|C1)|NSA[12]|"
            r"PES1|PNO1|PWP2|RCL1|RIOK[12]|RRP(?:1|7A|8|9|12|15)|RRS1|SBDS|"
            r"TCOF1|TSR1|UTP(?:3|4|6|11L|14A|14C|15|18|20|23|25)|WDR3|WDR36)$",
            re.IGNORECASE,
        ),
    ),
    "Translation": (
        re.compile(
            r"^(?:AARS1|AIMP[123]|CARS1|DARS1|DENR|EARS[12]|EEF[12][A-Z0-9]*|"
            r"EIF[1-6][A-Z0-9]*|EIF2AK[1-4]|EPRS1|ETF1|FARSA|FARSB|GARS1|"
            r"GSPT[12]|HARS1|IARS[12]|KARS1|LARS[12]|MARS1|MCTS1|NARS1|"
            r"PABP[A-Z0-9]*|PAIP[12]|PARS2?|QARS1|RARS[12]|SARS1|TARS1|"
            r"TMA[0-9A-Z]+|VARS[12]|WARS1|YARS[12])$",
            re.IGNORECASE,
        ),
    ),
    "RNA processing": (
        re.compile(
            r"^(?:CPSF[0-9A-Z]+|CSTF[0-9A-Z]+|DDX(?:23|39A|39B|41|46|47|49|5[0-9])|"
            r"GEMIN[0-9A-Z]+|HNRNP[A-Z0-9]+|LSM[0-9A-Z]+|PRPF[0-9A-Z]+|"
            r"RBM[0-9A-Z]+|SART[0-9]+|SF3[AB][0-9A-Z]+|SMN[12]?|"
            r"SNR(?:NP)?[A-Z0-9]+)$",
            re.IGNORECASE,
        ),
    ),
    "Transcription": (
        re.compile(
            r"^(?:CCNC|CCNH|CDK7|GTF[0-9A-Z]+|MED[0-9A-Z]+|NELF[ABCE]?|"
            r"POLR[123][A-Z0-9]+|SUPT[0-9A-Z]+|TAF[0-9A-Z]+|TBP|TF[23][A-Z0-9]+)$",
            re.IGNORECASE,
        ),
    ),
    "Transport": (
        re.compile(
            r"^(?:CSE1L|IPO[0-9A-Z]+|KPN[AB][0-9A-Z]*|NUP[0-9A-Z]+|NXF[1235]?|"
            r"NXT[12]|RAN(?:BP[0-9A-Z]*|GAP1|GEF1)?|SEC61[ABG]|SRP[0-9A-Z]+|"
            r"SSR[1-4]|TNPO[123]?|XPO[0-9A-Z]+)$",
            re.IGNORECASE,
        ),
    ),
    "Genome maintenance": (
        re.compile(
            r"^(?:ATM|ATR|BLM|BRCA[12]|CHEK[12]|FANC[A-Z0-9]+|MCM[2-9]|PCNA|"
            r"POLD[1-4]|POLE[2-4]?|RAD[0-9A-Z]+|RPA[1234]?|RFC[1-5]|RECQL[145]?|"
            r"TOPBP1|XRCC[0-9]+)$",
            re.IGNORECASE,
        ),
    ),
}

CLASS_PRIORITY = [
    "Ribosomal protein",
    "Ribosome biogenesis",
    "Translation",
    "RNA processing",
    "Transcription",
    "Genome maintenance",
    "Transport",
    "Development",
    "Signaling",
    "Other coding",
]

#Gene ontology biological process anchors
#Theses term are primarly used for classification
#   (BP children of these anchors will be classify in the same class)
CLASS_BP_ANCHORS = {
    "Ribosome biogenesis": {
        "GO:0042254",  # ribosome biogenesis
    },
    "Translation": {
        "GO:0006412",  # translation
    },
    "RNA processing": {
        "GO:0006396",  # RNA processing
    },
    "Transcription": {
        "GO:0006351",  # DNA-templated transcription
        "GO:0001172",  # RNA-templated transcription
    },
    "Transport": {
        "GO:0006810",  # transport
    },
    "Development": {
        "GO:0032502",  # developmental process
    },
    "Signaling": {
        "GO:0023052",  # signaling
        "GO:0007154",  # cell communication
    },
    "Genome maintenance": {
        "GO:0007049",  # cell cycle
        "GO:0006259",  #DNA metabolic process
    },
}

CLASS_MF_ANCHORS = {
    "Transcription": {
        "GO:0140110",  # DNA-binding transcription factor activity
    },
    "Genome maintenance": {
        "GO:0003678",  # DNA helicase activity
        "GO:0003684",  # damaged DNA binding
        "GO:0003697",  # single-stranded DNA binding
        "GO:0003887",  # DNA-directed DNA polymerase activity
        "GO:0017108",  # 5'-flap endonuclease activity
        "GO:0140658",  # ATP-dependent chromatin remodeler activity
    },
}
GO_UPPER_RELATIONSHIPS = {
    "part_of",
    "regulates",
    "positively_regulates",
    "negatively_regulates",
}

ENSEMBL_VERSION_RE = re.compile(r"^(ENS[A-Z0-9]+)\.\d+$", re.IGNORECASE)
INTERNAL_LOOKUP_COLUMN = "__host_gene_id_resolved"


@dataclass
class PredictionResult:
    host_type: str
    rule: str


def clean_string(value: object) -> str:
    """Return a stripped string while safely handling missing scalar values."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def resolve_ensembl_id(value: object) -> str:
    """Normalize an Ensembl identifier and drop transcript-style version suffixes."""
    text = clean_string(value)
    if not text or text.lower() in {"nan", "none"}:
        return ""
    if text == "Intergenic":
        return "Intergenic"
    match = ENSEMBL_VERSION_RE.match(text)
    if match:
        return match.group(1)
    return text


def split_pipe(value: object, prefix: str | None = None) -> list[str]:
    """Split a pipe-delimited field into unique cleaned items."""
    text = clean_string(value)
    if not text:
        return []
    result = []
    seen = set()
    for part in text.split("|"):
        part = part.strip()
        if not part:
            continue
        if prefix and not part.startswith(prefix):
            continue
        if part not in seen:
            seen.add(part)
            result.append(part)
    return result


def safe_literal_eval(value: object) -> object:
    """Parse JSON- or Python-like serialized values into native objects."""
    if isinstance(value, (list, dict)):
        return value
    text = clean_string(value)
    if not text:
        return None
    for parser in (json.loads, ast.literal_eval):
        try:
            return parser(text)
        except Exception:
            continue
    return None


def ensure_directory(path: Path) -> None:
    """Create the parent directory for a file path when it does not exist."""
    path.parent.mkdir(parents=True, exist_ok=True)


def download_if_needed(url: str, destination: Path, refresh: bool = False) -> Path:
    """Download a remote file unless a local cached copy can be reused."""
    ensure_directory(destination)
    if destination.exists() and not refresh:
        return destination
    request = Request(url, headers={"User-Agent": "Classification_genes/1.0"})
    with urlopen(request) as response, destination.open("wb") as handle:
        handle.write(response.read())
    return destination


def post_json(url: str, payload: dict[str, str]) -> object:
    """Send a form-encoded POST request and decode the JSON response."""
    data = urlencode(payload).encode("utf-8")
    request = Request(
        url,
        data=data,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "Classification_genes/1.0",
        },
    )
    with urlopen(request) as response:
        return json.loads(response.read().decode("utf-8"))


def chunked(values: list[str], size: int) -> Iterable[list[str]]:
    """Yield a list in fixed-size batches."""
    for start in range(0, len(values), size):
        yield values[start : start + size]


def serialize_json(value: object) -> str:
    """Serialize a nested value to JSON while keeping empty values blank."""
    if value in (None, "", [], {}):
        return ""
    return json.dumps(value, sort_keys=True)


def extract_aspect_entries(value: object, aspect: str | None = None) -> list[dict[str, object]]:
    """Extract GO annotation entries as a uniform list of dictionaries."""
    parsed = safe_literal_eval(value)
    if parsed is None:
        return []
    if aspect is not None and isinstance(parsed, dict):
        parsed = parsed.get(aspect)
    if isinstance(parsed, dict):
        return [parsed]
    if isinstance(parsed, list):
        return [entry for entry in parsed if isinstance(entry, dict)]
    return []


def extract_ensembl_entries(value: object) -> list[dict[str, object]]:
    """Extract Ensembl annotation entries as a uniform list of dictionaries."""
    parsed = safe_literal_eval(value)
    if isinstance(parsed, dict):
        return [parsed]
    if isinstance(parsed, list):
        return [entry for entry in parsed if isinstance(entry, dict)]
    return []


def extract_ensembl_type_of_gene(value: object, gene_id: object | None = None) -> str:
    """Return the best Ensembl type_of_gene value for one gene annotation."""
    gene_id_text = clean_string(gene_id)
    entries = extract_ensembl_entries(value)

    if gene_id_text:
        for entry in entries:
            if resolve_ensembl_id(entry.get("gene")) != gene_id_text:
                continue
            type_of_gene = clean_string(entry.get("type_of_gene"))
            if type_of_gene:
                return type_of_gene

    for entry in entries:
        type_of_gene = clean_string(entry.get("type_of_gene"))
        if type_of_gene:
            return type_of_gene
    return ""


def select_gene_biotype(ensembl_type_of_gene: object, fallback_biotype: object) -> str:
    """Prefer Ensembl biotype details when available because they better resolve lncRNAs."""
    return clean_string(ensembl_type_of_gene) or clean_string(fallback_biotype)


def is_lnc_biotype(value: object) -> bool:
    """Return True when a biotype explicitly marks a long non-coding RNA."""
    return clean_string(value).lower() in {"lncrna", "lincrna"}


def extract_go_ids_from_value(value: object) -> list[str]:
    """Collect unique GO identifiers from a serialized GO field."""
    ids: list[str] = []
    seen = set()
    for entry in extract_aspect_entries(value):
        go_id = entry.get("id")
        if isinstance(go_id, str) and go_id.startswith("GO:") and go_id not in seen:
            seen.add(go_id)
            ids.append(go_id)
    return ids


def extract_go_terms_from_value(value: object) -> list[str]:
    """Collect unique GO term names from a serialized GO field."""
    terms: list[str] = []
    seen = set()
    for entry in extract_aspect_entries(value):
        term = entry.get("term")
        if isinstance(term, str) and term and term not in seen:
            seen.add(term)
            terms.append(term)
    return terms


def first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate column that exists in a dataframe."""
    for column in candidates:
        if column in df.columns:
            return column
    return None


def default_annotation_path(input_path: Path) -> Path:
    """Return the default annotation file path stored next to this script."""
    return SCRIPT_DIR / f"{input_path.stem}.mygene_annotations.csv"


def normalize_annotation_table(annotation_df: pd.DataFrame) -> pd.DataFrame:
    """Standardize annotation columns into the schema expected by the classifier."""
    df = annotation_df.copy()

    gene_id_column = first_existing_column(
        df,
        ["gene_id", "query", "_id", "ensembl.gene", "ensembl_gene_id"],
    )
    if gene_id_column is None:
        raise ValueError(
            "Could not find the gene ID column in the annotation file. "
            "Please use a MyGene-style annotation file containing one of: "
            "gene_id, query, _id, ensembl.gene, ensembl_gene_id."
        )

    df["gene_id"] = df[gene_id_column].map(resolve_ensembl_id)

    if "ensembl" not in df.columns:
        df["ensembl"] = ""
    if "ensembl.type_of_gene" not in df.columns:
        df["ensembl.type_of_gene"] = ""

    missing_ensembl_biotype_mask = df["ensembl.type_of_gene"].map(clean_string).eq("")
    if missing_ensembl_biotype_mask.any():
        df.loc[missing_ensembl_biotype_mask, "ensembl.type_of_gene"] = df.loc[
            missing_ensembl_biotype_mask
        ].apply(
            lambda row: extract_ensembl_type_of_gene(row.get("ensembl"), row.get("gene_id")),
            axis=1,
        )

    symbol_column = first_existing_column(df, ["symbol", "gene_name", "Symbol"])
    if symbol_column is not None:
        df["symbol"] = df[symbol_column]
    elif "symbol" not in df.columns:
        df["symbol"] = ""

    gene_name_column = first_existing_column(df, ["gene_name", "symbol", "Symbol"])
    df["gene_name"] = df[gene_name_column] if gene_name_column is not None else ""

    long_name_column = first_existing_column(df, ["name", "Long_Name", "description"])
    df["name"] = df[long_name_column] if long_name_column is not None else ""

    fallback_biotype_column = first_existing_column(df, ["gene_biotype", "type_of_gene", "biotype"])
    fallback_biotype = df[fallback_biotype_column] if fallback_biotype_column is not None else ""
    df["gene_biotype"] = [
        select_gene_biotype(ensembl_type_of_gene, fallback_value)
        for ensembl_type_of_gene, fallback_value in zip(
            df["ensembl.type_of_gene"],
            fallback_biotype if isinstance(fallback_biotype, pd.Series) else [fallback_biotype] * len(df),
        )
    ]

    if "taxid" not in df.columns:
        df["taxid"] = pd.NA
    if "notfound" not in df.columns:
        df["notfound"] = False

    if "go" in df.columns:
        for aspect in ("BP", "MF", "CC"):
            column = f"go.{aspect}"
            if column not in df.columns:
                df[column] = df["go"].apply(
                    lambda value, current_aspect=aspect: serialize_json(
                        safe_literal_eval(value).get(current_aspect)
                        if isinstance(safe_literal_eval(value), dict)
                        else None
                    )
                )

    for aspect in ("BP", "MF", "CC"):
        raw_column = f"go.{aspect}"
        id_column = f"go.{aspect}.id"
        term_column = f"go.{aspect}.term"

        if raw_column not in df.columns:
            df[raw_column] = ""

        if id_column not in df.columns:
            df[id_column] = ""
        if term_column not in df.columns:
            df[term_column] = ""

        df[raw_column] = df[raw_column].astype("object")
        df[id_column] = df[id_column].astype("object")
        df[term_column] = df[term_column].astype("object")

        missing_id_mask = df[id_column].map(clean_string).eq("")
        if missing_id_mask.any():
            df.loc[missing_id_mask, id_column] = df.loc[missing_id_mask, raw_column].apply(
                lambda value: "|".join(extract_go_ids_from_value(value))
            )

        missing_term_mask = df[term_column].map(clean_string).eq("")
        if missing_term_mask.any():
            df.loc[missing_term_mask, term_column] = df.loc[missing_term_mask, raw_column].apply(
                lambda value: "|".join(extract_go_terms_from_value(value))
            )

    keep_columns = [
        "gene_id",
        "ensembl",
        "ensembl.type_of_gene",
        "symbol",
        "gene_name",
        "name",
        "gene_biotype",
        "taxid",
        "notfound",
        "go.BP",
        "go.BP.id",
        "go.BP.term",
        "go.MF",
        "go.MF.id",
        "go.MF.term",
        "go.CC",
        "go.CC.id",
        "go.CC.term",
    ]
    for column in keep_columns:
        if column not in df.columns:
            df[column] = ""

    normalized = df[keep_columns].copy()
    normalized = normalized[normalized["gene_id"].map(clean_string).ne("")]
    normalized = normalized.drop_duplicates(subset=["gene_id"], keep="first")
    return normalized


def build_annotation_rows_from_mygene(response_rows: list[dict[str, object]]) -> pd.DataFrame:
    """Convert raw MyGene API records into the normalized annotation table."""
    records = []
    for row in response_rows:
        go = row.get("go") if isinstance(row.get("go"), dict) else {}
        bp = go.get("BP")
        mf = go.get("MF")
        cc = go.get("CC")
        ensembl = row.get("ensembl")
        gene_id = resolve_ensembl_id(row.get("query") or row.get("_id"))
        ensembl_type_of_gene = extract_ensembl_type_of_gene(ensembl, gene_id)

        record = {
            "gene_id": gene_id,
            "ensembl": serialize_json(ensembl),
            "ensembl.type_of_gene": ensembl_type_of_gene,
            "symbol": row.get("symbol", ""),
            "gene_name": row.get("symbol", ""),
            "name": row.get("name", ""),
            "gene_biotype": select_gene_biotype(ensembl_type_of_gene, row.get("type_of_gene", "")),
            "taxid": row.get("taxid", pd.NA),
            "notfound": bool(row.get("notfound", False)),
            "go.BP": serialize_json(bp),
            "go.BP.id": "|".join(extract_go_ids_from_value(bp)),
            "go.BP.term": "|".join(extract_go_terms_from_value(bp)),
            "go.MF": serialize_json(mf),
            "go.MF.id": "|".join(extract_go_ids_from_value(mf)),
            "go.MF.term": "|".join(extract_go_terms_from_value(mf)),
            "go.CC": serialize_json(cc),
            "go.CC.id": "|".join(extract_go_ids_from_value(cc)),
            "go.CC.term": "|".join(extract_go_terms_from_value(cc)),
        }
        records.append(record)

    return normalize_annotation_table(pd.DataFrame(records))


def fetch_mygene_annotations(gene_ids: list[str], email: str | None = None) -> pd.DataFrame:
    """Download gene annotations from MyGene.info in batches."""
    if not gene_ids:
        return normalize_annotation_table(pd.DataFrame({"gene_id": []}))

    all_rows: list[dict[str, object]] = []
    total_batches = (len(gene_ids) + MYGENE_BATCH_SIZE - 1) // MYGENE_BATCH_SIZE
    for batch_index, batch in enumerate(chunked(gene_ids, MYGENE_BATCH_SIZE), start=1):
        print(f"Fetching MyGene annotations: batch {batch_index}/{total_batches} ({len(batch)} genes)")
        payload = {
            "ids": ",".join(batch),
            "fields": MYGENE_FIELDS,
        }
        if email:
            payload["email"] = email
        response = post_json(MYGENE_GENE_URL, payload)
        if not isinstance(response, list):
            raise RuntimeError("Unexpected MyGene response: expected a list of gene objects.")
        all_rows.extend(response)

    return build_annotation_rows_from_mygene(all_rows)


def load_or_update_annotation_file(
    gene_ids: list[str],
    annotation_path: Path,
    email: str | None = None,
    force_download: bool = False,
) -> pd.DataFrame:
    """Load one annotation file and update it only for genes still missing."""
    existing = normalize_annotation_table(pd.DataFrame({"gene_id": []}))
    annotation_has_ensembl_biotype = False
    should_save_annotations = force_download or not annotation_path.exists()
    missing_gene_ids = list(gene_ids)

    if annotation_path.exists() and not force_download:
        print(f"Using annotation file {annotation_path}")
        annotation_raw = pd.read_csv(annotation_path, low_memory=False)
        annotation_has_ensembl_biotype = (
            "ensembl.type_of_gene" in annotation_raw.columns or "ensembl" in annotation_raw.columns
        )
        existing = normalize_annotation_table(annotation_raw)
        cached_gene_ids = set(existing["gene_id"]) if not existing.empty else set()
        missing_gene_ids = sorted(set(gene_ids) - cached_gene_ids)
        if not existing.empty and not annotation_has_ensembl_biotype:
            print("Annotation file uses the legacy schema; refreshing it to recover Ensembl biotypes")
            missing_gene_ids = list(gene_ids)
            existing = normalize_annotation_table(pd.DataFrame({"gene_id": []}))
            should_save_annotations = True
    elif force_download:
        print(f"Force-downloading annotations into {annotation_path}")
    else:
        print(f"Creating annotation file {annotation_path}")

    if missing_gene_ids:
        fetched = fetch_mygene_annotations(sorted(set(missing_gene_ids)), email=email)
        combined = pd.concat([existing, fetched], ignore_index=True)
        combined = combined.drop_duplicates(subset=["gene_id"], keep="last")
        ensure_directory(annotation_path)
        combined.to_csv(annotation_path, index=False)
        print(f"Saved annotations to {annotation_path}")
        return normalize_annotation_table(combined)

    if should_save_annotations and existing.empty:
        ensure_directory(annotation_path)
        existing.to_csv(annotation_path, index=False)
        print(f"Saved annotations to {annotation_path}")

    print("Annotation file already covers all host genes")
    return existing


def get_annotation_ids(row: pd.Series, aspect: str) -> list[str]:
    """Return GO IDs for one aspect from a normalized annotation row."""
    ids = split_pipe(row.get(f"go.{aspect}.id"), prefix="GO:")
    if ids:
        return ids
    return extract_go_ids_from_value(row.get(f"go.{aspect}"))


def get_annotation_terms(row: pd.Series, aspect: str) -> list[str]:
    """Return GO terms for one aspect from a normalized annotation row."""
    terms = split_pipe(row.get(f"go.{aspect}.term"))
    if terms:
        return terms
    return extract_go_terms_from_value(row.get(f"go.{aspect}"))


class GOGroupClassifier:
    def __init__(self, godag: GODag) -> None:
        """Store the GO graph used to traverse GO ancestors during classification."""
        self.godag = godag

    def expand_go_ids(self, go_ids: list[str], namespace: str) -> set[str]:
        """Expand GO IDs through is_a plus selected upper GO relationships."""
        expanded: set[str] = set()
        root_by_namespace = {
            "biological_process": "GO:0008150",
            "molecular_function": "GO:0003674",
            "cellular_component": "GO:0005575",
        }
        for go_id in go_ids:
            expanded.update(self.collect_upper_go_ids(go_id, namespace))
        root_id = root_by_namespace.get(namespace)
        if root_id:
            expanded.discard(root_id)
        return expanded

    def ribosomal_protein_rule(
        self,
        gene_name: object,
        long_name: object,
        expanded_mf: set[str],
    ) -> bool:
        """Detect true structural ribosomal proteins with strict signals only."""
        if expanded_mf & RIBOSOMAL_PROTEIN_MF:
            return True
        if isinstance(gene_name, str) and STRUCTURAL_RIBOSOMAL_SYMBOL_RE.match(gene_name.strip()):
            return True
        if isinstance(long_name, str) and any(
            pattern.match(long_name.strip()) for pattern in STRUCTURAL_RIBOSOMAL_NAME_PATTERNS
        ):
            return True
        return False

    def gene_name_fallback(self, gene_name: object) -> str | None:
        """Classify by strict symbol pattern when GO evidence is missing or uninformative."""
        if not isinstance(gene_name, str):
            return None
        gene_name = gene_name.strip()
        if not gene_name:
            return None
        for host_type in CLASS_PRIORITY:
            patterns = NO_GO_GENE_NAME_KEYWORDS.get(host_type, ())
            if any(pattern.match(gene_name) for pattern in patterns):
                return host_type
        return None

    def collect_upper_go_ids(self, go_id: str, namespace: str) -> set[str]:
        """Collect is_a, part_of, and regulation-linked ancestors in one namespace."""
        if go_id not in self.godag:
            return set()
        start_term = self.godag[go_id]
        if start_term.namespace != namespace:
            return set()

        collected: set[str] = set()
        stack = [go_id]
        while stack:
            current_id = stack.pop()
            if current_id in collected or current_id not in self.godag:
                continue
            current_term = self.godag[current_id]
            if current_term.namespace != namespace:
                continue
            collected.add(current_id)
            stack.extend(parent.item_id for parent in current_term.parents)
            for relationship_type, related_terms in getattr(current_term, "relationship", {}).items():
                if relationship_type not in GO_UPPER_RELATIONSHIPS:
                    continue
                stack.extend(term.item_id for term in related_terms)
        return collected

    def build_term_ancestor_map(self, go_ids: list[str], namespace: str) -> dict[str, set[str]]:
        """Map each GO term to itself plus upper linked terms within one namespace."""
        term_ancestors: dict[str, set[str]] = {}
        for go_id in go_ids:
            ancestors = self.collect_upper_go_ids(go_id, namespace)
            if not ancestors:
                continue
            term_ancestors[go_id] = ancestors
        return term_ancestors

    def count_anchor_support(
        self,
        term_ancestors: dict[str, set[str]],
        anchor_ids: set[str],
    ) -> int:
        """Count GO terms whose ancestry reaches at least one class anchor."""
        if not anchor_ids:
            return 0
        return sum(1 for ancestors in term_ancestors.values() if ancestors & anchor_ids)

    def class_term_counts(self, bp_ids: list[str], mf_ids: list[str]) -> dict[str, int]:
        """Count supporting GO terms per class across BP and MF annotations."""
        bp_term_ancestors = self.build_term_ancestor_map(bp_ids, "biological_process")
        mf_term_ancestors = self.build_term_ancestor_map(mf_ids, "molecular_function")
        term_counts: dict[str, int] = {}
        for host_type in CLASS_PRIORITY:
            term_count = 0
            term_count += self.count_anchor_support(
                bp_term_ancestors,
                CLASS_BP_ANCHORS.get(host_type, set()),
            )
            term_count += self.count_anchor_support(
                mf_term_ancestors,
                CLASS_MF_ANCHORS.get(host_type, set()),
            )
            if term_count > 0:
                term_counts[host_type] = term_count
        return term_counts

    def predict(
        self,
        gene_id: object,
        gene_name: object,
        long_name: object,
        biotype: object,
        bp_ids: list[str],
        mf_ids: list[str],
        cc_ids: list[str],
    ) -> PredictionResult:
        """Assign a host-gene class from gene metadata and GO annotations."""
        gene_id_text = clean_string(gene_id)
        if not gene_id_text or gene_id_text == "Intergenic":
            return PredictionResult("Intergenic", "intergenic")

        if is_lnc_biotype(biotype):
            return PredictionResult("lncRNA", "biotype")

        expanded_mf = self.expand_go_ids(mf_ids, "molecular_function")
        has_go_annotation = bool(bp_ids or mf_ids or cc_ids)

        if self.ribosomal_protein_rule(gene_name, long_name, expanded_mf):
            return PredictionResult("Ribosomal protein", "ribosomal_protein_rule")

        term_counts = self.class_term_counts(bp_ids, mf_ids)
        if not term_counts:
            fallback_host_type = self.gene_name_fallback(gene_name)
            if fallback_host_type is not None:
                fallback_rule = (
                    "gene_name_keyword_fallback" if not has_go_annotation else "gene_name_keyword_rescue"
                )
                return PredictionResult(fallback_host_type, fallback_rule)
            return PredictionResult("Other coding", "no_anchor_match")

        ranked = sorted(
            term_counts.items(),
            key=lambda item: (
                -item[1],
                CLASS_PRIORITY.index(item[0]) if item[0] in CLASS_PRIORITY else len(CLASS_PRIORITY),
                item[0],
            ),
        )
        host_type, _ = ranked[0]
        return PredictionResult(host_type, "go_term_count")


def classify_annotation_row(row: pd.Series, classifier: GOGroupClassifier) -> PredictionResult:
    """Classify one normalized annotation row and preserve notfound status."""
    prediction = classifier.predict(
        gene_id=row.get("gene_id"),
        gene_name=row.get("gene_name", row.get("symbol")),
        long_name=row.get("name"),
        biotype=row.get("gene_biotype"),
        bp_ids=get_annotation_ids(row, "BP"),
        mf_ids=get_annotation_ids(row, "MF"),
        cc_ids=get_annotation_ids(row, "CC"),
    )
    if bool(row.get("notfound")) and prediction.rule == "no_anchor_match":
        return replace(prediction, rule="annotation_notfound")
    return prediction


def build_host_classification_table(
    input_df: pd.DataFrame,
    host_gene_column: str,
    annotation_df: pd.DataFrame,
    classifier: GOGroupClassifier,
) -> pd.DataFrame:
    """Build one row per unique host gene with classification and annotation fields."""
    unique_hosts = (
        input_df[[host_gene_column, INTERNAL_LOOKUP_COLUMN]]
        .drop_duplicates()
        .rename(columns={host_gene_column: "host_gene_id_original"})
    )

    merged = unique_hosts.merge(
        annotation_df,
        left_on=INTERNAL_LOOKUP_COLUMN,
        right_on="gene_id",
        how="left",
    )

    predictions = [classify_annotation_row(row, classifier) for _, row in merged.iterrows()]

    merged["host_gene_id_resolved"] = merged[INTERNAL_LOOKUP_COLUMN]
    merged["host_gene_symbol"] = merged["symbol"]
    merged["host_gene_description"] = merged["name"]
    merged["host_gene_biotype"] = merged["gene_biotype"]
    merged["host_gene_taxid"] = merged["taxid"]
    merged["host_gene_annotation_notfound"] = (
        merged["notfound"].where(merged["notfound"].notna(), False).astype(bool)
    )
    merged["host_gene_type"] = [prediction.host_type for prediction in predictions]
    merged["host_gene_type_rule"] = [prediction.rule for prediction in predictions]
    merged["host_gene_go_bp_ids"] = merged["go.BP.id"]
    merged["host_gene_go_bp_terms"] = merged["go.BP.term"]
    merged["host_gene_go_mf_ids"] = merged["go.MF.id"]
    merged["host_gene_go_mf_terms"] = merged["go.MF.term"]
    merged["host_gene_go_cc_ids"] = merged["go.CC.id"]
    merged["host_gene_go_cc_terms"] = merged["go.CC.term"]

    output_columns = [
        "host_gene_id_original",
        "host_gene_id_resolved",
        "host_gene_symbol",
        "host_gene_description",
        "host_gene_biotype",
        "host_gene_taxid",
        "host_gene_annotation_notfound",
        "host_gene_type",
        "host_gene_type_rule",
        "host_gene_go_bp_ids",
        "host_gene_go_bp_terms",
        "host_gene_go_mf_ids",
        "host_gene_go_mf_terms",
        "host_gene_go_cc_ids",
        "host_gene_go_cc_terms",
    ]
    return merged[output_columns].rename(columns={"host_gene_id_original": host_gene_column})


def merge_classification_into_input(
    input_df: pd.DataFrame,
    host_gene_column: str,
    host_classification_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge host-gene classifications back into the original input table."""
    output_columns = [
        "host_gene_id_resolved",
        "host_gene_symbol",
        "host_gene_description",
        "host_gene_biotype",
        "host_gene_taxid",
        "host_gene_annotation_notfound",
        "host_gene_type",
        "host_gene_type_rule",
        "host_gene_go_bp_ids",
        "host_gene_go_bp_terms",
        "host_gene_go_mf_ids",
        "host_gene_go_mf_terms",
        "host_gene_go_cc_ids",
        "host_gene_go_cc_terms",
    ]
    output_df = input_df.drop(columns=output_columns, errors="ignore").merge(
        host_classification_df,
        on=host_gene_column,
        how="left",
    )
    return output_df.drop(columns=[INTERNAL_LOOKUP_COLUMN], errors="ignore")


def build_argument_parser() -> argparse.ArgumentParser:
    """Define the command-line interface for the public classifier."""
    parser = argparse.ArgumentParser(
        description="Classify host genes from a CSV table using GO annotations and MyGene.info."
    )
    parser.add_argument("input_csv", type=Path, help="Input CSV file to classify.")
    parser.add_argument(
        "--host-gene-column",
        required=True,
        help="Column in the input CSV containing host gene Ensembl IDs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output CSV path. When omitted, the input CSV is overwritten in place.",
    )
    parser.add_argument(
        "--annotation-file",
        type=Path,
        help=(
            "Annotation CSV to reuse and update. "
            "Default: <input_stem>.mygene_annotations.csv next to this script."
        ),
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Ignore existing annotations in the annotation file and download them again.",
    )
    parser.add_argument(
        "--refresh-go",
        action="store_true",
        help="Re-download go-basic.obo even if it already exists next to this script.",
    )
    parser.add_argument(
        "--email",
        help="Optional email sent to MyGene.info, as encouraged in their documentation.",
    )
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()

    input_path = args.input_csv.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    input_df = pd.read_csv(input_path, low_memory=False).copy()
    if args.host_gene_column not in input_df.columns:
        raise ValueError(
            f"Column '{args.host_gene_column}' was not found in {input_path.name}. "
            f"Available columns: {', '.join(input_df.columns)}"
        )

    input_df[INTERNAL_LOOKUP_COLUMN] = input_df[args.host_gene_column].map(resolve_ensembl_id)
    gene_ids = sorted(
        {
            gene_id
            for gene_id in input_df[INTERNAL_LOOKUP_COLUMN].dropna().astype(str)
            if gene_id and gene_id != "Intergenic"
        }
    )
    print(f"Found {len(gene_ids):,} unique host gene IDs to classify")

    go_basic_path = download_if_needed(GO_BASIC_URL, GO_BASIC_PATH, refresh=args.refresh_go)
    print(f"Loading GO ontology from {go_basic_path}")
    godag = GODag(str(go_basic_path), optional_attrs={"relationship"})
    classifier = GOGroupClassifier(godag)

    annotation_file = (
        args.annotation_file.resolve()
        if args.annotation_file
        else default_annotation_path(input_path)
    )
    annotation_df = load_or_update_annotation_file(
        gene_ids=gene_ids,
        annotation_path=annotation_file,
        email=args.email,
        force_download=args.force_download,
    )

    host_classification_df = build_host_classification_table(
        input_df=input_df,
        host_gene_column=args.host_gene_column,
        annotation_df=annotation_df,
        classifier=classifier,
    )
    output_df = merge_classification_into_input(
        input_df=input_df,
        host_gene_column=args.host_gene_column,
        host_classification_df=host_classification_df,
    )

    output_path = args.output.resolve() if args.output else input_path
    output_df.to_csv(output_path, index=False)

    print(f"Wrote merged output to {output_path}")
    print("\nHost gene type counts:")
    print(host_classification_df["host_gene_type"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
