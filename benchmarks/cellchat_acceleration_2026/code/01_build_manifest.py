#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path


KEEP_EXT = {
    ".rds": "rds",
    ".h5ad": "h5ad",
    ".h5": "h5",
    ".hdf5": "h5",
    ".mtx": "mtx",
    ".csv": "csv",
    ".tsv": "tsv",
}


def human_bytes(n):
    units = ["B", "KB", "MB", "GB", "TB"]
    x = float(n)
    for u in units:
        if x < 1024 or u == units[-1]:
            return f"{x:.1f}{u}"
        x /= 1024


def infer_dataset_root(path, data_root):
    rel = path.relative_to(data_root)
    parts = rel.parts
    return parts[0] if parts else "."


def sample_id(path):
    name = path.name
    for suffix in [".rds", ".h5ad", ".hdf5", ".h5", ".mtx", ".csv", ".tsv"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def iter_files(data_root):
    for root, dirs, files in os.walk(data_root):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", "cellchat_result"}]
        for fn in files:
            p = Path(root) / fn
            ext = p.suffix.lower()
            if ext in KEEP_EXT or "cellchat" in fn.lower():
                try:
                    st = p.stat()
                except OSError:
                    continue
                yield p, st.st_size, KEEP_EXT.get(ext, "other")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", default="data")
    ap.add_argument("--out", default="results/data_manifest.tsv")
    ap.add_argument("--candidates", default="results/dataset_candidates.tsv")
    args = ap.parse_args()

    data_root = Path(args.data_root).resolve()
    rows = []
    for p, size, kind in iter_files(data_root):
        rows.append(
            {
                "dataset_root": infer_dataset_root(p, data_root),
                "sample_id": sample_id(p),
                "kind": kind,
                "bytes": size,
                "human_size": human_bytes(size),
                "path": str(p),
                "parent": str(p.parent),
            }
        )
    rows.sort(key=lambda r: (-int(r["bytes"]), r["path"]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [
            "dataset_root", "sample_id", "kind", "bytes", "human_size", "path", "parent"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    candidate_rows = []
    for r in rows:
        if r["kind"] == "rds":
            priority = 1
        elif r["kind"] in {"h5ad", "h5"}:
            priority = 2
        elif r["kind"] == "mtx":
            priority = 3
        else:
            continue
        candidate_rows.append({**r, "priority": priority})
    candidate_rows.sort(key=lambda r: (int(r["priority"]), -int(r["bytes"]), r["path"]))

    cand = Path(args.candidates)
    cand.parent.mkdir(parents=True, exist_ok=True)
    with cand.open("w", newline="") as fh:
        fieldnames = ["priority", "dataset_root", "sample_id", "kind", "bytes", "human_size", "path", "parent"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    print(f"manifest_rows={len(rows)} out={out}")
    print(f"candidate_rows={len(candidate_rows)} out={cand}")


if __name__ == "__main__":
    main()

