#!/usr/bin/env python3
"""Substitute the four path placeholders throughout the repository.

The scripts refer to data through ${PROJECT_ROOT}, ${DATA_ROOT}, ${WORK_ROOT}
and ${HOME}. These are written as plain text inside ordinary string literals,
so exporting shell variables does not reach them: Python and R both read
Path("${PROJECT_ROOT}/x") as a directory literally named "${PROJECT_ROOT}".
This script rewrites them in place, once, to the locations you give it.

Run it from anywhere:

  python config/set_paths.py --project-root /data/luad \\
                            --data-root /data/raw/LUAD \\
                            --work-root /data/work

Nothing is written until you add --apply; without it the script only reports
what it would change. Keep a pristine copy of the repository (or use git) if
you want to be able to point the code at a different machine later.

Paths are written with forward slashes, which both Python and R accept on
Windows as well.
"""
import argparse
import re
import sys
from pathlib import Path

PLACEHOLDERS = ("PROJECT_ROOT", "DATA_ROOT", "WORK_ROOT", "HOME")
EXTS = {".py", ".R", ".r", ".sh", ".md", ".yml", ".yaml", ".txt", ".tsv"}
SELF = "set_paths.py"


def repo_root():
    return Path(__file__).resolve().parent.parent


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-root", help="repository plus its results/ tree")
    ap.add_argument("--data-root", help="raw downloads: GEO, ArrayExpress, Visium")
    ap.add_argument("--work-root", help="large intermediates and figure tables")
    ap.add_argument("--home", help="user home; only a few helper paths use it")
    ap.add_argument("--apply", action="store_true", help="write the changes")
    args = ap.parse_args()

    values = {
        "PROJECT_ROOT": args.project_root,
        "DATA_ROOT": args.data_root,
        "WORK_ROOT": args.work_root,
        "HOME": args.home or str(Path.home()),
    }
    missing = [k for k, v in values.items() if not v]
    if missing:
        sys.exit(f"missing: {', '.join('--' + m.lower().replace('_', '-') for m in missing)}")
    values = {k: v.replace("\\", "/").rstrip("/") for k, v in values.items()}

    root = repo_root()
    pattern = re.compile(r"\$\{(" + "|".join(PLACEHOLDERS) + r")\}")
    touched, total = 0, 0
    per_key = {k: 0 for k in PLACEHOLDERS}

    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix not in EXTS or path.name == SELF:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        hits = pattern.findall(text)
        if not hits:
            continue
        touched += 1
        total += len(hits)
        for h in hits:
            per_key[h] += 1
        if args.apply:
            path.write_text(pattern.sub(lambda m: values[m.group(1)], text),
                            encoding="utf-8")

    print(f"{'rewrote' if args.apply else 'would rewrite'} {total} placeholders "
          f"in {touched} files")
    for k, n in per_key.items():
        print(f"  ${{{k}}}  {n:>4}  ->  {values[k]}")
    if not args.apply:
        print("\nnothing written; add --apply to make the change")
        return

    left = 0
    for path in root.rglob("*"):
        if path.is_file() and path.suffix in EXTS and path.name != SELF:
            left += len(pattern.findall(
                path.read_text(encoding="utf-8", errors="replace")))
    print(f"\nplaceholders remaining: {left}")
    if left:
        sys.exit("some placeholders were not substituted")


if __name__ == "__main__":
    main()
