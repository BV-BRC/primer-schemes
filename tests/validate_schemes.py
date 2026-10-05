#!/usr/bin/env python3
"""Validate the primer schemes in this repository.

Checks shipped scheme data against the published upstream artifacts and against
internal consistency rules. Pure data validation: no alignment tools, no network.

    python3 tests/validate_schemes.py            # validate everything
    python3 tests/validate_schemes.py -v         # list every check

Exit status is 0 if all checks pass, 1 otherwise.
"""

import argparse
import collections
import hashlib
import json
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Published md5s from quick-lab/primerschemes index.json. A scheme listed here is
# shipped verbatim and must stay byte-identical to upstream.
UPSTREAM_MD5 = {
    "nCoV-2019/V5.3.2/SARS-CoV-2.primer.bed": "46a3aeddc452678bd0cd33efc1c9558f",
    "nCoV-2019/V5.3.2/SARS-CoV-2.reference.fasta": "7f8995394dfc7d5ffeb9fe8322ade58c",
    "nCoV-2019/V5.4.2/SARS-CoV-2.primer.bed": "e7897c8be8e488836ce7777e7709ea53",
    "nCoV-2019/V5.4.2/SARS-CoV-2.reference.fasta": "d11d06b5d1eb1d85c69e341c3c026e08",
}

# Primers that deliberately do not match the reference, because they target a drifted
# lineage. Keyed by scheme dir -> {primer name: (0-based pos, reference base, primer base)}.
# The check below requires an EXACT match to this set, so a dropped spike-in fails just
# as loudly as an unexpected mismatch.
EXPECTED_MISMATCHES = {
    "nCoV-2019/V5.3.2": {
        "SARS-CoV-2_84_RIGHT_2": (26059, "C", "T"),
    },
    "nCoV-2019/V5.4.2": {
        "SARS-CoV-2_11_RIGHT_1": (3564, "T", "C"),
        "SARS-CoV-2_69_RIGHT_2": (21710, "C", "T"),
        "SARS-CoV-2_69_RIGHT_1": (21717, "G", "T"),
        "SARS-CoV-2_70_RIGHT_1": (21940, "G", "T"),
        "SARS-CoV-2_74_LEFT_1": (22769, "G", "A"),
        "SARS-CoV-2_84_RIGHT_2": (26059, "C", "T"),
    },
}

# v5.3.2 -> v5.4.2 is exactly these five spike-ins and nothing else.
EXPECTED_DELTA = {
    ("SARS-CoV-2_11_RIGHT_1", 3560, 3584, "1", "-"),
    ("SARS-CoV-2_69_RIGHT_1", 21696, 21722, "1", "-"),
    ("SARS-CoV-2_69_RIGHT_2", 21696, 21722, "1", "-"),
    ("SARS-CoV-2_70_RIGHT_1", 21927, 21960, "2", "-"),
    ("SARS-CoV-2_74_LEFT_1", 22742, 22774, "2", "+"),
}

# Schemes to run the full structural/topology suite over. Others in this repo predate
# these conventions (5-column beds, no sequence column) and are left alone.
FULL_CHECK = ["nCoV-2019/V5.3.2", "nCoV-2019/V5.4.2"]

PRIMER_RE = re.compile(r"^(?P<scheme>.+)_(?P<amp>\d+)_(?P<dir>LEFT|RIGHT)(?:_(?P<alt>\S+))?$")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


class Checker:
    def __init__(self, verbose):
        self.verbose = verbose
        self.failures = []
        self.count = 0

    def check(self, ok, label, detail=""):
        self.count += 1
        if ok:
            if self.verbose:
                print(f"  ok   {label}")
        else:
            self.failures.append(f"{label}{': ' + detail if detail else ''}")
            print(f"  FAIL {label}" + (f"\n         {detail}" if detail else ""))
        return ok


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def load_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    raise ValueError(f"{path}: expected a single record")
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq).upper()


def load_bed(path):
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            raw = line.rstrip("\n")
            if not raw.strip() or raw.startswith("#"):
                continue
            rows.append((lineno, raw, raw.split("\t")))
    return rows


def check_scheme(c, scheme_dir):
    """Structural, sequence and topology checks for one scheme directory."""
    print(f"\n{scheme_dir}")
    base = os.path.join(REPO, scheme_dir)
    bed_path = os.path.join(base, "SARS-CoV-2.primer.bed")
    ref_path = os.path.join(base, "SARS-CoV-2.reference.fasta")
    scheme_path = os.path.join(base, "SARS-CoV-2.scheme.bed")

    for p in (bed_path, ref_path, scheme_path):
        if not c.check(os.path.isfile(p), f"{os.path.basename(p)} exists"):
            return

    refname, ref = load_fasta(ref_path)
    rows = load_bed(bed_path)

    # --- file hygiene -----------------------------------------------------------
    with open(bed_path, "rb") as fh:
        blob = fh.read()
    c.check(b"\r" not in blob, "primer.bed has no CR characters")
    c.check(blob.endswith(b"\n"), "primer.bed ends with a newline")

    with open(scheme_path, "rb") as fh:
        sblob = fh.read()
    c.check(b"\r" not in sblob, "scheme.bed has no CR characters")
    c.check(sblob.endswith(b"\n"), "scheme.bed ends with a newline")

    # scheme.bed must be columns 1-6 of primer.bed; rewrite-primers.pl reads it and
    # dies if it is absent, so it must stay in step with primer.bed.
    expect_scheme = "".join("\t".join(f[:6]) + "\n" for _, _, f in rows)
    c.check(sblob.decode() == expect_scheme, "scheme.bed == columns 1-6 of primer.bed")

    # --- per-row structure ------------------------------------------------------
    bad_cols = [ln for ln, _, f in rows if len(f) != 7]
    c.check(not bad_cols, "every primer.bed row has 7 columns",
            f"rows {bad_cols[:5]}")

    names, amps = [], collections.defaultdict(lambda: {"LEFT": [], "RIGHT": []})
    bad_coord, bad_chrom, bad_len, bad_name, bad_pool, bad_strand = [], [], [], [], [], []
    for ln, _, f in rows:
        if len(f) != 7:
            continue
        chrom, start, end, name, pool, strand, seq = f
        start, end = int(start), int(end)
        names.append(name)
        if chrom != refname:
            bad_chrom.append((ln, chrom))
        if not (0 <= start < end <= len(ref)):
            bad_coord.append((ln, start, end))
        if len(seq) != end - start:
            bad_len.append((ln, name, len(seq), end - start))
        if pool not in ("1", "2"):
            bad_pool.append((ln, pool))
        if strand not in ("+", "-"):
            bad_strand.append((ln, strand))
        m = PRIMER_RE.match(name)
        if not m:
            bad_name.append((ln, name))
        else:
            amps[int(m.group("amp"))][m.group("dir")].append(
                dict(start=start, end=end, name=name, pool=pool, strand=strand, seq=seq))

    c.check(not bad_chrom, f"every chrom == reference id ({refname})", str(bad_chrom[:5]))
    c.check(not bad_coord, f"0 <= start < end <= {len(ref)}", str(bad_coord[:5]))
    c.check(not bad_len, "len(seq) == end - start", str(bad_len[:5]))
    c.check(not bad_pool, "pool is 1 or 2", str(bad_pool[:5]))
    c.check(not bad_strand, "strand is + or -", str(bad_strand[:5]))
    c.check(not bad_name, "every primer name parses as <scheme>_<amp>_<LEFT|RIGHT>[_alt]",
            str(bad_name[:5]))

    dupes = [n for n, k in collections.Counter(names).items() if k > 1]
    c.check(not dupes, "primer names are unique", str(dupes[:5]))

    # --- sequence vs reference coordinate ---------------------------------------
    # The strongest check here: catches off-by-one coordinates, wrong strand and
    # truncated sequence, none of which an md5 catches once a file is regenerated.
    found = {}
    for ln, _, f in rows:
        if len(f) != 7:
            continue
        _, start, end, name, _, strand, seq = f
        start, end = int(start), int(end)
        expect = ref[start:end]
        observed = revcomp(seq) if strand == "-" else seq
        if len(observed) != len(expect):
            continue  # already reported by the length check
        diffs = [(start + i, a, b) for i, (a, b) in enumerate(zip(expect, observed)) if a != b]
        if diffs:
            found[name] = diffs

    expected = EXPECTED_MISMATCHES.get(scheme_dir, {})
    single = {n: d[0] for n, d in found.items() if len(d) == 1}
    multi = {n: d for n, d in found.items() if len(d) > 1}
    c.check(not multi, "no primer differs from the reference at more than one position",
            str({k: v for k, v in list(multi.items())[:3]}))
    c.check(single == expected,
            f"exactly {len(expected)} primers differ from the reference, as pinned",
            f"unexpected={ {k: v for k, v in single.items() if expected.get(k) != v} } "
            f"missing={ {k: v for k, v in expected.items() if single.get(k) != v} }")

    # --- amplicon topology ------------------------------------------------------
    ids = sorted(amps)
    c.check(ids == list(range(1, max(ids) + 1)) if ids else False,
            f"amplicons numbered 1..{max(ids) if ids else 0} with no gaps")
    c.check(all(amps[i]["LEFT"] and amps[i]["RIGHT"] for i in ids),
            "every amplicon has at least one LEFT and one RIGHT primer",
            str([i for i in ids if not (amps[i]["LEFT"] and amps[i]["RIGHT"])][:5]))
    c.check(all(all(p["strand"] == "+" for p in amps[i]["LEFT"]) and
                all(p["strand"] == "-" for p in amps[i]["RIGHT"]) for i in ids),
            "LEFT primers are + and RIGHT primers are -")
    c.check(all({p["pool"] for p in amps[i]["LEFT"] + amps[i]["RIGHT"]} ==
                {"1" if i % 2 else "2"} for i in ids),
            "pool alternates with amplicon parity (odd->1, even->2)",
            str([i for i in ids
                 if {p["pool"] for p in amps[i]["LEFT"] + amps[i]["RIGHT"]} !=
                 {"1" if i % 2 else "2"}][:5]))
    c.check(all(max(p["end"] for p in amps[i]["LEFT"]) <
                min(p["start"] for p in amps[i]["RIGHT"]) for i in ids),
            "every amplicon has a non-empty insert (LEFT.end < RIGHT.start)")

    # Overlapping tiling: amplicon N+1 must start before amplicon N ends, otherwise a
    # variant under a primer could not be recovered from a neighbour.
    gaps = [(i, i + 1) for i in ids[:-1]
            if min(p["start"] for p in amps[i + 1]["LEFT"]) >=
            max(p["end"] for p in amps[i]["RIGHT"])]
    c.check(not gaps, "consecutive amplicons overlap (no coverage gap)", str(gaps[:5]))

    # --- the transform sars2-onecodex.pl applies --------------------------------
    # sars2-onecodex.pl:135 derives strand from the primer NAME, not the bed column:
    #   $x[3] =~ m/LEFT|(F$)/ ? "+" : "-"
    # If a scheme ever names primers such that this disagrees with the bed, reads would
    # be trimmed on the wrong end. Pin it here, next to the data.
    wrong = [f[3] for _, _, f in rows if len(f) == 7
             and ("+" if re.search(r"LEFT|(F$)", f[3]) else "-") != f[5]]
    c.check(not wrong, "/LEFT|(F$)/ on the primer name reproduces the bed strand column",
            str(wrong[:5]))

    return rows


def check_delta(c, old_rows, new_rows):
    print("\nnCoV-2019/V5.3.2 -> nCoV-2019/V5.4.2")
    def key(f):
        return (f[3], int(f[1]), int(f[2]), f[4], f[5])
    old = {key(f) for _, _, f in old_rows if len(f) == 7}
    new = {key(f) for _, _, f in new_rows if len(f) == 7}
    c.check(new - old == EXPECTED_DELTA, "added primers are exactly the five spike-ins",
            f"unexpected={sorted(new - old - EXPECTED_DELTA)} "
            f"missing={sorted(EXPECTED_DELTA - (new - old))}")
    c.check(not (old - new), "no primer was removed", str(sorted(old - new)[:5]))

    # The spike-ins are alternative oligos for binding sites that already exist, so they
    # must not introduce a new interval. This is why adding v5.4.2 cannot change trimming.
    def intervals(rows):
        return {(int(f[1]), int(f[2]), f[5]) for _, _, f in rows if len(f) == 7}
    c.check(intervals(old_rows) == intervals(new_rows),
            "no new trim interval: every spike-in duplicates an existing site",
            str(sorted(intervals(new_rows) - intervals(old_rows))))


def check_manifest(c):
    print("\nbvbrc_manifest.json")
    path = os.path.join(REPO, "bvbrc_manifest.json")
    try:
        with open(path) as fh:
            manifest = json.load(fh)
    except Exception as exc:
        c.check(False, "manifest is valid JSON", str(exc))
        return

    c.check(True, "manifest is valid JSON")

    # Duplicate keys parse as last-wins and silently mask earlier bodies. Re-read with a
    # hook that sees every pair so a reintroduced duplicate is caught.
    dupes = []

    def hook(pairs):
        seen = set()
        for k, _ in pairs:
            if k in seen:
                dupes.append(k)
            seen.add(k)
        return dict(pairs)

    with open(path) as fh:
        json.load(fh, object_pairs_hook=hook)
    c.check(not dupes, "no duplicate keys anywhere in the manifest", str(sorted(set(dupes))))

    org = manifest["organisms"][0]
    primers = org["primers"]
    artic = primers.get("ARTIC")
    c.check(artic is not None, "ARTIC entry present")
    if artic:
        versions = [s["version"] for s in artic["schemes"]]
        c.check("V5.4.2" in versions, "V5.4.2 is registered", str(versions))
        # sars2-onecodex.pl:113 uses $schemes->[-1] when no primer_version is given.
        c.check(versions[-1] == "V5.4.2",
                "V5.4.2 is last, so it is the default for ARTIC jobs", str(versions))

    # Every declared file must exist. This is what catches a stale or copy-pasted entry.
    missing = []
    for kit, body in sorted(primers.items()):
        for s in body["schemes"]:
            for field in ("primers", "reference"):
                rel = os.path.join(body["path"], s["version"], s[field])
                if not os.path.isfile(os.path.join(REPO, rel)):
                    missing.append(f"{kit}/{s['version']}: {rel}")
    c.check(not missing, "every path declared in the manifest exists on disk",
            "\n         ".join(missing))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="store_true", help="list passing checks too")
    args = ap.parse_args()

    c = Checker(args.verbose)

    print("upstream md5")
    for rel, want in sorted(UPSTREAM_MD5.items()):
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            c.check(False, f"{rel} exists")
            continue
        got = md5(path)
        c.check(got == want, f"{rel}", f"expected {want}, got {got}")

    rows = {}
    for scheme in FULL_CHECK:
        rows[scheme] = check_scheme(c, scheme)

    if rows.get("nCoV-2019/V5.3.2") and rows.get("nCoV-2019/V5.4.2"):
        check_delta(c, rows["nCoV-2019/V5.3.2"], rows["nCoV-2019/V5.4.2"])

    check_manifest(c)

    print()
    if c.failures:
        print(f"FAILED: {len(c.failures)} of {c.count} checks")
        return 1
    print(f"PASSED: all {c.count} checks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
