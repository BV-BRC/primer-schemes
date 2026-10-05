## v5.3.2

This is a minimal primer scheme, which uses primer replacement rather than spike-in.

Source: `quick-lab/primerschemes`, `artic-sars-cov-2/400/v5.3.2` (status `validated`,
CC-BY-4.0).

`SARS-CoV-2.primer.bed` and `SARS-CoV-2.reference.fasta` are the upstream files verbatim:

| file | md5 |
| --- | --- |
| `SARS-CoV-2.primer.bed` | `46a3aeddc452678bd0cd33efc1c9558f` |
| `SARS-CoV-2.reference.fasta` | `7f8995394dfc7d5ffeb9fe8322ade58c` |

`SARS-CoV-2.scheme.bed` is columns 1-6 of `primer.bed`, present because
`sars2_assembly/rewrite-primers.pl` requires a `.scheme.bed` in every version directory.

### Resynced to upstream

The bed files were originally imported from the older `quick-lab/SARS-CoV-2` repository
and had drifted from the authoritative `quick-lab/primerschemes` copy. They have been
replaced with the upstream files. Trimming behaviour is unchanged — every coordinate is
identical — but three things differ:

- **`SARS-CoV-2_84_RIGHT_3` is now present** (193 primers, up from 192). The old import
  expressed this oligo as a single entry carrying the IUPAC degenerate base `R`
  (`TGTTCAACACCARTGTCTGTACTC`); upstream ships the two explicit oligos it stands for,
  `..A..` and `..G..`, at identical coordinates. Only the second of those differs from
  the reference, so this scheme has exactly one primer that is not an exact reference
  match: `SARS-CoV-2_84_RIGHT_2` (position 26059, C in the reference, T in the primer).
- **Primer names no longer carry the `_400_` infix** (`SARS-CoV-2_84_RIGHT_2`, not
  `SARS-CoV-2_400_84_RIGHT_2`), matching upstream and the V4/V4.1/V5.4.2 directories.
  Nothing in the pipeline parses the name beyond `LEFT`/`RIGHT`.
- **`primer.bed` now carries the primer sequence** in column 7, and `scheme.bed` is a
  clean 6-column BED. The previous `scheme.bed` had CRLF line endings, a ragged
  19/21/39-column trailing-tab layout, and no final newline.

`SARS-CoV-2.reference.fasta` was already byte-identical to upstream and is unchanged.

Thanks go to the BC Centre for Disease Control and Public Health Laboratory.
