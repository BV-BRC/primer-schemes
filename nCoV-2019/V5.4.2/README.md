## v5.4.2

400 bp ARTIC scheme for SARS-CoV-2. An update to v5.3.2 by the BC Centre for Disease
Control and Public Health Laboratory (BCCDC), adding five spike-in primers that restore
amplification against mutations in the then-dominant JN.1 lineages.

Source: `quick-lab/primerschemes`, `artic-sars-cov-2/400/v5.4.2` (status `validated`,
CC-BY-4.0). Release note:
https://community.artic.network/t/scheme-release-artic-sars-cov2-400-v5-4-2/546

`SARS-CoV-2.primer.bed` and `SARS-CoV-2.reference.fasta` are the upstream files verbatim:

| file | md5 |
| --- | --- |
| `SARS-CoV-2.primer.bed` | `e7897c8be8e488836ce7777e7709ea53` |
| `SARS-CoV-2.reference.fasta` | `d11d06b5d1eb1d85c69e341c3c026e08` |

`SARS-CoV-2.scheme.bed` is columns 1-6 of `primer.bed`, present because
`sars2_assembly/rewrite-primers.pl` requires a `.scheme.bed` in every version directory.

### Changes from v5.3.2

Five primers added (198 total, up from 193). These are spike-ins, not replacements: each
sits at exactly the same coordinates as an existing primer for the same binding site and
differs only in sequence, so the set of intervals trimmed from reads is unchanged.

| chrom | start | end | primer | pool | strand | spike-in |
| --- | --- | --- | --- | --- | --- | --- |
| MN908947.3 | 3560 | 3584 | SARS-CoV-2_11_RIGHT_1 | 1 | - | 1.82 ul |
| MN908947.3 | 21696 | 21722 | SARS-CoV-2_69_RIGHT_1 | 1 | - | 5.08 ul |
| MN908947.3 | 21696 | 21722 | SARS-CoV-2_69_RIGHT_2 | 1 | - | 5.08 ul |
| MN908947.3 | 21927 | 21960 | SARS-CoV-2_70_RIGHT_1 | 2 | - | 3.56 ul |
| MN908947.3 | 22742 | 22774 | SARS-CoV-2_74_LEFT_1 | 2 | + | 1.78 ul |

Spike-in volume is the amount of 100 uM primer stock to add to 200 ul of the v5.3.2
100 uM stock.

Each new primer carries one substitution relative to Wuhan-Hu-1 (MN908947.3), which is
why it binds the drifted lineage. Reference-coordinate positions, 0-based:

| primer | position | Wuhan-Hu-1 | primer |
| --- | --- | --- | --- |
| SARS-CoV-2_11_RIGHT_1 | 3564 | T | C |
| SARS-CoV-2_69_RIGHT_2 | 21710 | C | T |
| SARS-CoV-2_69_RIGHT_1 | 21717 | G | T |
| SARS-CoV-2_70_RIGHT_1 | 21940 | G | T |
| SARS-CoV-2_74_LEFT_1 | 22769 | G | A |

Those five, plus `SARS-CoV-2_84_RIGHT_2` (carried over from v5.3.2), are the only primers
in this scheme that do not match the reference exactly. `tests/validate_schemes.py`
pins that list.

Thanks go to the BC Centre for Disease Control and Public Health Laboratory.
