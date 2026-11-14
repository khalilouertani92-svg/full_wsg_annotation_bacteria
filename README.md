

# genome_function_summary.py

Turns an eggNOG-mapper annotation file into COG, KEGG and GO functional
classifications: Excel tables, a per-protein table, and bar-chart figures
(PNG 600 dpi, PDF, SVG).

## 0. Setup Virtual enviroment (once)
python -m venv wsgvenv

## 1. Install (once)

Python 3.8 or newer, then:

```
pip install matplotlib openpyxl
```

## 2. Download the three reference files

| File | Where to get it | Save as |
|---|---|---|
| NCBI COG definitions | https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab | `cog-24.def.tab` |
| KEGG ko00001 hierarchy | https://rest.kegg.jp/get/br:ko00001 (open, then File > Save as) | `ko00001.keg` |
| Gene Ontology (basic) | https://current.geneontology.org/ontology/go-basic.obo | `go-basic.obo` |

From eggNOG-mapper, download `query.emapper.annotations` (the "Functional
annotations" file, not the Excel version).

## 3. Run

Put the four files and the script in one folder, then:

```
python genome_function_summary.py --annotations query.emapper.annotations --cog-def cog-24.def.tab --kegg ko00001.keg --go go-basic.obo --out results_Pp-BG1
```

## 4. Output (in the `--out` folder)

| File | Content |
|---|---|
| `COG_classification.xlsx` | proteins per COG category, with % |
| `KEGG_classification.xlsx` | proteins per KEGG level-B pathway category, with % |
| `GO_classification.xlsx` | proteins per level-2 GO term, with % |
| `per_protein_classification.tsv` | every protein with its COG, KEGG and GO categories (opens in Excel) |
| `b_COG_classification.png/.pdf/.svg` | panel (b) |
| `c_KEGG_classification.png/.pdf/.svg` | panel (c) |
| `d_GO_classification.png/.pdf/.svg` | panel (d) |

## 5. Check your run (Pp-BG1)

With the Pp-BG1 eggNOG file and the reference versions below, the screen
summary should read:

```
Annotated proteins:                1655
  with a COG category:             1526   (no COG: 129)
  with a KEGG KO:                  1256
  in KEGG pathway categories:      694
  with GO terms:                   677
  with a level-2 GO category:      649
```

Reference versions used for the article figures: NCBI COG 2024; KEGG ko00001
last updated November 13, 2026; GO releases/2026-07-26. Newer KEGG or GO
releases can shift a few counts slightly; report the versions you used.

## Method (for the Methods section)

- **COG:** COG IDs from eggNOG (`COG_category` column) were assigned to
  functional categories with the NCBI COG 2024 definitions (letters are used
  directly if the column already contains them, as in eggNOG-mapper v2).
- **KEGG:** KEGG orthologs were placed in the KEGG BRITE hierarchy ko00001;
  pathway categories are level-A groups 09100-09160.
- **GO:** GO terms were propagated to level 2 (direct `is_a` children of the
  three GO roots) through `is_a` and `part_of` relationships; alternative IDs
  were resolved and obsolete terms replaced where a replacement exists.
- Each protein is counted once per category.
