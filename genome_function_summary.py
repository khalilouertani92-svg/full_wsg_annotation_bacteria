#!/usr/bin/env python3
"""
genome_function_summary.py
==========================

Summarise eggNOG-mapper annotations into COG, KEGG and GO functional
classifications, write per-category count tables (Excel) and a per-protein
table, and draw publication bar charts.

Inputs (all plain-text files):
  1. eggNOG-mapper annotation table       query.emapper.annotations
  2. NCBI COG definitions                 cog-24.def.tab
     https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab
  3. KEGG Orthology hierarchy (ko00001)   ko00001.keg
     https://rest.kegg.jp/get/br:ko00001   (save the page as ko00001.keg)
  4. Gene Ontology (basic version)         go-basic.obo
     https://current.geneontology.org/ontology/go-basic.obo

Method
  COG  : each protein's COG ID (eggNOG v3 / eggNOG 7 "COG_category" column) is
         mapped to its functional category letter(s) with the NCBI COG
         definitions. If the column already holds letters (eggNOG-mapper v2),
         they are used directly. "S" = function unknown; "-" = no COG.
  KEGG : each KEGG ortholog (K number) is placed in the KEGG BRITE ko00001
         hierarchy (level A / level B). The pathway figure uses level-A groups
         09100-09160 (Metabolism ... Human diseases); BRITE (09180) and
         "Not included" (09190) are reported in the table only.
  GO   : each GO term is mapped to its level-2 ancestors (terms whose is_a
         parent is a root: biological_process, cellular_component,
         molecular_function), following is_a and part_of links. Alternative
         IDs are resolved; obsolete terms use their replacement if one exists,
         otherwise they are skipped.
  A protein is counted once per category; a protein in several categories is
  counted in each.

Requirements: Python >= 3.8, matplotlib, openpyxl
  pip install matplotlib openpyxl

Example:
  python genome_function_summary.py \
      --annotations query.emapper.annotations \
      --cog-def cog-24.def.tab --kegg ko00001.keg --go go-basic.obo \
      --out results_Pp-BG1
"""
import argparse
import os
import re
import sys
from collections import Counter, OrderedDict, defaultdict

# --------------------------------------------------------------------------- #
# Reference data that does not change: COG functional categories
# --------------------------------------------------------------------------- #
COG_GROUPS = [
    ('Information storage and processing', 'JAKLB'),
    ('Cellular processes and signaling', 'DYVTMNZWUOX'),
    ('Metabolism', 'CGEFHIPQ'),
    ('Poorly characterized', 'RS'),
]
COG_NAMES = {
    'J': 'Translation, ribosomal structure and biogenesis', 'A': 'RNA processing and modification',
    'K': 'Transcription', 'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning', 'Y': 'Nuclear structure',
    'V': 'Defense mechanisms', 'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis', 'N': 'Cell motility', 'Z': 'Cytoskeleton',
    'W': 'Extracellular structures', 'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'X': 'Mobilome: prophages, transposons', 'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism', 'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism', 'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism', 'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'R': 'General function prediction only', 'S': 'Function unknown',
}
GO_ROOTS = OrderedDict([('GO:0008150', 'Biological process'),
                        ('GO:0005575', 'Cellular component'),
                        ('GO:0003674', 'Molecular function')])
KEGG_PATHWAY_GROUPS = ('09100', '09120', '09130', '09140', '09150', '09160')


# --------------------------------------------------------------------------- #
# Readers
# --------------------------------------------------------------------------- #
def read_emapper(path):
    """Return list of dicts, one per annotated protein."""
    header, rows = None, []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#query'):
                header = line[1:].split('\t')
                continue
            if not line or line.startswith('#'):
                continue
            if header is None:
                sys.exit('ERROR: no "#query" header line found in ' + path)
            rows.append(dict(zip(header, line.split('\t'))))
    for col in ('query', 'COG_category', 'KEGG_ko', 'GOs'):
        if rows and col not in rows[0]:
            sys.exit(f'ERROR: column "{col}" missing from {path}')
    return rows


def read_cog_def(path):
    """COG ID -> (letters, name)."""
    d = {}
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3 and re.match(r'^COG\d{4}$', f[0]):
                d[f[0]] = (f[1], f[2])
    return d


def read_kegg_ko00001(path):
    """Return (ko -> set of B codes, B code -> (A code, A name, B name), ordered B list, release)."""
    ko2b = defaultdict(set)
    binfo = OrderedDict()
    a_code = a_name = b_code = None
    release = ''
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#Last updated'):
                release = line.split(':', 1)[1].strip()
            if not line or line[0] not in 'ABD':
                continue
            if line[0] == 'A':
                m = re.match(r'^A\s*(\d{5})\s+(.*)$', line)
                if m:
                    a_code, a_name = m.group(1), m.group(2).strip()
            elif line[0] == 'B':
                m = re.match(r'^B\s+(\d{5})\s+(.*)$', line)
                if m:
                    b_code = m.group(1)
                    binfo[b_code] = (a_code, a_name, m.group(2).strip())
            else:
                m = re.match(r'^D\s+(K\d{5})', line)
                if m and b_code:
                    ko2b[m.group(1)].add(b_code)
    if not binfo:
        sys.exit('ERROR: could not parse KEGG hierarchy from ' + path)
    return ko2b, binfo, release


def read_go_obo(path):
    """Return term dict, alt_id map and data-version."""
    terms, alt, version = {}, {}, ''
    cur = None
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('data-version:'):
                version = line.split(':', 1)[1].strip()
            if line == '[Term]':
                cur = {'is_a': [], 'part_of': [], 'alt': [], 'obsolete': False, 'replaced_by': None}
                continue
            if line.startswith('['):
                cur = None
                continue
            if cur is None or ': ' not in line:
                continue
            key, val = line.split(': ', 1)
            if key == 'id':
                cur['id'] = val.strip()
                terms[cur['id']] = cur
            elif key == 'name':
                cur['name'] = val.strip()
            elif key == 'namespace':
                cur['namespace'] = val.strip()
            elif key == 'is_a':
                cur['is_a'].append(val.split()[0])
            elif key == 'relationship' and val.startswith('part_of '):
                cur['part_of'].append(val.split()[1])
            elif key == 'alt_id':
                cur['alt'].append(val.strip())
            elif key == 'is_obsolete' and val.strip() == 'true':
                cur['obsolete'] = True
            elif key == 'replaced_by':
                cur['replaced_by'] = val.split()[0]
    for t in terms.values():
        for a in t['alt']:
            alt[a] = t['id']
    return terms, alt, version


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #
def classify(rows, cogdef, ko2b, binfo, terms, alt):
    level2 = sorted(t['id'] for t in terms.values()
                    if not t['obsolete'] and any(p in GO_ROOTS for p in t['is_a']))
    level2_set = set(level2)
    memo = {}

    def ancestors(tid):
        if tid in memo:
            return memo[tid]
        memo[tid] = {tid}  # guards against cycles
        res = {tid}
        t = terms.get(tid)
        if t:
            for p in t['is_a'] + t['part_of']:
                res |= ancestors(p)
        memo[tid] = res
        return res

    def resolve(goid):
        if goid not in terms and goid in alt:
            goid = alt[goid]
        t = terms.get(goid)
        if t is None:
            return None
        if t['obsolete']:
            return t['replaced_by'] if t['replaced_by'] in terms else None
        return goid

    out, stats = [], Counter()
    cog_unmapped, go_skipped = set(), set()
    for r in rows:
        rec = {'query': r['query'], 'COG_ID': '', 'COG_letters': '', 'COG_description': '',
               'KEGG_KO': '', 'KEGG_categories': [], 'GO_terms': '', 'GO_level2': []}
        # --- COG
        c = r['COG_category'].strip()
        if re.match(r'^COG\d{4}$', c):
            rec['COG_ID'] = c
            if c in cogdef:
                rec['COG_letters'], rec['COG_description'] = cogdef[c]
            else:
                cog_unmapped.add(c)
        elif c and c != '-' and all(ch in COG_NAMES for ch in c):
            rec['COG_letters'] = c
            rec['COG_description'] = '; '.join(COG_NAMES[ch] for ch in c)
        # --- KEGG
        k = r['KEGG_ko'].strip()
        if k and k != '-':
            kos = [x.replace('ko:', '').strip() for x in k.split(',') if x.strip() not in ('', '-')]
            rec['KEGG_KO'] = ','.join(kos)
            cats = set()
            for ko in kos:
                cats |= ko2b.get(ko, set())
            rec['KEGG_categories'] = sorted(cats)
            stats['kegg_with_ko'] += 1
            if any(binfo[b][0] in KEGG_PATHWAY_GROUPS for b in cats):
                stats['kegg_in_pathway_groups'] += 1
        # --- GO
        g = r['GOs'].strip()
        if g and g != '-':
            ids = [x.strip() for x in g.split(',') if x.strip()]
            rec['GO_terms'] = ','.join(ids)
            l2 = set()
            for gid in ids:
                rid = resolve(gid)
                if rid is None:
                    go_skipped.add(gid)
                    continue
                l2 |= ancestors(rid) & level2_set
            rec['GO_level2'] = sorted(l2)
            stats['go_with_terms'] += 1
            if l2:
                stats['go_with_level2'] += 1
        if rec['COG_letters']:
            stats['cog_assigned'] += 1
        elif c == '-' or not c:
            stats['cog_none'] += 1
        out.append(rec)
    return out, level2, stats, cog_unmapped, go_skipped


def count_categories(recs, binfo, terms):
    cog = Counter()
    for r in recs:
        for ch in set(r['COG_letters']):
            cog[ch] += 1
    kegg = Counter(b for r in recs for b in r['KEGG_categories'])
    go = Counter(t for r in recs for t in r['GO_level2'])
    return cog, kegg, go


# --------------------------------------------------------------------------- #
# Excel output
# --------------------------------------------------------------------------- #
def _excel_table(path, sheet, headers, rows, total_label, total_value, pct_header, notes, widths):
    """rows: list of lists; the last column before % is the count. % is a live formula."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill, Alignment
    wb = Workbook()
    ws = wb.active
    ws.title = sheet
    F = 'Arial'
    hdr_font = Font(name=F, bold=True, color='FFFFFF')
    fill = PatternFill('solid', fgColor='1F4E78')
    ws.append(headers + [pct_header])
    for c in ws[1]:
        c.font, c.fill = hdr_font, fill
        c.alignment = Alignment(vertical='center', wrap_text=True)
    ncol = len(headers)
    cnt_col = chr(ord('A') + ncol - 1)
    pct_col = chr(ord('A') + ncol)
    first, last = 2, len(rows) + 1
    tot_row = last + 2
    for i, row in enumerate(rows, start=first):
        ws.append(row)
        ws[f'{pct_col}{i}'] = f'=IF(${cnt_col}${tot_row}=0,0,{cnt_col}{i}/${cnt_col}${tot_row})'
        ws[f'{pct_col}{i}'].number_format = '0.0%'
    ws.cell(tot_row, ncol - 1, total_label).font = Font(name=F, bold=True)
    ws.cell(tot_row, ncol, total_value).font = Font(name=F, bold=True)
    r = tot_row + 2
    for n in notes:
        ws.cell(r, 1, n).font = Font(name=F, italic=True, color='52514E')
        r += 1
    for row in ws.iter_rows(min_row=2, max_row=last):
        for c in row:
            c.font = Font(name=F)
    for i, w in enumerate(widths + [14]):
        ws.column_dimensions[chr(ord('A') + i)].width = w
    ws.freeze_panes = 'A2'
    wb.save(path)


def write_cog_excel(path, cog, n_assigned, n_none, notes):
    rows = []
    for gname, letters in COG_GROUPS:
        for L in letters:
            rows.append([gname, L, COG_NAMES[L], cog.get(L, 0)])
    notes = [f'Proteins without a COG assignment: {n_none}.',
             'Percentages are relative to proteins with a COG assignment (total above).',
             'Proteins with multi-letter COG categories are counted under each letter, '
             'so the column sum exceeds the total.'] + notes
    _excel_table(path, 'COG classification', ['COG group', 'Letter', 'Functional category', 'Proteins'],
                 rows, 'Proteins with a COG assignment', n_assigned, '% of COG-assigned proteins',
                 notes, [36, 8, 58, 12])


def write_kegg_excel(path, kegg, binfo, n_pathway, n_ko, notes):
    rows = []
    for b, (a_code, a_name, b_name) in binfo.items():
        if a_code in KEGG_PATHWAY_GROUPS:
            rows.append([f'{a_code} {a_name}', f'{b} {b_name}', kegg.get(b, 0)])
    extra = [f'{b} {binfo[b][2]} ({binfo[b][1]}): {kegg.get(b, 0)}'
             for b in binfo if binfo[b][0] not in KEGG_PATHWAY_GROUPS]
    notes = [f'Proteins with at least one KEGG ortholog (KO): {n_ko}.',
             'Percentages are relative to proteins assigned to at least one pathway category (total above).',
             'A protein is counted once per category and may appear in several categories.',
             'Not shown in the table (BRITE / unclassified): ' + '; '.join(extra)] + notes
    _excel_table(path, 'KEGG classification', ['KEGG level A', 'KEGG level B', 'Proteins'],
                 rows, 'Proteins in pathway categories', n_pathway, '% of pathway-assigned proteins',
                 notes, [44, 52, 12])


def write_go_excel(path, go, level2, terms, n_l2, n_go, notes):
    rows = []
    for root, dom in GO_ROOTS.items():
        sub = [t for t in level2 if terms[t].get('namespace', '').replace('_', ' ') == dom.lower()]
        sub.sort(key=lambda t: (-go.get(t, 0), t))
        for t in sub:
            rows.append([dom, t, terms[t]['name'], go.get(t, 0)])
    notes = [f'Proteins with at least one GO term: {n_go}.',
             'Percentages are relative to proteins with at least one level-2 GO category (total above).',
             'A protein is counted once per category and may appear in several categories.'] + notes
    _excel_table(path, 'GO level 2', ['GO domain', 'GO ID', 'GO term (level 2)', 'Proteins'],
                 rows, 'Proteins with a level-2 GO category', n_l2, '% of GO-assigned proteins',
                 notes, [22, 13, 62, 12])


def write_per_protein(path, recs, binfo, terms):
    header = ['query', 'COG_ID', 'COG_letters', 'COG_description', 'KEGG_KO', 'KEGG_categories',
              'GO_terms', 'GO_level2']
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(header) + '\n')
        for r in recs:
            kc = '; '.join(f'{b} {binfo[b][2]}' for b in r['KEGG_categories'])
            gl = '; '.join(f'{t} {terms[t]["name"]}' for t in r['GO_level2'])
            fh.write('\t'.join([r['query'], r['COG_ID'], r['COG_letters'], r['COG_description'],
                                r['KEGG_KO'], kc, r['GO_terms'], gl]) + '\n')


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300']
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def hbar(rows, title, fname, legend_title, pct_total=None, width=7.2, legend_loc='lower right'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.patches import Patch
    fonts = {f.name for f in font_manager.fontManager.ttflist}
    font = next((f for f in ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'] if f in fonts), 'DejaVu Sans')
    plt.rcParams.update({'font.family': font, 'font.size': 8, 'axes.linewidth': 0.6,
                         'xtick.major.width': 0.6, 'ytick.major.width': 0,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    n = len(rows)
    fig, ax = plt.subplots(figsize=(width, 0.9 + n * 0.19))
    ys = list(range(n))[::-1]
    vmax = max(r[1] for r in rows)
    for y, (lab, v, gi, gn) in zip(ys, rows):
        ax.barh(y, v, height=0.72, color=SLOTS[gi], edgecolor='white', linewidth=0.5, zorder=3)
        txt = f'{v}' if pct_total is None else f'{v} ({100 * v / pct_total:.1f}%)'
        ax.text(v + vmax * 0.012, y, txt, va='center', ha='left', fontsize=7, color=INK2, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5, color=INK)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlim(0, vmax * 1.2)
    ax.set_xlabel('Number of proteins', fontsize=8, color=INK)
    ax.xaxis.grid(True, color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color('#9a9993')
    ax.tick_params(axis='x', colors=INK2, labelsize=7)
    groups = []
    for r in rows:
        if (r[2], r[3]) not in groups:
            groups.append((r[2], r[3]))
    ax.legend(handles=[Patch(facecolor=SLOTS[gi], edgecolor='none', label=gn) for gi, gn in groups],
              title=legend_title, loc=legend_loc, frameon=False, fontsize=7, title_fontsize=7.5,
              handlelength=1.1, handleheight=0.9, borderaxespad=0.2)
    ax.set_title(title, loc='left', fontsize=9.5, fontweight='bold', color=INK, pad=8)
    fig.tight_layout()
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(fname + '.' + ext, dpi=600 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig)


def draw_figures(outdir, cog, n_assigned, kegg, binfo, go, level2, terms):
    rows = []
    for gi, (gname, letters) in enumerate(COG_GROUPS):
        for L in letters:
            if cog.get(L, 0) > 0:
                rows.append((f'[{L}] {COG_NAMES[L]}', cog[L], gi, gname))
    hbar(rows, '(b) COG functional classification', os.path.join(outdir, 'b_COG_classification'),
         'COG group', pct_total=n_assigned, width=7.4, legend_loc='center right')

    rows, groups = [], []
    for b, (a_code, a_name, b_name) in binfo.items():
        if a_code in KEGG_PATHWAY_GROUPS and kegg.get(b, 0) > 0:
            if a_name not in groups:
                groups.append(a_name)
            rows.append((b_name, kegg[b], groups.index(a_name), a_name[0].upper() + a_name[1:].lower()))
    hbar(rows, '(c) KEGG pathway classification', os.path.join(outdir, 'c_KEGG_classification'),
         'KEGG category', width=7.0)

    rows = []
    for gi, (root, dom) in enumerate(GO_ROOTS.items()):
        sub = [t for t in level2 if terms[t].get('namespace', '').replace('_', ' ') == dom.lower()
               and go.get(t, 0) > 0]
        sub.sort(key=lambda t: (-go[t], t))
        for t in sub:
            name = terms[t]['name'].replace('biological process involved in ', '')
            rows.append((name, go[t], gi, dom))
    hbar(rows, '(d) GO classification (level 2)', os.path.join(outdir, 'd_GO_classification'),
         'GO domain', width=6.8)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--annotations', required=True, help='eggNOG-mapper *.emapper.annotations file')
    ap.add_argument('--cog-def', required=True, help='NCBI cog-24.def.tab')
    ap.add_argument('--kegg', required=True, help='KEGG ko00001 hierarchy (rest.kegg.jp/get/br:ko00001)')
    ap.add_argument('--go', required=True, help='go-basic.obo')
    ap.add_argument('--out', default='functional_summary', help='output folder')
    ap.add_argument('--no-figures', action='store_true', help='skip drawing figures')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows = read_emapper(a.annotations)
    cogdef = read_cog_def(a.cog_def)
    ko2b, binfo, kegg_release = read_kegg_ko00001(a.kegg)
    terms, alt, go_version = read_go_obo(a.go)

    recs, level2, st, cog_unmapped, go_skipped = classify(rows, cogdef, ko2b, binfo, terms, alt)
    cog, kegg, go = count_categories(recs, binfo, terms)

    src = [f'Source annotation file: {os.path.basename(a.annotations)} ({len(rows)} annotated proteins).',
           f'References: NCBI COG 2024 ({os.path.basename(a.cog_def)}); KEGG ko00001 '
           f'(last updated {kegg_release or "n/a"}); GO {go_version or os.path.basename(a.go)}.',
           'Generated with genome_function_summary.py.']
    write_cog_excel(os.path.join(a.out, 'COG_classification.xlsx'), cog, st['cog_assigned'], st['cog_none'], src)
    write_kegg_excel(os.path.join(a.out, 'KEGG_classification.xlsx'), kegg, binfo,
                     st['kegg_in_pathway_groups'], st['kegg_with_ko'], src)
    write_go_excel(os.path.join(a.out, 'GO_classification.xlsx'), go, level2, terms,
                   st['go_with_level2'], st['go_with_terms'], src)
    write_per_protein(os.path.join(a.out, 'per_protein_classification.tsv'), recs, binfo, terms)
    if not a.no_figures:
        draw_figures(a.out, cog, st['cog_assigned'], kegg, binfo, go, level2, terms)

    print(f'Annotated proteins:                {len(rows)}')
    print(f'  with a COG category:             {st["cog_assigned"]}   (no COG: {st["cog_none"]})')
    print(f'  with a KEGG KO:                  {st["kegg_with_ko"]}')
    print(f'  in KEGG pathway categories:      {st["kegg_in_pathway_groups"]}')
    print(f'  with GO terms:                   {st["go_with_terms"]}')
    print(f'  with a level-2 GO category:      {st["go_with_level2"]}')
    if cog_unmapped:
        print(f'WARNING: {len(cog_unmapped)} COG IDs not in {a.cog_def}: {sorted(cog_unmapped)[:10]}')
    if go_skipped:
        print(f'Note: {len(go_skipped)} obsolete/unknown GO terms skipped: {sorted(go_skipped)}')
    print(f'Results written to: {os.path.abspath(a.out)}')


if __name__ == '__main__':
    main()
