/* C10+ cis-pQTL evidence-tier portal — vanilla JS, no dependencies, no network calls except the local data file. MIT License. */
(function () {
  'use strict';
  var DATA_URL = 'data/c10plus_tiers_v2.2.json';
  var D = null, E = null, S = null, N = 0;       // data columns, enums, studies, n rows
  var studyKeys = [];                             // per study: [SYMBOL, ENSG, UNIPROT, OID]
  var symList = [];                               // unique symbols (upper) -> study indices
  var symToStudies = {};
  var state = { q: '', tiers: [1, 2, 3, 4], onlyflag: false, beta: 'any', fuzzy: false, sort: null, dir: 1, page: 0, size: 100 };
  var current = [];                               // [{i, r, m}] filtered rows (r = rank, m = match type)
  var $ = function (id) { return document.getElementById(id); };

  var TIER_SHORT = { 1: 'anchored, no PAV', 2: 'anchored + PAV / isoform-only', 3: 'unresolved: no anchor, no PAV', 4: 'caution: no anchor + PAV' };
  var MS_LABEL = { supported: 'supported', refuted: 'refuted', indeterminate_underpowered: 'indeterminate', indeterminate: 'indeterminate', not_assessed: '–' };
  var PEP_LABEL = { supported: 'peptide: supported', refuted: 'peptide: refuted', 'intermediate_0.2_0.8': 'peptide: intermediate', not_assessed: '' };
  var DISC_LABEL = { overlap: 'discovery: overlap', no_overlap: 'discovery: no overlap', not_assessed: '' };
  var MATCH_NAME = ['exact', 'prefix', 'contains', 'fuzzy'];

  function esc(s) { return String(s).replace(/[&<>"']/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); }
  function v(f, i) { return D[f][i]; }
  function ev(f, i) { return E[f][D[f][i]]; }                  // enum value
  function st(i) { return S[D.s[i]]; }                          // study row
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function raw(x) { return x === null || x === undefined ? '' : String(x); }
  function span(f, val, text, cls) { return '<span data-f="' + f + '" data-v="' + esc(raw(val)) + '"' + (cls ? ' class="' + cls + '"' : '') + '>' + esc(text) + '</span>'; }
  function fmtNum(x, kind) {
    if (x === null || x === undefined) return 'n/a';
    if (kind === 'p') return x === 0 ? '0' : (Math.abs(x) < 1e-3 ? x.toExponential(2) : x.toPrecision(3));
    if (kind === 'beta') return x.toPrecision(4);
    if (kind === 'mlp') return x.toFixed(1);
    if (kind === 'h4') return x.toFixed(2);
    return String(x);
  }
  function repCell(f, i, dirf) {
    var r = v(f, i), t = r === null ? 'n/a' : (r === 1 ? 'yes' : 'no');
    var h = span(f, r, t, r === null ? 'na' : (r === 1 ? 'yes' : 'no'));
    if (dirf) {
      var d = ev(dirf, i);
      if (d !== 'NA') h += '<span class="sub">' + span(dirf, d, d, d === 'opposite' ? 'opp' : '') + '</span>';
    }
    return h;
  }
  function yn(f, i, yesText, noText) { var x = v(f, i); return span(f, x, x === 1 ? yesText : noText, x === 1 ? 'yes' : 'na'); }

  /* ---------- column definitions (main table) ---------- */
  var COLS = [
    { k: 'tier', h: 'Tier', tip: 'Pre-registered rule tier. T1 eQTL-anchored, no PAV; T2 anchored with PAV, or isoform-only (sQTL/tuQTL); T3 no anchor, no PAV; T4 no anchor with PAV (regardless of platform concordance). T4 is not a false-positive call.',
      sort: function (i) { return v('tier', i); },
      cell: function (i) { var t = v('tier', i); return '<span class="badge t' + t + '">' + span('tier', t, 'T' + t) + '</span><span class="sub">' + span('tier_label', D._tl[t], TIER_SHORT[t]) + '</span>'; } },
    { k: 'symbol', h: 'Gene', tip: 'HGNC symbol of the cognate gene (Olink target) and its Ensembl gene ID. The chip shows how the row matched your query.',
      sort: function (i) { return st(i)[1]; },
      cell: function (i, m) { var s = st(i); return span('symbol', s[1], s[1]) + (m !== undefined && m !== null ? ' <span class="match">' + MATCH_NAME[m] + '</span>' : '') + '<span class="sub mono">' + span('gene_id', s[2], s[2]) + '</span>'; } },
    { k: 'uniprot', h: 'UniProt / OID', tip: 'UniProt accession and Olink assay OID of the UKB-PPP Olink Explore assay (from the Open Targets studyId).',
      sort: function (i) { return st(i)[3]; },
      cell: function (i) { var s = st(i); return span('uniprot', s[3], s[3], 'mono') + '<span class="sub mono">' + span('olink_oid', s[4], s[4]) + '</span>'; } },
    { k: 'lead_variant_id', h: 'Lead variant', tip: 'Lead variant of the 95% credible set, chrom_pos_ref_alt on GRCh38 (Open Targets identifier); below, the number of variants in the 95% credible set.',
      sort: function (i) { var p = v('lead_variant_id', i).split('_'); var c = p[0] === 'X' ? 23 : +p[0]; return c * 1e9 + (+p[1]); },
      cell: function (i) { return span('lead_variant_id', v('lead_variant_id', i), v('lead_variant_id', i), 'mono') + '<span class="sub">CS95 n=' + span('n_variants_cs95', v('n_variants_cs95', i), v('n_variants_cs95', i)) + '</span>'; } },
    { k: 'lead_beta', h: 'β', num: true, tip: 'Open Targets lead β for the alt allele (Olink NPX scale), 4 significant digits.',
      sort: function (i) { return v('lead_beta', i); },
      cell: function (i) { var x = v('lead_beta', i); return span('lead_beta', x, fmtNum(x, 'beta')); } },
    { k: 'lead_mlog10p', h: '−log10 P', num: true, tip: '−log10 P of the lead variant in UKB-PPP (capped near 322 for the strongest signals).',
      sort: function (i) { return v('lead_mlog10p', i); },
      cell: function (i) { var x = v('lead_mlog10p', i); return span('lead_mlog10p', x, fmtNum(x, 'mlp')); } },
    { k: 'pav', h: 'PAV', tip: 'Primary PAV flag: the 95% credible set contains a protein-altering variant of the cognate gene (VEP 116, any transcript). “lead” = the lead variant itself is a cognate PAV; PIP = summed PIP of cognate PAVs.',
      sort: function (i) { return v('pav_cognate_cs95_vep116_alltx', i) * 10 + v('pav_lead_variant', i); },
      cell: function (i) {
        var p = v('pav_cognate_cs95_vep116_alltx', i);
        var h = span('pav_cognate_cs95_vep116_alltx', p, p === 1 ? 'yes' : 'no', p === 1 ? 'yes' : 'na');
        if (p === 1) h += '<span class="sub">' + (v('pav_lead_variant', i) === 1 ? span('pav_lead_variant', 1, 'lead') + ' · ' : '') + 'PIP ' + span('pav_pip_sum', v('pav_pip_sum', i), v('pav_pip_sum', i).toFixed(2)) + '</span>';
        return h; } },
    { k: 'anchor', h: 'eQTL anchor', tip: 'Cognate eQTL anchor: GTEx v10 or eQTL Catalogue eQTL colocalisation PP.H4 ≥ 0.8; below, the number of eQTL datasets with PP.H4 ≥ 0.8 and whether an sQTL/tuQTL (isoform) anchor exists.',
      sort: function (i) { return v('anchor_eqtl', i) * 1000 + v('n_anchor_eqtl_datasets', i); },
      cell: function (i) {
        var a = v('anchor_eqtl', i);
        var h = span('anchor_eqtl', a, a === 1 ? 'yes' : 'no', a === 1 ? 'yes' : 'na');
        var sub = [];
        if (a === 1) sub.push('n=' + span('n_anchor_eqtl_datasets', v('n_anchor_eqtl_datasets', i), v('n_anchor_eqtl_datasets', i)));
        if (v('anchor_isoform', i) === 1) sub.push(span('anchor_isoform', 1, 'isoform'));
        return h + (sub.length ? '<span class="sub">' + sub.join(' · ') + '</span>' : ''); } },
    { k: 'rep_sun2018', h: 'Sun2018', tip: 'SomaScan v3 (Sun et al. 2018): colocalisation PP.H4 ≥ 0.8 with a Sun 2018 cis credible set. n/a = gene not assayed. Below: effect direction at the top shared variant.',
      sort: function (i) { return v('rep_sun2018', i); }, cell: function (i) { return repCell('rep_sun2018', i, 'dir_sun2018'); } },
    { k: 'rep_decode_v4', h: 'deCODE', tip: 'SomaScan v4 deCODE (Eldjarn et al. 2023): a deCODE cis signal (lead or LD-class member) lies in the UKB 95% credible set. n/a = gene not in the deCODE tables. Below: direction at the deCODE lead.',
      sort: function (i) { return v('rep_decode_v4', i); }, cell: function (i) { return repCell('rep_decode_v4', i, 'dir_decode_v4'); } },
    { k: 'rep_aric_ea', h: 'ARIC-EA', tip: 'SomaScan v4 ARIC European ancestry (Zhang et al. 2022): any ARIC signal variant (lead, conditional or fine-mapped) lies in the UKB 95% set. n/a = no ARIC cis signal for the gene. Below: direction.',
      sort: function (i) { return v('rep_aric_ea', i); }, cell: function (i) { return repCell('rep_aric_ea', i, 'dir_aric_ea'); } },
    { k: 'rep_aric_aa', h: 'ARIC-AA', tip: 'SomaScan v4 ARIC African ancestry (Zhang et al. 2022), same rule as ARIC-EA. n/a = no ARIC-AA cis signal for the gene.',
      sort: function (i) { return v('rep_aric_aa', i); }, cell: function (i) { return repCell('rep_aric_aa', i, null); } },
    { k: 'rep_fenland', h: 'Fenland', tip: 'SomaScan v4 Fenland (Pietzner et al. 2021; GWAS Catalog curated associations, name-mapped): the Fenland cis lead lies in the UKB 95% set. n/a = no Fenland cis signal. Below: direction. Note: deCODE, ARIC and Fenland are the same platform.',
      sort: function (i) { return v('rep_fenland', i); }, cell: function (i) { return repCell('rep_fenland', i, 'dir_fenland'); } },
    { k: 'ms_primary', h: 'MS (MSPA)', tip: 'Pre-registered MS label from the MS pQTL assessment (Suhre et al. 2025 tables): supported = Bonferroni MS replication; refuted = power > 99% but not replicated; indeterminate = underpowered; – = not assessed. Below: MSPA peptide label. Not ground truth.',
      sort: function (i) { return ['supported', 'indeterminate_underpowered', 'refuted', 'not_assessed'].indexOf(ev('ms_primary', i)); },
      cell: function (i) {
        var m = ev('ms_primary', i), p = ev('ms_mspa_peptide', i);
        return span('ms_primary', m, MS_LABEL[m], m === 'supported' ? 'yes' : (m === 'refuted' ? 'no' : 'na')) + (PEP_LABEL[p] ? '<span class="sub">' + span('ms_mspa_peptide', p, PEP_LABEL[p]) + '</span>' : ''); } },
    { k: 'ms_ext_primary', h: 'MS (Seer full)', tip: 'Addendum A1 label from Seer Proteograph full summary statistics: supported = index variant P < 0.05/1,779 with the same direction; refuted = power ≥ 0.99 and not supported. Same MS study as MSPA (κ = 0.89) — not independent. Below: overlap with published Seer discovery pQTLs.',
      sort: function (i) { return ['supported', 'indeterminate', 'refuted', 'not_assessed'].indexOf(ev('ms_ext_primary', i)); },
      cell: function (i) {
        var m = ev('ms_ext_primary', i), d = ev('ms_seer_discovery', i);
        return span('ms_ext_primary', m, MS_LABEL[m], m === 'supported' ? 'yes' : (m === 'refuted' ? 'no' : 'na')) + (DISC_LABEL[d] ? '<span class="sub">' + span('ms_seer_discovery', d, DISC_LABEL[d]) + '</span>' : ''); } },
    { k: 'xa_replicated', h: 'Cross-ancestry', tip: 'Olink UKB-PPP non-European ancestries (fixed-effect meta): replicated = β aligned to the EUR direction > 0 and P < 0.05. n/a = not evaluable. Below: number of contributing ancestries and meta P. Same platform, so epitope effects are not excluded.',
      sort: function (i) { return v('xa_replicated', i); },
      cell: function (i) {
        var h = repCell('xa_replicated', i, null), k = v('xa_n_ancestries', i);
        if (k !== null) h += '<span class="sub">k=' + span('xa_n_ancestries', k, k) + ' · P ' + span('xa_meta_p', v('xa_meta_p', i), fmtNum(v('xa_meta_p', i), 'p')) + '</span>';
        return h; } },
    { k: 'flag_concordant_unanchored_PAV', h: 'Concordant + unanchored + PAV', tip: 'Key caution class: replicated on at least one SomaScan layer (platform-concordant) but without an eQTL anchor and with a cognate PAV. Agreement may reflect a shared epitope effect; confirmation needed.',
      sort: function (i) { return v('flag_concordant_unanchored_PAV', i); },
      cell: function (i) { var x = v('flag_concordant_unanchored_PAV', i); return span('flag_concordant_unanchored_PAV', x, x === 1 ? 'yes' : '–', x === 1 ? 'flag' : 'na'); } },
    { k: 'flag_beta_inconsistent', h: 'β-inconsistent', tip: 'Post hoc flag: |β| is more than 5× the value expected from P, MAF and N = 34,000 (Open Targets β inconsistent with P/MAF). Treat the β of these rows with caution.',
      sort: function (i) { return v('flag_beta_inconsistent', i); },
      cell: function (i) { var x = v('flag_beta_inconsistent', i); return span('flag_beta_inconsistent', x, x === 1 ? 'yes' : '–', x === 1 ? 'flag' : 'na'); } }
  ];

  /* ---------- CSV / detail fields (resource field names) ---------- */
  var EXPORT = ['studyLocusId', 'studyId', 'symbol', 'gene_id', 'uniprot', 'olink_oid', 'chrom', 'lead_variant_id', 'lead_beta', 'lead_mlog10p',
    'lead_maf_gnomad_nfe', 'n_variants_cs95', 'tier', 'tier_label', 'anchor_eqtl', 'n_anchor_eqtl_datasets', 'h4max_eqtl_gtex_v10', 'h4max_eqtl_eqtlcatalogue',
    'anchor_isoform', 'h4max_sqtl_tuqtl', 'pav_cognate_cs95_vep116_alltx', 'pav_cognate_cs95_canonical', 'pav_lead_variant', 'pav_pip_sum',
    'rep_sun2018', 'dir_sun2018', 'rep_decode_v4', 'dir_decode_v4', 'rep_aric_ea', 'dir_aric_ea', 'rep_aric_aa', 'rep_fenland', 'dir_fenland',
    'platform_concordance', 'ms_primary', 'ms_mspa_peptide', 'ms_seer_discovery', 'ms_ext_primary', 'ms_ext_nominal', 'seer_ext_p', 'seer_ext_opposite_sig',
    'xa_n_ancestries', 'xa_meta_beta_aligned', 'xa_meta_p', 'xa_replicated', 'flag_concordant_unanchored_PAV', 'lead_beta_to_expected_ratio',
    'flag_beta_inconsistent', 'n_other_genes_trans_cs_shared'];
  var GROUPS = [
    ['Identifiers', ['studyLocusId', 'studyId', 'symbol', 'gene_id', 'uniprot', 'olink_oid', 'chrom', 'lead_variant_id']],
    ['Lead signal', ['lead_beta', 'lead_mlog10p', 'lead_maf_gnomad_nfe', 'n_variants_cs95', 'lead_beta_to_expected_ratio', 'flag_beta_inconsistent', 'n_other_genes_trans_cs_shared']],
    ['Tier and transcript anchoring', ['tier', 'tier_label', 'anchor_eqtl', 'n_anchor_eqtl_datasets', 'h4max_eqtl_gtex_v10', 'h4max_eqtl_eqtlcatalogue', 'anchor_isoform', 'h4max_sqtl_tuqtl']],
    ['Protein-altering variants', ['pav_cognate_cs95_vep116_alltx', 'pav_cognate_cs95_canonical', 'pav_lead_variant', 'pav_pip_sum']],
    ['Affinity replication (SomaScan)', ['rep_sun2018', 'dir_sun2018', 'rep_decode_v4', 'dir_decode_v4', 'rep_aric_ea', 'dir_aric_ea', 'rep_aric_aa', 'rep_fenland', 'dir_fenland', 'platform_concordance', 'flag_concordant_unanchored_PAV']],
    ['Mass spectrometry (one Seer study; two label sets)', ['ms_primary', 'ms_mspa_peptide', 'ms_seer_discovery', 'ms_ext_primary', 'ms_ext_nominal', 'seer_ext_p', 'seer_ext_opposite_sig']],
    ['Cross-ancestry (UKB-PPP non-EUR)', ['xa_n_ancestries', 'xa_meta_beta_aligned', 'xa_meta_p', 'xa_replicated']]
  ];
  var STUDY_F = { studyId: 0, symbol: 1, gene_id: 2, uniprot: 3, olink_oid: 4, chrom: 5 };
  function rawField(f, i) {
    if (f in STUDY_F) return st(i)[STUDY_F[f]];
    if (f === 'tier_label') return D._tl[v('tier', i)];
    if (E[f]) return ev(f, i);
    var x = D[f][i];
    return x === null || x === undefined ? '' : x;
  }

  /* ---------- search ---------- */
  function osa(a, b, max) {                       // optimal string alignment distance with early exit
    var la = a.length, lb = b.length;
    if (Math.abs(la - lb) > max) return max + 1;
    var d = [], i, j;
    for (i = 0; i <= la; i++) { d[i] = [i]; }
    for (j = 0; j <= lb; j++) d[0][j] = j;
    for (i = 1; i <= la; i++) {
      var rowMin = Infinity;
      for (j = 1; j <= lb; j++) {
        var c = a[i - 1] === b[j - 1] ? 0 : 1;
        var x = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + c);
        if (i > 1 && j > 1 && a[i - 1] === b[j - 2] && a[i - 2] === b[j - 1]) x = Math.min(x, d[i - 2][j - 2] + 1);
        d[i][j] = x; if (x < rowMin) rowMin = x;
      }
      if (rowMin > max) return max + 1;
    }
    return d[la][lb];
  }
  function normTerm(t) {
    t = t.trim().toUpperCase();
    var m = t.match(/^(?:CHR)?([0-9]{1,2}|X|Y|MT?)[:_\-]([0-9]+)(?:[:_\-\/>]?([ACGTN]*)(?:[:_\-\/>]([ACGTN]*))?)?$/);
    if (m) {
      var s = m[1] + '_' + m[2];
      if (m[3]) s += '_' + m[3];
      if (m[4] !== undefined && m[4] !== '') s += '_' + m[4];
      return { t: s, variant: true };
    }
    if (/^ENSG\d+\.\d+$/.test(t)) t = t.replace(/\.\d+$/, '');
    return { t: t, variant: false };
  }
  function search(q) {
    var hits = new Map();                         // row -> rank
    var notes = [];
    function add(i, r) { var o = hits.get(i); if (o === undefined || r < o) hits.set(i, r); }
    function addStudy(s, r) { var rows = D._byStudy[s]; for (var k = 0; k < rows.length; k++) add(rows[k], r); }
    var terms = q.split(/[,;\s]+/).filter(function (x) { return x.length; });
    if (!terms.length) { for (var i = 0; i < N; i++) hits.set(i, null); return { hits: hits, notes: notes }; }
    terms.forEach(function (raw) {
      var nt = normTerm(raw), t = nt.t, n0 = hits.size, found = false;
      if (/^RS\d+$/.test(t)) { notes.push('rsIDs are not in the resource; search “' + esc(raw) + '” by chr_pos_ref_alt (GRCh38).'); return; }
      if (!nt.variant) {
        for (var s = 0; s < studyKeys.length; s++) {
          var keys = studyKeys[s], best = 9;
          for (var k = 0; k < 4; k++) {
            var key = keys[k];
            if (key === t) { best = 0; break; }
            if (key.indexOf(t) === 0) best = Math.min(best, 1);
            else if (t.length >= 3 && key.indexOf(t) > 0) best = Math.min(best, 2);
          }
          if (best < 9) { addStudy(s, best); found = true; }
        }
      }
      // row-level identifiers: lead variant, studyLocusId
      if (/^[0-9XYM]/.test(t) && t.length >= 3) {
        var vv = D._varU;
        for (var r = 0; r < N; r++) {
          if (vv[r] === t) { add(r, 0); found = true; }
          else if (vv[r].indexOf(t) === 0) { add(r, 1); found = true; }
        }
      }
      if (/^[0-9A-F]{6,32}$/.test(t)) {
        var ids = D.studyLocusId, tl = t.toLowerCase();
        for (var r2 = 0; r2 < N; r2++) {
          if (ids[r2] === tl) { add(r2, 0); found = true; } else if (ids[r2].indexOf(tl) === 0) { add(r2, 1); found = true; }
        }
      }
      if (!nt.variant && t.length >= 3 && (!found || state.fuzzy)) {
        var max = t.length <= 5 ? 1 : 2, fz = [];
        for (var u = 0; u < symList.length; u++) {
          var sy = symList[u];
          if (sy === t) continue;
          if (osa(t, sy, max) <= max) { fz.push(sy); symToStudies[sy].forEach(function (s2) { addStudy(s2, 3); }); }
        }
        if (!found) notes.push(fz.length ? 'No exact or prefix match for “' + esc(raw) + '”; showing fuzzy matches: ' + fz.slice(0, 12).map(function (x) { return '<button type="button" class="sugg" data-q="' + esc(x) + '">' + esc(x) + '</button>'; }).join(' ') : 'No match for “' + esc(raw) + '”.');
      }
      if (!found && hits.size === n0 && nt.variant) notes.push('No lead variant starts with “' + esc(t) + '”. Only lead variants are indexed.');
    });
    return { hits: hits, notes: notes };
  }

  /* ---------- filtering, sorting, rendering ---------- */
  var lastNotes = [];
  function run(resetPage) {
    if (!D) return;
    var res = search(state.q);
    lastNotes = res.notes;
    var tiers = {}; state.tiers.forEach(function (t) { tiers[t] = 1; });
    var out = [];
    res.hits.forEach(function (r, i) {
      if (!tiers[v('tier', i)]) return;
      if (state.onlyflag && v('flag_concordant_unanchored_PAV', i) !== 1) return;
      var b = v('flag_beta_inconsistent', i);
      if (state.beta === 'exclude' && b === 1) return;
      if (state.beta === 'only' && b !== 1) return;
      out.push({ i: i, r: r });
    });
    // default order: match rank, symbol, tier, -log10P desc
    out.sort(function (a, b) {
      var ra = a.r === null ? 0 : a.r, rb = b.r === null ? 0 : b.r;
      if (ra !== rb) return ra - rb;
      var sa = st(a.i)[1], sb = st(b.i)[1];
      if (sa !== sb) return sa < sb ? -1 : 1;
      var ta = v('tier', a.i) - v('tier', b.i); if (ta) return ta;
      return v('lead_mlog10p', b.i) - v('lead_mlog10p', a.i);
    });
    if (state.sort) {
      var col = COLS.filter(function (c) { return c.k === state.sort; })[0];
      if (col) {
        out.forEach(function (o, k) { o.o = k; o.sv = col.sort(o.i); });
        var dir = state.dir;
        out.sort(function (a, b) {
          var x = a.sv, y = b.sv, xn = x === null || x === undefined || x === -1, yn2 = y === null || y === undefined || y === -1;
          if (xn !== yn2) return xn ? 1 : -1;         // missing always last
          if (!xn && x !== y) return (x < y ? -1 : 1) * dir;
          return a.o - b.o;
        });
      }
    }
    current = out;
    if (resetPage) state.page = 0;
    render();
  }
  function renderHead() {
    $('thead').innerHTML = COLS.map(function (c) {
      var ind = state.sort === c.k ? (state.dir === 1 ? '▲' : '▼') : '';
      var aria = state.sort === c.k ? (state.dir === 1 ? 'ascending' : 'descending') : 'none';
      return '<th scope="col" aria-sort="' + aria + '" data-k="' + c.k + '"' + (c.num ? ' class="num"' : '') + '><button type="button" class="hbtn" data-sort="' + c.k + '">' + esc(c.h) + '<span class="sortind">' + ind + '</span></button>' +
        '<button type="button" class="info" data-tip="' + esc(c.tip) + '" aria-label="About ' + esc(c.h) + ': ' + esc(c.tip) + '">?</button></th>';
    }).join('');
  }
  function render() {
    renderHead();
    var size = state.size, pages = Math.max(1, Math.ceil(current.length / size));
    if (state.page >= pages) state.page = pages - 1;
    var a = state.page * size, b = Math.min(current.length, a + size), html = [];
    for (var k = a; k < b; k++) {
      var o = current[k], i = o.i, t = v('tier', i);
      html.push('<tr class="t' + t + '" data-row="' + i + '" data-id="' + D.studyLocusId[i] + '" tabindex="0">' + COLS.map(function (c) {
        return '<td' + (c.num ? ' class="num"' : '') + '>' + c.cell(i, o.r) + '</td>';
      }).join('') + '</tr>');
    }
    $('tbody').innerHTML = html.join('');
    var nGenes = {}; current.forEach(function (o) { nGenes[D.s[o.i]] = 1; });
    var nAssay = Object.keys(nGenes).length;
    $('rescount').textContent = fmtInt(current.length) + ' credible set' + (current.length === 1 ? '' : 's') + ' · ' + fmtInt(nAssay) + ' assay' + (nAssay === 1 ? '' : 's');
    $('pageinfo').textContent = current.length ? 'page ' + (state.page + 1) + ' / ' + pages + ' (rows ' + fmtInt(a + 1) + '–' + fmtInt(b) + ')' : '';
    $('prev').disabled = state.page <= 0; $('next').disabled = state.page >= pages - 1; $('csv').disabled = !current.length;
    var msg = state.q ? 'Query “' + esc(state.q) + '”.' : 'No query: browsing all credible sets.';
    if (!current.length) msg += ' No rows match the query and filters.';
    $('status').innerHTML = msg + (lastNotes.length ? ' ' + lastNotes.join(' ') : '');
  }

  /* ---------- detail dialog ---------- */
  function showDetail(i) {
    var s = st(i), h = '<table>';
    GROUPS.forEach(function (g) {
      h += '<tr class="grp"><th colspan="2">' + esc(g[0]) + '</th></tr>';
      g[1].forEach(function (f) {
        var x = rawField(f, i);
        h += '<tr><th>' + esc(f) + '</th><td>' + span(f, x, x === '' ? 'NA' : x, 'mono') + '</td></tr>';
      });
    });
    h += '</table><p class="small">Values are rounded as in the portal data; see the full TSV for all 76 fields at full precision. External links open third-party sites: ' +
      '<a rel="noopener noreferrer" target="_blank" href="https://platform.opentargets.org/credible-set/' + esc(D.studyLocusId[i]) + '">Open Targets credible set</a> · ' +
      '<a rel="noopener noreferrer" target="_blank" href="https://platform.opentargets.org/target/' + esc(s[2]) + '">Open Targets target</a> · ' +
      '<a rel="noopener noreferrer" target="_blank" href="https://www.uniprot.org/uniprotkb/' + esc(s[3]) + '/entry">UniProt</a></p>';
    $('detail-h').textContent = s[1] + ' · ' + D.lead_variant_id[i] + ' · T' + v('tier', i);
    $('dbody').innerHTML = h;
    var dlg = $('detail');
    if (typeof dlg.showModal === 'function') dlg.showModal(); else dlg.setAttribute('open', '');
  }

  /* ---------- CSV ---------- */
  function csvCell(x) { x = raw(x); return /[",\n\r]/.test(x) ? '"' + x.replace(/"/g, '""') + '"' : x; }
  function exportCsv() {
    var head = EXPORT.concat(['match_type']);
    var lines = [head.join(',')];
    current.forEach(function (o) {
      var row = EXPORT.map(function (f) { return csvCell(rawField(f, o.i)); });
      row.push(o.r === null || o.r === undefined ? '' : MATCH_NAME[o.r]);
      lines.push(row.join(','));
    });
    var blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' });
    var url = URL.createObjectURL(blob), a = document.createElement('a');
    var d = new Date(), stamp = d.getFullYear() + String(d.getMonth() + 1).padStart(2, '0') + String(d.getDate()).padStart(2, '0');
    var tag = (state.q || 'all').replace(/[^A-Za-z0-9_\-]+/g, '_').slice(0, 40);
    a.href = url; a.download = 'c10plus_tiers_v' + D._meta.version + '_' + tag + '_' + stamp + '.csv';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }

  /* ---------- URL hash state (no cookies / storage) ---------- */
  function writeHash() {
    var p = new URLSearchParams();
    if (state.q) p.set('q', state.q);
    if (state.tiers.length !== 4) p.set('tier', state.tiers.join(''));
    if (state.onlyflag) p.set('flag', '1');
    if (state.beta !== 'any') p.set('beta', state.beta);
    if (state.fuzzy) p.set('fuzzy', '1');
    if (state.sort) p.set('sort', (state.dir === -1 ? '-' : '') + state.sort);
    var h = p.toString();
    if (('#' + h) !== location.hash && !(h === '' && location.hash === '')) history.replaceState(null, '', h ? '#' + h : location.pathname + location.search);
  }
  function readHash() {
    var p = new URLSearchParams(location.hash.replace(/^#/, ''));
    state.q = p.get('q') || '';
    var t = p.get('tier'); state.tiers = t ? t.split('').map(Number).filter(function (x) { return x >= 1 && x <= 4; }) : [1, 2, 3, 4];
    state.onlyflag = p.get('flag') === '1';
    state.beta = ['any', 'exclude', 'only'].indexOf(p.get('beta')) >= 0 ? p.get('beta') : 'any';
    state.fuzzy = p.get('fuzzy') === '1';
    var s = p.get('sort'); state.sort = s ? s.replace(/^-/, '') : null; state.dir = s && s[0] === '-' ? -1 : 1;
    $('q').value = state.q;
    document.querySelectorAll('input[name=tier]').forEach(function (el) { el.checked = state.tiers.indexOf(+el.value) >= 0; });
    $('onlyflag').checked = state.onlyflag; $('betaflag').value = state.beta; $('fuzzy').checked = state.fuzzy;
  }
  function readForm() {
    state.q = $('q').value.trim();
    state.tiers = []; document.querySelectorAll('input[name=tier]').forEach(function (el) { if (el.checked) state.tiers.push(+el.value); });
    state.onlyflag = $('onlyflag').checked; state.beta = $('betaflag').value; state.fuzzy = $('fuzzy').checked;
    state.size = +$('pagesize').value;
  }

  /* ---------- tooltips ---------- */
  var tip = null, tipOwner = null, tipPinned = false;
  function showTip(el) {
    tip.textContent = el.getAttribute('data-tip'); tip.hidden = false; tipOwner = el;
    var r = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    var x = Math.min(Math.max(8, r.left + r.width / 2 - w / 2), window.innerWidth - w - 8);
    var y = r.bottom + 6; if (y + h > window.innerHeight - 8) y = Math.max(8, r.top - h - 6);
    tip.style.left = x + 'px'; tip.style.top = y + 'px';
  }
  function hideTip() { tip.hidden = true; tipOwner = null; tipPinned = false; }

  /* ---------- init ---------- */
  function fillMeta(meta) {
    document.querySelectorAll('[data-meta]').forEach(function (el) {
      var x = meta[el.getAttribute('data-meta')];
      if (x !== undefined) el.textContent = typeof x === 'number' ? fmtInt(x) : x;
    });
    document.querySelectorAll('[data-tier-count]').forEach(function (el) {
      var t = el.getAttribute('data-tier-count'), n = meta.tier_counts['T' + t];
      el.textContent = 'n = ' + fmtInt(n) + ' (' + (100 * n / meta.n_rows).toFixed(1) + '%)';
    });
  }
  function prepare(J) {
    D = J.data; E = J.enums; S = J.studies.rows; N = D.studyLocusId.length;
    D._meta = J.meta; D._tl = J.meta.tier_labels;
    D._byStudy = S.map(function () { return []; });
    for (var i = 0; i < N; i++) D._byStudy[D.s[i]].push(i);
    D._varU = D.lead_variant_id.map(function (x) { return x.toUpperCase(); });
    studyKeys = S.map(function (s) { return [s[1].toUpperCase(), s[2].toUpperCase(), s[3].toUpperCase(), s[4].toUpperCase()]; });
    S.forEach(function (s, k) { var u = s[1].toUpperCase(); (symToStudies[u] = symToStudies[u] || []).push(k); });
    symList = Object.keys(symToStudies).sort();
    $('genelist').innerHTML = symList.map(function (x) { return '<option value="' + esc(x) + '">'; }).join('');
    fillMeta(J.meta);
  }
  function bind() {
    $('qform').addEventListener('submit', function (e) { e.preventDefault(); readForm(); writeHash(); run(true); });
    $('clear').addEventListener('click', function () { $('q').value = ''; readForm(); writeHash(); run(true); $('q').focus(); });
    document.querySelectorAll('#qform input[type=checkbox], #qform select').forEach(function (el) {
      el.addEventListener('change', function () { readForm(); writeHash(); run(true); });
    });
    document.addEventListener('click', function (e) {
      var ex = e.target.closest('.ex, .sugg');
      if (ex) { $('q').value = ex.getAttribute('data-q'); readForm(); writeHash(); run(true); return; }
      var info = e.target.closest('.info');
      if (info) { e.stopPropagation(); if (tipOwner === info && tipPinned) hideTip(); else { showTip(info); tipPinned = true; } return; }
      if (!tip.hidden) hideTip();
      var sb = e.target.closest('[data-sort]');
      if (sb) { var k = sb.getAttribute('data-sort'); if (state.sort === k) state.dir = -state.dir; else { state.sort = k; state.dir = 1; } writeHash(); run(false); return; }
      var tr = e.target.closest('#tbody tr');
      if (tr && !e.target.closest('a')) showDetail(+tr.getAttribute('data-row'));
    });
    document.addEventListener('mouseover', function (e) { var info = e.target.closest && e.target.closest('.info'); if (info && !tipPinned && window.matchMedia('(hover: hover)').matches) showTip(info); });
    document.addEventListener('mouseout', function (e) { var info = e.target.closest && e.target.closest('.info'); if (info && !tipPinned && window.matchMedia('(hover: hover)').matches) hideTip(); });
    document.addEventListener('focusin', function (e) { if (e.target.classList && e.target.classList.contains('info')) showTip(e.target); });
    document.addEventListener('focusout', function (e) { if (e.target.classList && e.target.classList.contains('info') && !tipPinned) hideTip(); });
    $('tbody').addEventListener('keydown', function (e) { if (e.key === 'Enter') { var tr = e.target.closest('tr'); if (tr) showDetail(+tr.getAttribute('data-row')); } });
    $('dclose').addEventListener('click', function () { $('detail').close ? $('detail').close() : $('detail').removeAttribute('open'); });
    $('detail').addEventListener('click', function (e) { if (e.target === $('detail')) $('detail').close(); });
    $('prev').addEventListener('click', function () { state.page--; render(); $('tablewrap').scrollTop = 0; });
    $('next').addEventListener('click', function () { state.page++; render(); $('tablewrap').scrollTop = 0; });
    $('csv').addEventListener('click', exportCsv);
    $('tablewrap').addEventListener('scroll', function () { if (!tip.hidden) hideTip(); });
    window.addEventListener('hashchange', function () { readHash(); run(true); });
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape' && !tip.hidden) hideTip(); });
  }
  function init() {
    tip = $('tip');
    renderHead();
    bind();
    readHash();
    fetch(DATA_URL, { credentials: 'omit', cache: 'default' })
      .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
      .then(function (J) { prepare(J); readForm(); run(true); document.documentElement.setAttribute('data-ready', '1'); })
      .catch(function (err) {
        $('status').textContent = 'Could not load the data file (' + err.message + '). If you opened this page from disk, serve the folder over HTTP (e.g. python -m http.server) or use the hosted portal.';
      });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
