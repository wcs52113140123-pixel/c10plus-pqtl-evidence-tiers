/* pQTL Evidence Tiers — result page. Vanilla JS, no dependencies; the only request is the local data file. MIT License. */
(function () {
  'use strict';
  var DATA_URL = 'data/tiers_v2.3.json';
  var PAGE_CARDS = 40;
  var D, E, S, N, META;
  var byStudy, varU, rsMap, studyKeys, symList, symToStudies;
  var state = { q: '', tiers: [1, 2, 3, 4], flag: false, sort: 'p' };
  var current = [], groups = [], shown = 0, notes = [];
  var $ = function (id) { return document.getElementById(id); };

  var TIER_LABEL = { 1: 'Anchored, no PAV', 2: 'Anchored with PAV, or isoform-only', 3: 'Unresolved: no anchor, no PAV', 4: 'Unanchored, with PAV' };
  var MS_TEXT = { supported: 'Supported', refuted: 'Refuted', indeterminate_underpowered: 'Indeterminate', indeterminate: 'Indeterminate', not_assessed: 'Not assessed', not_evaluable: 'Not evaluable' };
  var MATCH = ['exact', 'prefix', 'contains', 'fuzzy'];
  var BINARY = { anchor_eqtl: 1, anchor_isoform: 1, pav_cognate_cs95_vep116_alltx: 1, pav_cognate_cs95_canonical: 1, pav_lead_variant: 1, rep_sun2018: 1, rep_decode_v4: 1,
    rep_aric_ea: 1, rep_aric_aa: 1, rep_fenland: 1, xa_replicated: 1, flag_concordant_unanchored_PAV: 1, flag_beta_inconsistent: 1 };
  var GROUPS = [
    ['Identifiers', [['studyLocusId', 'Credible set ID'], ['studyId', 'Study ID'], ['gene_id', 'Ensembl'], ['uniprot', 'UniProt'], ['olink_oid', 'Olink ID'], ['chrom', 'Chromosome'],
      ['lead_variant_id', 'Lead variant (GRCh38)'], ['rsIds', 'rsID'], ['tier', 'Tier'], ['tier_label', 'Tier label']]],
    ['Lead signal', [['lead_beta', 'β (NPX)'], ['lead_mlog10p', '−log10 P'], ['lead_maf_gnomad_nfe', 'MAF (gnomAD NFE)'], ['n_variants_cs95', 'Variants in 95% set'],
      ['lead_beta_to_expected_ratio', '|β| / expected'], ['flag_beta_inconsistent', 'β inconsistent flag'], ['n_other_genes_trans_cs_shared', 'Trans-hotspot proteins']]],
    ['Transcript anchoring', [['anchor_eqtl', 'eQTL anchor'], ['n_anchor_eqtl_datasets', 'Anchoring eQTL datasets'], ['h4max_eqtl_gtex_v10', 'PP.H4, GTEx v10'],
      ['h4max_eqtl_eqtlcatalogue', 'PP.H4, eQTL Catalogue'], ['anchor_isoform', 'Isoform anchor'], ['h4max_sqtl_tuqtl', 'PP.H4, sQTL/tuQTL']]],
    ['Protein-altering variants', [['pav_cognate_cs95_vep116_alltx', 'PAV, any transcript'], ['pav_cognate_cs95_canonical', 'PAV, canonical'], ['pav_lead_variant', 'Lead variant is PAV'], ['pav_pip_sum', 'PAV PIP sum']]],
    ['Affinity replication', [['rep_sun2018', 'Sun 2018 (SomaScan v3)'], ['dir_sun2018', 'Sun 2018 direction'], ['rep_decode_v4', 'deCODE (v4)'], ['dir_decode_v4', 'deCODE direction'],
      ['rep_aric_ea', 'ARIC-EA (v4)'], ['dir_aric_ea', 'ARIC-EA direction'], ['rep_aric_aa', 'ARIC-AA (v4)'], ['rep_fenland', 'Fenland (v4)'], ['dir_fenland', 'Fenland direction'],
      ['platform_concordance', 'Platform concordance'], ['flag_concordant_unanchored_PAV', 'Concordant + unanchored + PAV']]],
    ['Mass spectrometry', [['ms_primary', 'MSPA label'], ['ms_mspa_peptide', 'MSPA peptide label'], ['ms_seer_discovery', 'Seer discovery overlap'], ['ms_ext_primary', 'Seer full-statistics label'],
      ['ms_ext_nominal', 'Seer nominal label'], ['seer_ext_p', 'Seer P'], ['seer_ext_opposite_sig', 'Seer opposite and significant'],
      ['ms_gs_primary', 'Generation Scotland label'], ['ms_gs_nominal', 'Generation Scotland nominal label'], ['gs_p', 'Generation Scotland P'], ['gs_opposite_sig', 'Generation Scotland opposite and significant'],
      ['hb_evaluable', 'HOLBAEK statistics released'], ['ms_hb_primary', 'HOLBAEK label'], ['ms_hb_nominal', 'HOLBAEK nominal label'], ['hb_p', 'HOLBAEK P'], ['hb_opposite_sig', 'HOLBAEK opposite and significant'],
      ['ms_independent_combined', 'Independent MS layers combined']]],
    ['Cross-ancestry (UKB-PPP)', [['xa_n_ancestries', 'Ancestries'], ['xa_meta_beta_aligned', 'Meta β, aligned'], ['xa_meta_p', 'Meta P'], ['xa_replicated', 'Replicated']]]
  ];
  var EXPORT = ['studyLocusId', 'studyId', 'symbol', 'gene_id', 'uniprot', 'olink_oid', 'chrom', 'lead_variant_id', 'rsIds', 'lead_beta', 'lead_mlog10p',
    'lead_maf_gnomad_nfe', 'n_variants_cs95', 'tier', 'tier_label', 'anchor_eqtl', 'n_anchor_eqtl_datasets', 'h4max_eqtl_gtex_v10', 'h4max_eqtl_eqtlcatalogue',
    'anchor_isoform', 'h4max_sqtl_tuqtl', 'pav_cognate_cs95_vep116_alltx', 'pav_cognate_cs95_canonical', 'pav_lead_variant', 'pav_pip_sum',
    'rep_sun2018', 'dir_sun2018', 'rep_decode_v4', 'dir_decode_v4', 'rep_aric_ea', 'dir_aric_ea', 'rep_aric_aa', 'rep_fenland', 'dir_fenland',
    'platform_concordance', 'ms_primary', 'ms_mspa_peptide', 'ms_seer_discovery', 'ms_ext_primary', 'ms_ext_nominal', 'seer_ext_p', 'seer_ext_opposite_sig',
    'ms_gs_primary', 'ms_gs_nominal', 'gs_p', 'gs_opposite_sig', 'hb_evaluable', 'ms_hb_primary', 'ms_hb_nominal', 'hb_p', 'hb_opposite_sig', 'ms_independent_combined',
    'xa_n_ancestries', 'xa_meta_beta_aligned', 'xa_meta_p', 'xa_replicated', 'flag_concordant_unanchored_PAV', 'lead_beta_to_expected_ratio',
    'flag_beta_inconsistent', 'n_other_genes_trans_cs_shared'];
  var STUDY_F = { studyId: 0, symbol: 1, gene_id: 2, uniprot: 3, olink_oid: 4, chrom: 5 };

  function esc(s) { return String(s).replace(/[&<>"']/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]; }); }
  function v(f, i) { return D[f][i]; }
  function ev(f, i) { return E[f][D[f][i]]; }
  function st(i) { return S[D.s[i]]; }
  function fmtInt(n) { return Number(n).toLocaleString('en-US'); }
  function raw(x) { return x === null || x === undefined ? '' : String(x); }
  function rawField(f, i) {
    if (f in STUDY_F) return st(i)[STUDY_F[f]];
    if (f === 'tier_label') return META.tier_labels[v('tier', i)];
    if (E[f]) return ev(f, i);
    var x = D[f][i];
    return x === null || x === undefined ? '' : x;
  }
  function sp(f, val, text, cls) { return '<span data-f="' + f + '" data-v="' + esc(raw(val)) + '"' + (cls ? ' class="' + cls + '"' : '') + '>' + esc(text) + '</span>'; }
  function fmtP(x) { return x === 0 ? '0' : (Math.abs(x) < 1e-3 ? x.toExponential(2) : x.toPrecision(3)); }

  /* ---------- search ---------- */
  function osa(a, b, max) {
    var la = a.length, lb = b.length, i, j;
    if (Math.abs(la - lb) > max) return max + 1;
    var d = [];
    for (i = 0; i <= la; i++) d[i] = [i];
    for (j = 0; j <= lb; j++) d[0][j] = j;
    for (i = 1; i <= la; i++) {
      var m = Infinity;
      for (j = 1; j <= lb; j++) {
        var c = a[i - 1] === b[j - 1] ? 0 : 1;
        var x = Math.min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + c);
        if (i > 1 && j > 1 && a[i - 1] === b[j - 2] && a[i - 2] === b[j - 1]) x = Math.min(x, d[i - 2][j - 2] + 1);
        d[i][j] = x; if (x < m) m = x;
      }
      if (m > max) return max + 1;
    }
    return d[la][lb];
  }
  function normTerm(t) {
    t = t.trim().toUpperCase();
    var m = t.match(/^(?:CHR)?([0-9]{1,2}|X|Y|MT?)[:_\-]([0-9]+)(?:[:_\-\/>]?([ACGTN]*)(?:[:_\-\/>]([ACGTN]*))?)?$/);
    if (m) {
      var s = m[1] + '_' + m[2];
      if (m[3]) s += '_' + m[3];
      if (m[4]) s += '_' + m[4];
      return { t: s, kind: 'variant' };
    }
    if (/^RS\d+$/.test(t)) return { t: t.toLowerCase(), kind: 'rs' };
    if (/^ENSG\d+\.\d+$/.test(t)) t = t.replace(/\.\d+$/, '');
    return { t: t, kind: 'id' };
  }
  function search(q) {
    var hits = new Map(); notes = [];
    function add(i, r) { var o = hits.get(i); if (o === undefined || r < o) hits.set(i, r); }
    function addStudy(s, r) { byStudy[s].forEach(function (i) { add(i, r); }); }
    var terms = q.split(/[,;\s]+/).filter(function (x) { return x.length; });
    if (!terms.length) { if (state.flag) for (var i = 0; i < N; i++) hits.set(i, null); return hits; }
    terms.forEach(function (rawT) {
      var nt = normTerm(rawT), t = nt.t, found = false, r, s, k;
      if (nt.kind === 'rs') {
        (rsMap[t] || []).forEach(function (i) { add(i, 0); found = true; });
        if (!found) notes.push('No lead variant has rsID “' + esc(rawT) + '”. Only lead variants of credible sets are indexed.');
        return;
      }
      if (nt.kind === 'id') {
        for (s = 0; s < studyKeys.length; s++) {
          var keys = studyKeys[s], best = 9;
          for (k = 0; k < 4; k++) {
            var key = keys[k];
            if (key === t) { best = 0; break; }
            if (key.indexOf(t) === 0) best = Math.min(best, 1);
            else if (t.length >= 3 && key.indexOf(t) > 0) best = Math.min(best, 2);
          }
          if (best < 9) { addStudy(s, best); found = true; }
        }
        if (/^[0-9A-F]{6,32}$/.test(t)) {
          var tl = t.toLowerCase(), ids = D.studyLocusId;
          for (r = 0; r < N; r++) { if (ids[r] === tl) { add(r, 0); found = true; } else if (ids[r].indexOf(tl) === 0) { add(r, 1); found = true; } }
        }
      }
      if (nt.kind === 'variant' || (/^[0-9XYM]/.test(t) && t.length >= 3)) {
        for (r = 0; r < N; r++) {
          if (varU[r] === t) { add(r, 0); found = true; } else if (varU[r].indexOf(t) === 0) { add(r, 1); found = true; }
        }
        if (!found && nt.kind === 'variant') notes.push('No lead variant starts with “' + esc(t) + '”.');
      }
      if (nt.kind === 'id' && !found && t.length >= 3) {
        var max = t.length <= 5 ? 1 : 2, fz = [];
        symList.forEach(function (sy) { if (osa(t, sy, max) <= max) { fz.push(sy); symToStudies[sy].forEach(function (s2) { addStudy(s2, 3); }); } });
        notes.push(fz.length ? 'No exact match for “' + esc(rawT) + '”. Showing close matches: ' + fz.slice(0, 8).map(function (x) {
          return '<button type="button" class="sugg" data-q="' + esc(x) + '">' + esc(x) + '</button>'; }).join(', ') + '.' : 'No match for “' + esc(rawT) + '”.');
      }
    });
    return hits;
  }

  /* ---------- state <-> URL ---------- */
  function readURL() {
    var p = new URLSearchParams(location.search);
    state.q = (p.get('q') || '').trim();
    var t = p.get('tier');
    state.tiers = t ? t.split('').map(Number).filter(function (x) { return x >= 1 && x <= 4; }) : [1, 2, 3, 4];
    if (!state.tiers.length) state.tiers = [1, 2, 3, 4];
    state.flag = p.get('flag') === '1';
    state.sort = p.get('sort') === 'tier' ? 'tier' : 'p';
  }
  function urlFor() {
    var p = new URLSearchParams();
    if (state.q) p.set('q', state.q);
    if (state.tiers.length !== 4) p.set('tier', state.tiers.join(''));
    if (state.flag) p.set('flag', '1');
    if (state.sort !== 'p') p.set('sort', state.sort);
    var s = p.toString();
    return location.pathname + (s ? '?' + s : '');
  }
  function syncControls() {
    $('q').value = state.q;
    document.querySelectorAll('[data-tier]').forEach(function (b) { b.setAttribute('aria-pressed', state.tiers.indexOf(+b.getAttribute('data-tier')) >= 0 ? 'true' : 'false'); });
    $('flagonly').setAttribute('aria-pressed', state.flag ? 'true' : 'false');
    $('sort').value = state.sort;
  }

  /* ---------- run + render ---------- */
  function run() {
    syncControls();
    var hits = search(state.q), tiers = {}, out = [];
    state.tiers.forEach(function (t) { tiers[t] = 1; });
    var nAll = hits.size;
    hits.forEach(function (r, i) {
      if (!tiers[v('tier', i)]) return;
      if (state.flag && v('flag_concordant_unanchored_PAV', i) !== 1) return;
      out.push({ i: i, r: r === null ? 0 : r, m: r });
    });
    var g = {};
    out.forEach(function (o) { var s = D.s[o.i]; if (!g[s]) g[s] = { s: s, r: o.r, rows: [] }; g[s].r = Math.min(g[s].r, o.r); g[s].rows.push(o); });
    groups = Object.keys(g).map(function (k) { return g[k]; });
    groups.sort(function (a, b) { return a.r - b.r || (S[a.s][1] < S[b.s][1] ? -1 : S[a.s][1] > S[b.s][1] ? 1 : 0); });
    groups.forEach(function (gr) {
      gr.rows.sort(state.sort === 'tier'
        ? function (a, b) { return v('tier', a.i) - v('tier', b.i) || v('lead_mlog10p', b.i) - v('lead_mlog10p', a.i); }
        : function (a, b) { return v('lead_mlog10p', b.i) - v('lead_mlog10p', a.i) || Math.abs(v('lead_beta', b.i)) - Math.abs(v('lead_beta', a.i)); });
    });
    current = []; groups.forEach(function (gr, gi) { gr.rows.forEach(function (o) { o.g = gi; current.push(o); }); });
    $('results').innerHTML = ''; shown = 0; lastG = -1; lastSec = null;
    var hasQuery = !!state.q || state.flag;
    $('toolbar').hidden = !hasQuery;
    if (!hasQuery) {
      $('status').hidden = false;
      $('status').innerHTML = 'Enter a gene symbol, UniProt accession, rsID, variant (chr_pos_ref_alt) or Olink ID. Examples: ' +
        ['PCSK9', 'APOE', 'IL6R', 'LPA'].map(function (x) { return '<a href="?q=' + x + '">' + x + '</a>'; }).join(', ') + '.';
      $('more').hidden = true; $('count').textContent = '';
      document.title = 'Search · pQTL Evidence Tiers';
      return;
    }
    var msg = current.length ? '' : (nAll ? 'No credible sets match the selected filters.' : '');
    $('status').innerHTML = (msg + ' ' + notes.join(' ')).trim();
    $('status').hidden = !$('status').innerHTML;
    $('count').textContent = fmtInt(current.length) + ' credible set' + (current.length === 1 ? '' : 's') + ' · ' + fmtInt(groups.length) + ' protein' + (groups.length === 1 ? '' : 's');
    $('csv').disabled = !current.length;
    document.title = (state.q || 'Flagged') + ' · pQTL Evidence Tiers';
    renderMore();
    if (current.length === 1) { var d = document.querySelector('#results details.more'); if (d) d.open = true; }
  }
  var lastG = -1, lastSec = null;
  function renderMore() {
    var end = Math.min(current.length, shown + PAGE_CARDS), buf = [], k;
    function flush() { if (buf.length && lastSec) { lastSec.insertAdjacentHTML('beforeend', buf.join('')); buf = []; } }
    for (k = shown; k < end; k++) {
      var o = current[k];
      if (o.g !== lastG) {
        flush();
        $('results').insertAdjacentHTML('beforeend', groupHead(groups[o.g]));
        lastSec = $('results').lastElementChild; lastG = o.g;
        lastSec.querySelectorAll('.tierbar [data-w]').forEach(function (el) { el.style.width = el.getAttribute('data-w') + '%'; });
      }
      buf.push(card(o));
    }
    flush();
    shown = end;
    var rest = current.length - shown;
    $('more').hidden = rest <= 0;
    $('more').textContent = 'Show more (' + fmtInt(rest) + ' remaining)';
  }
  function groupHead(gr) {
    var s = S[gr.s], cnt = { 1: 0, 2: 0, 3: 0, 4: 0 }, all = byStudy[gr.s];
    all.forEach(function (i) { cnt[v('tier', i)]++; });
    var n = all.length, bar = [1, 2, 3, 4].map(function (t) { return cnt[t] ? '<span class="t' + t + '" data-w="' + (100 * cnt[t] / n).toFixed(2) + '"></span>' : ''; }).join('');
    var sum = [1, 2, 3, 4].filter(function (t) { return cnt[t]; }).map(function (t) { return 'Tier ' + t + ': ' + cnt[t]; }).join(' · ');
    var shownN = gr.rows.length;
    return '<section class="protein" aria-labelledby="p-' + gr.s + '"><div class="protein-head"><h2 id="p-' + gr.s + '">' + esc(s[1]) + '</h2>' +
      '<span class="ids"><span>UniProt <a href="https://www.uniprot.org/uniprotkb/' + esc(s[3]) + '/entry" rel="noopener noreferrer">' + esc(s[3]) + '</a></span>' +
      '<span>' + esc(s[2]) + '</span><span>' + esc(s[4]) + '</span><span>chr' + esc(s[5]) + '</span></span></div>' +
      '<div class="tierbar" role="img" aria-label="All ' + n + ' credible sets of this protein by tier: ' + sum + '">' + bar + '</div>' +
      '<p class="tiersum">' + n + ' credible set' + (n === 1 ? '' : 's') + ' (' + sum + ')' + (shownN !== n ? ' · ' + shownN + ' shown' : '') + '</p></section>';
  }
  function repRow(name, f, dirf, i) {
    var r = v(f, i), t = r === null ? 'n/a' : (r === 1 ? 'Yes' : 'No');
    var h = sp(f, r, t, r === null ? 'na' : (r === 1 ? 'yes' : 'no'));
    if (dirf && r === 1) { var d = ev(dirf, i); if (d !== 'NA') h += ' <span class="' + (d === 'opposite' ? 'opp' : 'no') + '">· ' + sp(dirf, d, d) + '</span>'; }
    return '<div class="row"><span class="src">' + name + '</span><span>' + h + '</span></div>';
  }
  function msRow(name, f, i) {
    var m = ev(f, i);
    return '<div class="row"><span class="src">' + name + '</span>' + sp(f, m, MS_TEXT[m], m === 'supported' ? 'yes' : (m === 'not_assessed' || m === 'not_evaluable' ? 'na' : 'no')) + '</div>';
  }
  function card(o) {
    var i = o.i, t = v('tier', i), id = D.studyLocusId[i], rs = D.rsIds[i];
    var pav = v('pav_cognate_cs95_vep116_alltx', i), anc = v('anchor_eqtl', i);
    var pavH = sp('pav_cognate_cs95_vep116_alltx', pav, pav === 1 ? 'Yes' : 'No', pav === 1 ? 'yes' : 'no') +
      (pav === 1 && v('pav_lead_variant', i) === 1 ? ' <span class="no">· ' + sp('pav_lead_variant', 1, 'lead variant') + '</span>' : '');
    var ancH = sp('anchor_eqtl', anc, anc === 1 ? 'Yes' : 'No', anc === 1 ? 'yes' : 'no') +
      (anc === 1 ? ' <span class="no">· ' + sp('n_anchor_eqtl_datasets', v('n_anchor_eqtl_datasets', i), v('n_anchor_eqtl_datasets', i)) + ' dataset' + (v('n_anchor_eqtl_datasets', i) === 1 ? '' : 's') + '</span>' : '') +
      (v('anchor_isoform', i) === 1 ? '<div class="no">' + sp('anchor_isoform', 1, 'isoform anchor') + '</div>' : '');
    var flags = [];
    if (v('flag_concordant_unanchored_PAV', i) === 1) flags.push('<span class="flagchip">' + sp('flag_concordant_unanchored_PAV', 1, 'Concordant + unanchored + PAV') + '</span>');
    if (v('flag_beta_inconsistent', i) === 1) flags.push('<span class="flagchip">' + sp('flag_beta_inconsistent', 1, 'β inconsistent') + '</span>');
    var mtag = o.m !== null && o.m > 0 ? ' <span class="no">(' + MATCH[o.m] + ' match)</span>' : '';
    return '<article class="card t' + t + '" id="cs-' + id + '" data-id="' + id + '" aria-label="' + esc(st(i)[1] + ' credible set, lead variant ' + D.lead_variant_id[i] + ', tier ' + t) + '">' +
      '<div class="card-top"><span class="badge t' + t + '">' + sp('tier', t, 'Tier ' + t) + '</span><span class="tlabel">' + sp('tier_label', META.tier_labels[t], TIER_LABEL[t]) + '</span>' +
      sp('lead_variant_id', D.lead_variant_id[i], D.lead_variant_id[i], 'variant') + (rs ? sp('rsIds', rs, rs.replace(/;/g, ', '), 'rs') : '') + mtag +
      '<span class="sig">−log10 P ' + sp('lead_mlog10p', v('lead_mlog10p', i), v('lead_mlog10p', i).toFixed(1)) + ' · β ' + sp('lead_beta', v('lead_beta', i), v('lead_beta', i).toPrecision(3)) +
      ' · <a class="plink" href="?q=' + id + '" aria-label="Permalink to this credible set">permalink</a></span></div>' +
      '<div class="evidence">' +
      '<div class="ev"><span class="k">PAV</span><div class="v">' + pavH + '</div></div>' +
      '<div class="ev"><span class="k">eQTL anchor</span><div class="v">' + ancH + '</div></div>' +
      '<div class="ev"><span class="k">Replication</span>' + repRow('deCODE', 'rep_decode_v4', 'dir_decode_v4', i) + repRow('ARIC-EA', 'rep_aric_ea', 'dir_aric_ea', i) + repRow('Fenland', 'rep_fenland', 'dir_fenland', i) + '</div>' +
      '<div class="ev ms"><span class="k">MS support</span>' + msRow('MSPA', 'ms_primary', i) + msRow('Seer', 'ms_ext_primary', i) + msRow('Gen. Scotland', 'ms_gs_primary', i) + msRow('HOLBAEK', 'ms_hb_primary', i) +
        (st(i)[5] === 'X' ? '<div class="row"><span class="na">chrX: no MS layer releases chrX data</span></div>' : '') + '</div>' +
      '<div class="ev"><span class="k">Flags</span><div class="v">' + (flags.length ? flags.join('') : '<span class="na">None</span>') + '</div></div>' +
      '</div><details class="more"><summary>Details</summary><div class="kv" data-row="' + i + '"></div></details></article>';
  }
  function detailHTML(i) {
    return GROUPS.map(function (g) {
      return '<section><h3>' + esc(g[0]) + '</h3><dl>' + g[1].map(function (p) {
        var f = p[0], x = rawField(f, i), txt;
        if (x === '' || x === 'NA') txt = '—';
        else if (BINARY[f]) txt = +x === 1 ? 'yes' : 'no';
        else txt = String(x);
        return '<div><dt title="' + f + '">' + esc(p[1]) + '</dt><dd>' + sp(f, x, txt) + '</dd></div>';
      }).join('') + '</dl></section>';
    }).join('');
  }

  /* ---------- CSV ---------- */
  function csvCell(x) { x = raw(x); return /[",\n\r]/.test(x) ? '"' + x.replace(/"/g, '""') + '"' : x; }
  function exportCsv() {
    var lines = [EXPORT.concat(['match_type']).join(',')];
    current.forEach(function (o) {
      var row = EXPORT.map(function (f) { return csvCell(rawField(f, o.i)); });
      row.push(o.m === null ? '' : MATCH[o.m]);
      lines.push(row.join(','));
    });
    var blob = new Blob([lines.join('\n') + '\n'], { type: 'text/csv;charset=utf-8' });
    var a = document.createElement('a'), url = URL.createObjectURL(blob);
    var tag = (state.q || (state.flag ? 'flagged' : 'all')).replace(/[^A-Za-z0-9_\-]+/g, '_').slice(0, 40);
    a.href = url; a.download = 'pqtl_evidence_tiers_v' + META.version + '_' + tag + '.csv';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 2000);
  }
  function copyLink() {
    var url = location.href, b = $('copy');
    function done(ok) { b.textContent = ok ? 'Link copied' : 'Copy failed'; setTimeout(function () { b.textContent = 'Copy link'; }, 1800); }
    if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(url).then(function () { done(true); }, function () { done(false); });
    else { var ta = document.createElement('textarea'); ta.value = url; document.body.appendChild(ta); ta.select(); var ok = false; try { ok = document.execCommand('copy'); } catch (e) { ok = false; } ta.remove(); done(ok); }
  }

  /* ---------- init ---------- */
  function prepare(J) {
    D = J.data; E = J.enums; S = J.studies.rows; META = J.meta; N = D.studyLocusId.length;
    byStudy = S.map(function () { return []; });
    for (var i = 0; i < N; i++) byStudy[D.s[i]].push(i);
    varU = D.lead_variant_id.map(function (x) { return x.toUpperCase(); });
    rsMap = {};
    D.rsIds.forEach(function (x, i) { if (x) x.split(';').forEach(function (r) { (rsMap[r] = rsMap[r] || []).push(i); }); });
    studyKeys = S.map(function (s) { return [s[1].toUpperCase(), s[2].toUpperCase(), s[3].toUpperCase(), s[4].toUpperCase()]; });
    symToStudies = {};
    S.forEach(function (s, k) { var u = s[1].toUpperCase(); (symToStudies[u] = symToStudies[u] || []).push(k); });
    symList = Object.keys(symToStudies).sort();
  }
  function navigate(push) {
    var u = urlFor();
    if (push) history.pushState(null, '', u); else history.replaceState(null, '', u);
    run();
  }
  function bind() {
    $('qform').addEventListener('submit', function (e) { e.preventDefault(); state.q = $('q').value.trim(); navigate(true); });
    document.querySelectorAll('[data-tier]').forEach(function (b) {
      b.addEventListener('click', function () {
        var t = +b.getAttribute('data-tier'), k = state.tiers.indexOf(t);
        if (k >= 0 && state.tiers.length > 1) state.tiers.splice(k, 1); else if (k < 0) { state.tiers.push(t); state.tiers.sort(); }
        navigate(false);
      });
    });
    $('flagonly').addEventListener('click', function () { state.flag = !state.flag; navigate(false); });
    $('sort').addEventListener('change', function () { state.sort = $('sort').value; navigate(false); });
    $('csv').addEventListener('click', exportCsv);
    $('copy').addEventListener('click', copyLink);
    $('more').addEventListener('click', function () { renderMore(); });
    $('results').addEventListener('toggle', function (e) {
      var d = e.target; if (d.tagName !== 'DETAILS' || !d.open) return;
      var kv = d.querySelector('.kv'); if (kv && !kv.innerHTML) kv.innerHTML = detailHTML(+kv.getAttribute('data-row'));
    }, true);
    document.addEventListener('click', function (e) {
      var s = e.target.closest('.sugg'); if (s) { state.q = s.getAttribute('data-q'); navigate(true); }
    });
    window.addEventListener('popstate', function () { readURL(); run(); });
  }
  function init() {
    readURL(); syncControls(); bind();
    fetch(DATA_URL, { credentials: 'omit' })
      .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
      .then(function (J) { prepare(J); run(); document.documentElement.setAttribute('data-ready', '1'); })
      .catch(function (err) { $('status').textContent = 'The data file could not be loaded (' + err.message + '). Serve this folder over HTTP, or use the hosted site.'; });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
