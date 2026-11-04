/* pQTL Evidence Tiers — home page: lazy gene-symbol suggestions from a small local index. MIT License. */
(function () {
  'use strict';
  var input = document.getElementById('q'), list = document.getElementById('genes'), loaded = false;
  if (!input || !list) return;
  input.addEventListener('focus', function () {
    if (loaded) return; loaded = true;
    fetch('data/genes_v2.3.json', { credentials: 'omit' })
      .then(function (r) { return r.ok ? r.json() : []; })
      .then(function (g) { list.innerHTML = g.map(function (x) { return '<option value="' + String(x).replace(/[&<>"']/g, '') + '">'; }).join(''); })
      .catch(function () { /* suggestions are optional */ });
  });
})();
