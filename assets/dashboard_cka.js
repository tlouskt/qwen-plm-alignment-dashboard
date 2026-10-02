/* CKA main report: saved full-panel scores and signed row allocations only. */
(() => {
  'use strict';
  const D = window.ALIGNMENT_DATA, $ = id => document.getElementById(id);
  const ids = D.panel.protein_ids, ann = D.panel.annotations, N = ids.length;
  const indices = ids.map((_, i) => i), byId = new Map(ids.map((id, i) => [id, i]));
  const qs = [...D.qwen].sort((a, b) => a.size_b - b.size_b), plms = D.plms;
  const modelNames = {prot_t5_xl_uniref50: 'ProtT5 XL', esmc_300m: 'ESMC 300M', esmc_600m: 'ESMC 600M', esm3_sm_open_v1: 'ESM3 small', esm2_t6_8m: 'ESM2 8M', esm2_t12_35m: 'ESM2 35M', esm2_t30_150m: 'ESM2 150M', esm2_t33_650m: 'ESM2 650M', esmfold_v1: 'ESMFold v1'};
  const colors = ['#0072b2', '#d55e00', '#009e73', '#cc79a7', '#7057a3', '#7a6a31', '#4b7f86', '#b34b59', '#696969'];
  const fullName = id => plms.find(p => p.id === id)?.label || id;
  const shortName = id => modelNames[id] || fullName(id);
  const color = id => colors[plms.findIndex(p => p.id === id) % colors.length];
  const esc = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const number = value => Number.isFinite(value) ? value.toLocaleString('en-US', {maximumSignificantDigits: 3}) : 'Undefined';
  const signed = value => (value > 0 ? '+' : '') + number(value);
  const sum = values => values.reduce((total, value) => total + value, 0);
  const genes = indices.map(i => String(ann.gene_name?.[i] ?? '').trim()), geneCounts = new Map();
  genes.forEach(gene => { if (gene) geneCounts.set(gene, (geneCounts.get(gene) || 0) + 1); });
  const label = i => genes[i] ? geneCounts.get(genes[i]) > 1 ? `${genes[i]} (${ids[i]})` : genes[i] : ids[i];
  let latest, active = false, pending = false, task = Promise.resolve();
  let heatmapCohort = 'table';
  const api = window.CKA_REPORT = {ready: false, error: null, render};

  function family(i) {
    const raw = String(ann.group?.[i] ?? '').trim();
    const id = !raw || /^(unknown|none|null|nan)$/i.test(raw) ? 'unknown' : raw;
    const entry = D.panel.group_names?.[id];
    const short = entry?.name || ann.group_name?.[i] || '';
    const name = id === 'unknown' ? 'No recorded Pfam family' : entry?.description || entry?.name || ann.group_description?.[i] || ann.group_name?.[i] || 'Unlabeled Pfam family';
    return {id, name, short, label: name};
  }

  function topData(state) {
    const record = D.comparisons[`${state.qwen}|${state.a}`]?.cka;
    if (!record || !Number.isFinite(record.global) || record.contribution?.length !== N || record.contribution.some(value => !Number.isFinite(value))) throw Error('Complete saved CKA scores and signed contributions are required.');
    const requested = [10, 20, 50, 100].includes(Number(state.topN)) ? Number(state.topN) : 20;
    const positive = indices.filter(i => record.contribution[i] > 0).sort((i, j) => record.contribution[j] - record.contribution[i] || ids[i].localeCompare(ids[j]));
    const cutoff = positive.length ? record.contribution[positive[Math.min(requested, positive.length) - 1]] : null;
    const chosen = cutoff === null ? [] : positive.filter(i => record.contribution[i] >= cutoff);
    const groups = new Map();
    const rows = chosen.map((i, j) => {
      const f = family(i);
      if (!groups.has(f.id)) groups.set(f.id, {...f, n: 0});
      groups.get(f.id).n++;
      const rank = j && record.contribution[i] === record.contribution[chosen[j - 1]] ? null : j + 1;
      return {i, contribution: record.contribution[i], familyId: f.id, rank};
    });
    rows.forEach((row, j) => { if (row.rank === null) row.rank = rows[j - 1].rank; });
    const repeatedNames = new Map();
    for (const group of groups.values()) repeatedNames.set(group.label, (repeatedNames.get(group.label) || 0) + 1);
    for (const group of groups.values()) if (repeatedNames.get(group.label) > 1 && group.short) group.label += ` — ${group.short}`;
    const families = [...groups.values()].sort((a, b) => b.n - a.n || a.id.localeCompare(b.id));
    if (sum(families.map(group => group.n)) !== rows.length) throw Error('Family counts must cover the exact CKA top-protein cohort.');
    const total = sum(record.contribution), tolerance = Number(record.validation?.tolerance ?? 1e-10);
    if (Math.abs(total - record.global) > tolerance) throw Error('Saved signed CKA contributions do not reconstruct the global score.');
    const validShare = record.global > 0 && record.raw_denominator > 0;
    return {record, requested, cutoff, rows, families, groups, positive: positive.length,
      negative: record.contribution.filter(value => value < 0).length,
      zero: record.contribution.filter(value => value === 0).length,
      cohortSum: sum(rows.map(row => row.contribution)), validShare};
  }

  function differenceData(state) {
    const a = D.comparisons[`${state.qwen}|${state.a}`]?.cka;
    const b = D.comparisons[`${state.qwen}|${state.b}`]?.cka;
    if (![a, b].every(record => record && Number.isFinite(record.global) && record.raw_denominator > 0 && record.contribution?.length === N && record.contribution.every(Number.isFinite))) {
      return {available: false, reason: 'Complete finite saved CKA allocations with nonzero centered-matrix norms are required for both checkpoints.', rows: []};
    }
    const differences = indices.map(i => a.contribution[i] - b.contribution[i]);
    const globalGap = a.global - b.global, allocatedGap = sum(differences);
    const tolerance = Number(a.validation?.tolerance ?? 1e-10) + Number(b.validation?.tolerance ?? 1e-10);
    if (!Number.isFinite(tolerance) || tolerance < 0 || Math.abs(allocatedGap - globalGap) > tolerance) {
      return {available: false, reason: 'Saved protein allocations do not reconstruct the full-panel global CKA difference within the recorded tolerances.', rows: []};
    }
    const positive = indices.filter(i => differences[i] > 0).sort((i, j) => differences[j] - differences[i] || ids[i].localeCompare(ids[j]));
    const negative = indices.filter(i => differences[i] < 0).sort((i, j) => differences[i] - differences[j] || ids[i].localeCompare(ids[j]));
    const cutsTie = order => order.length > 10 && differences[order[9]] === differences[order[10]];
    const chosen = [...positive.slice(0, 10), ...negative.slice(0, 10)].sort((i, j) => differences[j] - differences[i] || ids[i].localeCompare(ids[j]));
    return {available: true, globalGap, allocatedGap, tolerance, positiveCount: positive.length, negativeCount: negative.length,
      positiveTieCut: cutsTie(positive), negativeTieCut: cutsTie(negative),
      rows: chosen.map(i => ({i, id: ids[i], difference: differences[i], a: a.contribution[i], b: b.contribution[i]}))};
  }

  function differenceMarkup(state, data) {
    const q = qs.find(model => model.id === state.qwen);
    const modelKey = `Qwen ${q.size_b}B · ${fullName(state.a)} versus ${fullName(state.b)}`;
    const selector = `<div class="simple-controls"><label>Compare with protein model<select id="cka-b">${plms.map(model => `<option value="${esc(model.id)}" ${model.id === state.b ? 'selected' : ''}>${esc(shortName(model.id))}</option>`).join('')}</select></label></div>`;
    let note;
    if (!data.available) note = esc(data.reason);
    else if (!data.rows.length) note = `All ${N.toLocaleString()} saved protein allocations are equal for these two selections; there are no positive or negative bars to display.`;
    else note = `Showing ${data.rows.length} proteins: up to 10 largest positive and 10 most negative differences, ranked at this Qwen checkpoint. Exact ties are ordered by protein ID and the limit can split a tied group.${data.positiveTieCut ? ' The positive cutoff splits an exact tie.' : ''}${data.negativeTieCut ? ' The negative cutoff splits an exact tie.' : ''}`;
    const status = data.available ? `<p class="simple-explainer">Full-panel global CKA difference (${esc(shortName(state.a))} − ${esc(shortName(state.b))}): <strong>${signed(data.globalGap)}</strong>. Contributions from all ${N.toLocaleString()} proteins sum to this difference; the bars show only the stated extremes.</p>` : '';
    return `<section class="simple-section" id="cka-difference-section"><h2>Which proteins contribute to the CKA difference between PLMs?</h2>${selector}
      ${card(`Which proteins contribute most to the global CKA difference: ${shortName(state.a)} − ${shortName(state.b)}?`, modelKey,
        `<p class="simple-selection-note" id="cka-difference-cohort">${note}</p>${status}<div class="simple-actions"><button id="cka-difference-track" ${data.available && data.rows.length ? '' : 'disabled'}>Track these ${data.rows.length} proteins</button></div><div class="simple-heatmap-scroll cka-difference-scroll"><div id="cka-difference-plot" class="simple-plot" aria-label="Protein allocations of the global CKA difference between the selected PLMs"></div></div>`,
        `Right/solid means a higher allocation with ${shortName(state.a)}; left/hatched means higher with ${shortName(state.b)}. Zero means equal allocations. Values allocate the global CKA difference; they are not percentages, local CKA scores, or deletion effects. Click a bar to focus; Track replaces the tracked set.`)}
      <details class="simple-formula-details" id="cka-difference-formula"><summary>CKA difference formula and definitions</summary><div class="simple-equation" role="math">d<sub>i</sub> = c<sub>i,A</sub> − c<sub>i,B</sub><br>∑<sub>i=1</sub><sup>N</sup> d<sub>i</sub> = CKA<sub>A</sub> − CKA<sub>B</sub></div><p class="simple-math-key">i identifies one protein; N = ${N.toLocaleString()} is the complete fixed panel. c<sub>i,A</sub> and c<sub>i,B</sub> are its saved signed row allocations for the same Qwen checkpoint paired with A = ${esc(fullName(state.a))} and B = ${esc(fullName(state.b))}. d<sub>i</sub> is their difference, in dimensionless CKA score units. ∑ adds all N protein differences. These are raw contributions; no per-model percentage normalization is subtracted. The two saved CKA normalization denominators remain those of their respective full-panel comparisons.</p></details></section>`;
  }

  function functionCell(i) {
    const value = ann.function?.[i];
    if (!value) return 'Not recorded';
    return value.length > 130 ? `<details class="simple-function"><summary>${esc(value.slice(0, 130))}…</summary><p>${esc(value)}</p></details>` : esc(value);
  }

  function card(title, key, body, caption) {
    return `<article class="simple-card"><h3>${esc(title)}</h3><div class="simple-plot-key">${esc(key)}</div>${body}<p class="simple-read"><strong>How to read this:</strong> ${esc(caption)}</p></article>`;
  }

  function heatmapRows(state, data) {
    return heatmapCohort === 'tracked' ? [...new Set(state.selected || [])].filter(id => byId.has(id)).map(id => byId.get(id)) : data.rows.map(row => row.i);
  }

  function markup(state, data, difference) {
    const q = qs.find(model => model.id === state.qwen), modelKey = `Qwen ${q.size_b}B × ${fullName(state.a)}`;
    const tracked = new Set(state.selected || []), focus = state.focus;
    const share = value => data.validShare ? number(100 * value / data.record.global) + '%' : 'Undefined';
    const heatRows = heatmapRows(state, data);
    const cohortNote = heatmapCohort === 'table'
      ? `These ${heatRows.length} proteins are ranked at Qwen ${q.size_b}B with ${fullName(state.a)}. The same proteins and row order appear at every size; changing the ranking checkpoint or PLM updates this cohort.`
      : `All ${heatRows.length} tracked proteins appear in their tracked order. Their identities stay fixed when Qwen or the PLM changes. ${state.selectionOrigin || ''}`;
    const rows = data.rows.map(row => `<tr data-cka-protein-id="${esc(ids[row.i])}" data-cka-contribution="${row.contribution}" data-cka-family="${esc(row.familyId)}"><td>${row.rank}</td><td><button class="protein-link" data-cka-focus="${esc(ids[row.i])}" aria-pressed="${ids[row.i] === focus}" title="${esc(ann.protein_name?.[row.i] || ids[row.i])}">${esc(label(row.i))}</button>${tracked.has(ids[row.i]) ? '<span class="small-note"> · tracked</span>' : ''}</td><td>${share(row.contribution)}</td><td title="${esc(row.familyId)}">${esc(data.groups.get(row.familyId).label)}</td><td>${ann.length?.[row.i] ?? 'Not recorded'}</td><td>${functionCell(row.i)}</td></tr>`).join('');
    return `<style>
      #cka-view .cka-controls{display:grid;grid-template-columns:1fr 2fr 1fr;gap:15px}
      #cka-top-table th:first-child,#cka-top-table td:first-child{min-width:56px;white-space:nowrap;overflow-wrap:normal}
      #cka-top-table th:nth-child(2),#cka-top-table td:nth-child(2){min-width:175px;white-space:nowrap}
      #cka-top-table th:nth-child(4),#cka-top-table td:nth-child(4){min-width:210px;max-width:310px;white-space:normal}
      #cka-top-table th:last-child,#cka-top-table td:last-child{min-width:280px;max-width:420px;white-space:normal}
      #cka-top-table td:nth-child(3){font-variant-numeric:tabular-nums;white-space:nowrap}
      #cka-top-table [aria-pressed=true]{font-weight:700;text-decoration:underline}
      #cka-view .cka-family-scroll{max-height:640px;overflow:auto}#cka-view .cka-family-scroll>.simple-plot{min-width:760px}
      #cka-view .cka-protein-scroll{max-height:650px;overflow:auto}#cka-view .cka-protein-scroll>.simple-plot{min-width:850px}
      #cka-view .cka-difference-scroll{max-height:740px;overflow:auto}#cka-view .cka-difference-scroll>.simple-plot{min-width:760px}
      #cka-view .cka-heatmap-controls{display:flex;flex-wrap:wrap;gap:15px;margin-bottom:12px}#cka-view .cka-heatmap-controls label{min-width:220px}
      @media(max-width:620px){#cka-view .cka-controls{grid-template-columns:1fr}}
    </style>
    <section class="simple-section" id="cka-global-section"><h2>CKA · How similar are the models’ protein relationships?</h2><p class="simple-explainer">CKA compares centered cosine-similarity patterns across all ${N.toLocaleString()} proteins. Higher global CKA means the two models have more similar patterns on this panel.</p>
      ${card('How does global CKA compare across all models?', 'Columns: Qwen sizes · Rows: separate PLM checkpoints · Color and printed numbers: global CKA', '<div class="simple-heatmap-scroll"><div id="cka-global-heatmap" class="simple-plot" aria-label="Global CKA for every Qwen and protein model pair"></div></div>', 'Each cell uses the complete fixed protein panel. Read across a row to compare Qwen sizes, or down a column to compare PLMs. Click a cell to choose the pair below; tracked proteins stay fixed.')}
    </section>
    <section class="simple-section" id="cka-top-section"><h2>Which proteins contribute most to the global Qwen–PLM CKA score?</h2>
      <div class="simple-controls cka-controls"><label>Qwen checkpoint<select id="cka-qwen">${qs.map(model => `<option value="${esc(model.id)}" ${model.id === state.qwen ? 'selected' : ''}>Qwen ${model.size_b}B</option>`).join('')}</select></label><label>Protein model<select id="cka-a">${plms.map(model => `<option value="${esc(model.id)}" ${model.id === state.a ? 'selected' : ''}>${esc(shortName(model.id))}</option>`).join('')}</select></label><label>Show top positive contributors<select id="cka-top-n">${[10, 20, 50, 100].map(n => `<option value="${n}" ${n === data.requested ? 'selected' : ''}>${n}</option>`).join('')}</select></label></div>
      <p class="simple-plot-key"><strong>${esc(modelKey)}</strong> · Global CKA: <strong>${number(data.record.global)}</strong></p>
      <p class="simple-selection-note" id="cka-top-cohort"><strong>Top ${data.requested} + exact boundary ties = ${data.rows.length} positive contributors.</strong> ${data.cutoff === null ? 'No positive signed contributions are recorded.' : `Cutoff share: ${share(data.cutoff)}. Ranked by the stored full-precision contribution; equal values share a rank.`}</p>
      <div class="simple-actions"><button id="cka-top-track" ${data.rows.length ? '' : 'disabled'}>Track these ${data.rows.length} proteins</button><span class="simple-selection-note">${tracked.size} proteins currently tracked across reports.</span></div>
      <details id="cka-exact-values"><summary>Values and protein annotations</summary>${card('Which proteins account for the largest shares of global CKA?', modelKey,
        data.rows.length ? `<p class="simple-explainer"><strong>Share of global CKA</strong> is the percentage allocated to this protein. A share of 1% accounts for one percent of this model pair’s full-panel CKA score. Positive and negative shares sum to 100% across all ${N.toLocaleString()} proteins.</p><div class="simple-table-wrap"><table class="simple-table" id="cka-top-table"><thead><tr><th>Rank</th><th>Gene / protein name</th><th title="100 × this protein’s signed contribution / the full-panel global CKA score. This percentage can be negative or exceed 100% when contributions cancel.">Share of global CKA</th><th>Pfam family</th><th>Length (amino acids)</th><th>Function</th></tr></thead><tbody>${rows}</tbody></table></div>` : '<p class="empty">No positive contributors to display.</p>',
        'Each row gives a protein’s signed percentage of global CKA, not local CKA or a deletion effect. Parenthesized accessions distinguish repeated gene names; a missing gene name falls back to its accession. Click a gene to focus it without changing the tracked set.')}</details>
      <p class="simple-explainer">These ${data.rows.length} proteins account for ${share(data.cohortSum)} of global CKA. The remaining ${N - data.rows.length} account for ${share(data.record.global - data.cohortSum)}. Across the full panel: ${data.positive} positive, ${data.negative} negative and ${data.zero} zero contributions.</p>
      ${card('How does each protein’s share of global CKA change across Qwen sizes?', `${fullName(state.a)} · Columns: all Qwen sizes · Rows: the same proteins`,
        `<div class="cka-heatmap-controls"><label>Heatmap proteins<select id="cka-heatmap-cohort"><option value="table" ${heatmapCohort === 'table' ? 'selected' : ''}>Current top proteins</option><option value="tracked" ${heatmapCohort === 'tracked' ? 'selected' : ''}>Tracked proteins</option></select></label></div><p class="simple-explainer">A cell of 1% accounts for one percent of that column’s global CKA; −1% subtracts one percent. Across the full panel, signed shares sum to 100%.</p><p class="simple-selection-note" id="cka-heatmap-cohort-note">${esc(cohortNote)}</p><div class="simple-heatmap-scroll cka-protein-scroll" tabindex="0" aria-label="Scroll through every protein in the chosen heatmap cohort"><div id="cka-protein-heatmap" class="simple-plot" aria-label="The same proteins’ percentage shares of global CKA across all Qwen sizes"></div></div><p class="small-note">Scroll to see every row; no proteins in this heatmap cohort are omitted. Click a cell to focus its protein without changing the tracked set.</p>`,
        'Each column uses that model pair’s global CKA score. Positive shares run from pale blue to dark navy; negative shares are red, and zero is near-white. Darker shades mean a larger share in that direction; printed signs also show the direction. The color range is symmetric around zero for the displayed values. Changes show relative share, not necessarily growth in the raw contribution, because global CKA also changes. Shares are undefined when global CKA or the unregularized denominator is zero.')}
      ${card('How many top CKA contributors belong to each Pfam family?', `${modelKey} · Ranked table cohort: ${data.rows.length} proteins`, '<div class="simple-heatmap-scroll cka-family-scroll"><div id="cka-family-plot" class="simple-plot" aria-label="Pfam family counts among the exact CKA top-contributor cohort"></div></div>', `Each bar counts a Pfam group in the ranked table cohort, even when the heatmap shows tracked proteins. Counts add to ${data.rows.length}. Each protein uses its first sorted original Pfam accession; missing annotations remain a separate group.`)}
      <details class="simple-formula-details" id="cka-definition"><summary>CKA definition and signed contribution formula</summary><p><strong>K</strong> contains the cosine similarities between proteins in the chosen Qwen checkpoint; <strong>L</strong> contains them in the chosen PLM checkpoint. Both are N × N matrices with the same ordered proteins on their rows and columns.</p><p>H = I − 11ᵀ/N centers each matrix: K<sub>c</sub> = HKH and L<sub>c</sub> = HLH.</p><p class="simple-math-key"><strong>Centering symbols:</strong> N = ${N.toLocaleString()} is the full-panel protein count, not the number selected. H is the centering matrix; I is the N × N identity matrix (ones on the diagonal, zeros elsewhere); 1 is an N-entry column of ones. Superscript T means transpose, turning that column into a row, so 11ᵀ is an N × N matrix of ones. Adjacent matrices in HKH and HLH are multiplied. Subscript c means centered: each cell has its row and column means subtracted and the overall mean added back.</p><div class="simple-equation" role="math">c<sub>i</sub> = ∑<sub>j=1</sub><sup>N</sup> K<sub>c,ij</sub>L<sub>c,ij</sub> / (‖K<sub>c</sub>‖<sub>F</sub>‖L<sub>c</sub>‖<sub>F</sub> + ε)</div><p class="simple-math-key"><strong>Contribution symbols:</strong> i is the protein being attributed; j runs over every protein, including i. K<sub>c,ij</sub> and L<sub>c,ij</sub> are the centered similarities for that pair; their adjacent scalar values are multiplied. ∑ means add the terms for all indicated indices. ‖K<sub>c</sub>‖<sub>F</sub> and ‖L<sub>c</sub>‖<sub>F</sub> are Frobenius norms: the square root of the sum of squared entries in each entire matrix. ε is the fixed 10<sup>−12</sup> added to their norm product. The slash divides by this whole denominator. c<sub>i</sub> is protein i’s signed allocation of the global score; its letter c is distinct from the centered-matrix subscript c.</p><div class="simple-equation" role="math">CKA = ∑<sub>i=1</sub><sup>N</sup> c<sub>i</sub> &nbsp; · &nbsp; Share<sub>i</sub> = 100 × c<sub>i</sub> / CKA</div><p class="simple-math-key"><strong>Score and share:</strong> CKA is the full-panel alignment score, found by adding all N contributions. Share<sub>i</sub> is protein i’s percentage of that score; × means multiplication and 100 converts the fraction to percent. The centered dot between these two equations is only a separator.</p><p>This is the saved biased cosine-CKA estimator, including self pairs. Each c<sub>i</sub> can be positive or negative; cancellation means shares can exceed 100%. A share is undefined when global CKA is zero or the unregularized denominator is zero. The row allocation is signed accounting, not a standalone local CKA or the effect of deleting a protein. Profile and deletion analyses are available under All analyses.</p></details>
    </section>${differenceMarkup(state, difference)}`;
  }

  async function draw(id, traces, layout) {
    const element = $(id);
    if (!element?.getClientRects().length) return;
    element.style.height = `${layout.height}px`;
    await Plotly.newPlot(element, traces, {paper_bgcolor: '#fff', plot_bgcolor: '#fff', font: {family: 'Arial, sans-serif', size: 12, color: '#193243'}, showlegend: false, ...layout}, {responsive: false, displaylogo: false, scrollZoom: false});
  }

  async function change(patch) {
    api.ready = false;
    try {
      await window.DESCRIPTIVE.update(patch);
      await render(window.DESCRIPTIVE.getState());
    } catch (error) {
      api.error = error.message;
      api.ready = true;
      $('cka-view').innerHTML = `<p class="notice error">${esc(error.message)}</p>`;
    }
  }

  async function drawProteinHeatmap(state, data) {
    const rows = heatmapRows(state, data);
    if (!rows.length) {
      $('cka-protein-heatmap').innerHTML = '<p class="empty">No proteins in this cohort. Use “Track these proteins” above or choose “Current top proteins.”</p>';
      return;
    }
    const records = qs.map(q => D.comparisons[`${q.id}|${state.a}`]?.cka);
    const values = rows.map(i => records.map(record => Number.isFinite(record?.contribution?.[i]) && Number.isFinite(record?.global) && record.global !== 0 && record.raw_denominator > 0 ? 100 * record.contribution[i] / record.global : null));
    const finiteValues = values.flat().filter(Number.isFinite);
    const extent = finiteValues.reduce((largest, value) => Math.max(largest, Math.abs(value)), 0) || 1;
    const customdata = rows.map((i, row) => qs.map((q, column) => [ids[i], q.id, values[row][column], records[column]?.global ?? null]));
    const focusedRow = rows.findIndex(i => ids[i] === state.focus);
    await draw('cka-protein-heatmap', [{type: 'heatmap', x: qs.map(q => q.id), y: rows.map(i => ids[i]), z: values,
      customdata, text: values.map(row => row.map(value => Number.isFinite(value) ? signed(value) + '%' : 'Undefined')), texttemplate: '%{text}', textfont: {size: 11, color: 'auto'},
      hovertext: rows.map((i, row) => qs.map((q, column) => `${esc(label(i))}<br>Protein ID: ${esc(ids[i])}<br>Qwen ${q.size_b}B · ${esc(fullName(state.a))}<br>Share of global CKA: ${Number.isFinite(values[row][column]) ? signed(values[row][column]) + '%' : 'Undefined'}<br>Global CKA: ${number(records[column]?.global)}`)),
      hovertemplate: '%{hovertext}<extra></extra>', hoverongaps: false, zmin: -extent, zmax: extent, zmid: 0,
      colorscale: [[0, '#67001f'], [.1, '#b2182b'], [.2, '#d6604d'], [.3, '#f4a582'], [.4, '#fddbc7'], [.5, '#f7f7f7'], [.6, '#d1e5f0'], [.7, '#92c5de'], [.8, '#4393c3'], [.9, '#2166ac'], [1, '#053061']], xgap: 1, ygap: 1,
      colorbar: {title: {text: 'Share of<br>global CKA (%)'}, tickmode: 'array', tickvals: [-extent, 0, extent], ticktext: [-extent, 0, extent].map(number), ticks: 'outside', ticklen: 5, tickcolor: '#193243', thickness: 13, lenmode: 'pixels', len: 230, yanchor: 'top', y: 1}}],
      {height: Math.max(350, rows.length * 26 + 135),
        xaxis: {title: {text: 'Qwen model size'}, type: 'category', tickvals: qs.map(q => q.id), ticktext: qs.map(q => `${q.size_b}B`), side: 'top', automargin: true},
        yaxis: {title: {text: 'Gene / protein name'}, type: 'category', tickvals: rows.map(i => ids[i]), ticktext: rows.map(i => esc(label(i))), categoryorder: 'array', categoryarray: rows.map(i => ids[i]), autorange: 'reversed', automargin: true},
        margin: {l: 200, r: 150, t: 85, b: 35},
        shapes: focusedRow < 0 ? [] : [{type: 'rect', xref: 'paper', x0: 0, x1: 1, yref: 'y', y0: focusedRow - .48, y1: focusedRow + .48, line: {color: '#193243', width: 2}, fillcolor: 'rgba(0,0,0,0)'}]});
    $('cka-protein-heatmap').on?.('plotly_click', event => {
      const id = event.points?.[0]?.customdata?.[0];
      if (typeof id === 'string' && byId.has(id)) { api.ready = false; setTimeout(() => change({focus: id}), 0); }
    });
  }

  async function drawDifference(state, data) {
    if (!data.available || !data.rows.length) {
      $('cka-difference-plot').innerHTML = `<p class="empty">${esc(data.available ? state.a === state.b ? `${shortName(state.a)} is selected twice. Its contribution differences are all zero.` : 'All saved protein allocations are equal; no nonzero differences to show.' : data.reason)}</p>`;
      return;
    }
    const extent = Math.max(...data.rows.map(row => Math.abs(row.difference))) * 1.12;
    const traces = [true, false].map(positive => {
      const rows = data.rows.filter(row => positive ? row.difference > 0 : row.difference < 0);
      const modelId = positive ? state.a : state.b;
      return {type: 'bar', orientation: 'h', name: `${positive ? 'Positive' : 'Negative'}: ${shortName(modelId)} higher allocation`,
        x: rows.map(row => row.difference), y: rows.map(row => row.id), customdata: rows.map(row => row.id),
        hovertext: rows.map(row => `${esc(label(row.i))}<br>Protein ID: ${esc(row.id)}<br>Qwen ${qs.find(q => q.id === state.qwen).size_b}B<br>${esc(fullName(state.a))} contribution: ${signed(row.a)}<br>${esc(fullName(state.b))} contribution: ${signed(row.b)}<br>${esc(shortName(state.a))} − ${esc(shortName(state.b))} contribution: ${signed(row.difference)}`),
        hovertemplate: '%{hovertext}<extra></extra>', textposition: 'none',
        marker: {color: color(modelId), line: {color: color(modelId), width: 1}, pattern: {shape: positive ? '' : '/', solidity: .25}}};
    });
    await draw('cka-difference-plot', traces, {height: Math.max(360, data.rows.length * 27 + 155), showlegend: true,
      legend: {orientation: 'h', x: 0, y: 1.15, font: {size: 11}},
      xaxis: {title: {text: 'Contribution to global CKA difference<br>(dimensionless)'}, range: [-extent, extent], tickformat: '.3~r', exponentformat: 'none', showexponent: 'none', nticks: 5, zeroline: true, zerolinewidth: 1.5, zerolinecolor: '#193243'},
      yaxis: {title: {text: 'Gene / protein name'}, type: 'category', categoryorder: 'array', categoryarray: data.rows.map(row => row.id), autorange: 'reversed', tickvals: data.rows.map(row => row.id), ticktext: data.rows.map(row => esc(label(row.i))), automargin: true},
      margin: {l: 210, r: 30, t: 75, b: 75}});
    $('cka-difference-plot').on?.('plotly_click', event => {
      const id = event.points?.[0]?.customdata;
      if (typeof id === 'string' && byId.has(id)) { api.ready = false; setTimeout(() => change({focus: id}), 0); }
    });
  }

  async function renderNow(state) {
    if (!qs.some(q => q.id === state.qwen) || !plms.some(p => p.id === state.a)) throw Error('Choose a recorded Qwen checkpoint and PLM.');
    const data = topData(state), difference = differenceData(state), host = $('cka-view');
    if (!host) throw Error('CKA report host #cka-view is missing.');
    const open = new Set([...host.querySelectorAll('details[open][id]')].map(node => node.id));
    for (const node of host.querySelectorAll('.simple-plot')) if (node.data) Plotly.purge(node);
    host.innerHTML = markup(state, data, difference);
    for (const id of open) if ($(id)) $(id).open = true;
    $('cka-qwen').onchange = event => change({qwen: event.target.value});
    $('cka-a').onchange = event => change({a: event.target.value});
    $('cka-b').onchange = event => change({b: event.target.value});
    $('cka-difference-track').onclick = () => { if (difference.available && difference.rows.length) change({selected: difference.rows.map(row => row.id), focus: difference.rows[0].id, selectionOrigin: `Up to 10 positive and 10 negative raw CKA allocation differences for ${state.qwen}: ${fullName(state.a)} minus ${fullName(state.b)}. Exact cutoff ties use protein ID; these raw IDs stay fixed across views.`}); };
    $('cka-top-n').onchange = event => change({topN: Number(event.target.value)});
    $('cka-heatmap-cohort').onchange = event => { heatmapCohort = event.target.value; render(window.DESCRIPTIVE.getState()); };
    $('cka-top-track').onclick = () => change({selected: data.rows.map(row => ids[row.i]), focus: ids[data.rows[0].i], selectionOrigin: `Top ${data.requested} positive CKA contributions plus exact boundary ties: ${data.rows.length} proteins for ${state.qwen} and ${fullName(state.a)}. These raw IDs remain fixed across reports and models.`});
    for (const button of host.querySelectorAll('[data-cka-focus]')) button.onclick = () => change({focus: button.dataset.ckaFocus});
    const values = plms.map(p => qs.map(q => D.comparisons[`${q.id}|${p.id}`].cka.global));
    const qi = qs.findIndex(q => q.id === state.qwen), pi = plms.findIndex(p => p.id === state.a);
    await draw('cka-global-heatmap', [{type: 'heatmap', x: qs.map(q => q.id), y: plms.map(p => p.id), z: values, text: values.map(row => row.map(number)), texttemplate: '%{text}', textfont: {size: 13, color: 'auto'}, customdata: plms.map(p => qs.map(q => `${q.id}|${p.id}`)), hovertext: plms.map(p => qs.map(q => `${esc(q.label)}<br>${esc(p.label)}`)), hovertemplate: '%{hovertext}<br>Global CKA: %{z:.3g}<extra></extra>', zmin: 0, zmax: 1, colorscale: 'Viridis', colorbar: {title: {text: 'Global CKA<br>(dimensionless)'}, thickness: 13, tickformat: '.3~g'}, xgap: 2, ygap: 2}],
      {height: 520, xaxis: {title: {text: 'Qwen model size'}, type: 'category', tickvals: qs.map(q => q.id), ticktext: qs.map(q => `${q.size_b}B`), automargin: true}, yaxis: {title: {text: 'Protein model checkpoint'}, type: 'category', tickvals: plms.map(p => p.id), ticktext: plms.map(p => shortName(p.id)), autorange: 'reversed', automargin: true}, margin: {l: 160, r: 145, t: 16, b: 65}, shapes: [{type: 'rect', x0: qi - .48, x1: qi + .48, y0: pi - .48, y1: pi + .48, line: {color: color(state.a), width: 3}, fillcolor: 'rgba(0,0,0,0)'}]});
    $('cka-global-heatmap').on?.('plotly_click', event => {
      const key = event.points?.[0]?.customdata;
      if (typeof key !== 'string') return;
      const [qwen, a] = key.split('|');
      api.ready = false;
      setTimeout(() => change({qwen, a}), 0);
    });
    await drawProteinHeatmap(state, data);
    await drawDifference(state, difference);
    if (!data.families.length) { $('cka-family-plot').innerHTML = '<p class="empty">No positive-contributor cohort to count.</p>'; return; }
    const wrap = text => {
      const lines = [''];
      for (const word of text.split(/\s+/)) { const j = lines.length - 1; if (lines[j] && lines[j].length + word.length + 1 > 36) lines.push(word); else lines[j] += (lines[j] ? ' ' : '') + word; }
      return lines.map(esc).join('<br>');
    };
    const ticks = data.families.map(group => wrap(group.label));
    const rowHeight = Math.max(29, ...ticks.map(text => text.split('<br>').length * 13 + 6));
    const largest = Math.max(...data.families.map(group => group.n));
    await draw('cka-family-plot', [{type: 'bar', orientation: 'h', x: data.families.map(group => group.n), y: data.families.map(group => group.id), text: data.families.map(group => `${esc(group.label)}<br>${group.n} of ${data.rows.length} top contributors<br>Recorded Pfam: ${esc(group.id)}`), texttemplate: '%{x}', textposition: 'outside', cliponaxis: false, hovertemplate: '%{text}<extra></extra>', marker: {color: data.families.map(group => group.id === 'unknown' ? '#667985' : color(state.a))}}],
      {height: Math.max(280, data.families.length * rowHeight + 90), xaxis: {title: {text: 'Number of top contributors'}, side: 'top', range: [0, largest * 1.16 + .2], dtick: Math.max(1, Math.ceil(largest / 5)), tickformat: 'd'}, yaxis: {title: {text: 'Pfam family'}, type: 'category', categoryorder: 'array', categoryarray: data.families.map(group => group.id), autorange: 'reversed', tickvals: data.families.map(group => group.id), ticktext: ticks}, margin: {l: 270, r: 35, t: 68, b: 20}});
  }

  function render(sharedState) {
    latest = {...sharedState, selected: [...(sharedState.selected || [])]};
    if (latest.focus && !byId.has(latest.focus)) throw Error('Focused protein is outside the fixed panel.');
    if (active) { pending = true; return task; }
    active = true; api.ready = false; api.error = null;
    task = (async () => {
      try { do { pending = false; await renderNow(latest); } while (pending); }
      catch (error) { api.error = error.message; if ($('cka-view')) $('cka-view').innerHTML = `<p class="notice error">${esc(error.message)}</p>`; console.error(error); }
      finally { active = false; api.ready = true; }
    })();
    return task;
  }
})();
