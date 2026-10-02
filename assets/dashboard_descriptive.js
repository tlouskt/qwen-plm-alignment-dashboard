/* Descriptive kNN views of the existing fixed-panel counts. No refits or inference. */
(() => {
  'use strict';
  const D=window.ALIGNMENT_DATA, $=id=>document.getElementById(id);
  const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const ids=D.panel.protein_ids,N=ids.length,ix=ids.map((_,i)=>i),idIndex=new Map(ids.map((id,i)=>[id,i]));
  const genes=ix.map(i=>String(D.panel.annotations.gene_name?.[i]??'').trim()),geneCounts=new Map();
  genes.forEach(g=>{if(g)geneCounts.set(g,(geneCounts.get(g)||0)+1);});
  const label=i=>genes[i]?(geneCounts.get(genes[i])>1?`${genes[i]} (${ids[i]})`:genes[i]):ids[i];
  const plms=D.plms,byPLM=Object.fromEntries(plms.map(p=>[p.id,p]));
  const qs=[...D.qwen].sort((a,b)=>a.size_b-b.size_b),first=qs[0],last=qs[qs.length-1];
  const find=(pattern,fallback)=>plms.find(p=>pattern.test(p.id))?.id||fallback;
  const state={section:'proteins',qwen:last.id,topN:20,proteinOrder:'high',proteinPage:0,a:find(/prot.?t5/i,plms[0].id),b:find(/esmc_300/i,plms[1].id),k:Number(D.settings.knn_ks[0]),search:'',pattern:'all',sort:'change',page:0,pageSize:12,focus:0,selected:new Set(),selectionOrigin:''};
  const cache=new Map(),plots=new Set(),A='#0072b2',B='#d55e00';
  const sum=a=>a.reduce((s,x)=>s+x,0),number=(x,n=3)=>Number.isFinite(x)?(n<=2?x.toLocaleString('en-US',{maximumFractionDigits:n,minimumFractionDigits:n}):Math.abs(x)>0&&Math.abs(x)<1e-6?x.toExponential(2):x.toLocaleString('en-US',{maximumSignificantDigits:3})):'—';
  const signed=(x,n=3)=>(x>0?'+':'')+number(x,n),model=id=>byPLM[id]?.label||id;
  const patternName=p=>({gain:'Gain',loss:'Loss',retained:'Same endpoints'})[p];
  function series(p) {
    const key=`${p}|${state.k}`;if(cache.has(key))return cache.get(key);
    const records=qs.map(q=>D.comparisons[`${q.id}|${p}`]?.knn?.[String(state.k)]);
    if(records.some(r=>!r||r.count?.length!==N||r.count.some(c=>!Number.isInteger(c)||c<0||c>state.k)))throw Error('A complete fixed-panel count series is required.');
    const rows=ix.map(i=>{const counts=records.map(r=>r.count[i]),deltaCount=counts.at(-1)-counts[0],steps=counts.slice(1).map((c,j)=>c-counts[j]);return {i,counts,values:counts.map(c=>c/state.k),deltaCount,delta:deltaCount/state.k,contribution:deltaCount/(N*state.k),pattern:deltaCount>0?'gain':deltaCount<0?'loss':'retained',steadyGain:deltaCount>0&&steps.every(x=>x>=0),steadyLoss:deltaCount<0&&steps.every(x=>x<=0),flat:counts.every(c=>c===counts[0])};});
    const start=sum(records[0].count)/(N*state.k),end=sum(records.at(-1).count)/(N*state.k);
    const result={rows,start,end,net:sum(rows.map(r=>r.deltaCount))/(N*state.k),positive:sum(rows.map(r=>Math.max(0,r.deltaCount)))/(N*state.k),negative:sum(rows.map(r=>Math.min(0,r.deltaCount)))/(N*state.k),archiveDelta:records.at(-1).global-records[0].global};
    if(Math.abs(result.net-result.archiveDelta)>5e-8)throw Error('Endpoint accounting differs from the recorded global change.');
    cache.set(key,result);return result;
  }
  function filtered() {
    const a=series(state.a).rows,b=series(state.b).rows,search=state.search.trim().toLowerCase();
    const list=ix.filter(i=>(!search||[label(i),ids[i],D.panel.annotations.protein_name?.[i]].some(x=>String(x??'').toLowerCase().includes(search)))&&(state.pattern==='all'||a[i].pattern===state.pattern));
    return list.sort((i,j)=>{let d=0;if(state.sort==='gene')d=label(i).localeCompare(label(j));else if(state.sort==='loss')d=a[i].deltaCount-a[j].deltaCount;else if(state.sort==='contrast')d=Math.abs(a[j].deltaCount-b[j].deltaCount)-Math.abs(a[i].deltaCount-b[i].deltaCount);else d=a[j].deltaCount-a[i].deltaCount;return d||ids[i].localeCompare(ids[j]);});
  }
  function chooseExamples() {
    const a=series(state.a).rows,ordered=[...ix].sort((i,j)=>a[j].deltaCount-a[i].deltaCount||ids[i].localeCompare(ids[j]));
    const gains=ordered.filter(i=>a[i].deltaCount>0).slice(0,6),losses=[...ix].filter(i=>a[i].deltaCount<0).sort((i,j)=>a[i].deltaCount-a[j].deltaCount||ids[i].localeCompare(ids[j])).slice(0,6);
    state.selected=new Set([...gains,...losses]);state.focus=gains[0]??0;state.selectionOrigin=`Six largest gains and six largest losses with ${model(state.a)}, k=${state.k}; tied changes ordered by protein ID. IDs stay fixed when models or k change.`;
  }
  function chooseDifferenceExamples() {
    const a=series(state.a).rows,b=series(state.b).rows;
    const magnitude=ix.map(i=>a[i].counts.map((v,j)=>Math.abs(v-b[i].counts[j])));
    const order=[...ix].sort((i,j)=>Math.max(...magnitude[j])-Math.max(...magnitude[i])||sum(magnitude[j])-sum(magnitude[i])||ids[i].localeCompare(ids[j]));
    state.selected=new Set(order.slice(0,40));state.focus=order[0];
    state.selectionOrigin=`Illustrative 40 proteins with the largest absolute ${model(state.a)} minus ${model(state.b)} difference across sizes (then mean absolute difference, then ID), k=${state.k}. This selection does not represent the whole panel.`;
  }
  function accounting(p) {const s=series(p),selected=[...state.selected];const pos=sum(selected.map(i=>Math.max(0,s.rows[i].deltaCount)))/(N*state.k),neg=sum(selected.map(i=>Math.min(0,s.rows[i].deltaCount)))/(N*state.k),net=pos+neg;return {...s,selectedPositive:pos,selectedNegative:neg,selectedNet:net,remaining:s.net-net,selectedShare:Math.abs(s.net)>1e-15?net/s.net:null};}
  function csv(rows) {if(!rows.length)return '';const keys=Object.keys(rows[0]),cell=v=>{let s=String(v??'');if(typeof v==='string'&&/^[=+@\-]/.test(s))s="'"+s;return /[",\r\n]/.test(s)?'"'+s.replace(/"/g,'""')+'"':s;};return keys.join(',')+'\n'+rows.map(r=>keys.map(k=>cell(r[k])).join(',')).join('\n')+'\n';}
  function download(name,rows) {const url=URL.createObjectURL(new Blob([csv(rows)],{type:'text/csv;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  function exportRows() {return ix.flatMap(i=>[...new Set([state.a,state.b])].flatMap(p=>{const r=series(p).rows[i];return qs.map((q,j)=>({panel_id:D.panel.id,protein_id:ids[i],gene_name:genes[i],protein_label:label(i),pfam_family:familyLabel(i),sequence_length:D.panel.annotations.length?.[i]??'',function:D.panel.annotations.function?.[i]??'',plm_checkpoint:p,qwen_checkpoint:q.id,qwen_size_b:q.size_b,k:state.k,panel_n:N,shared_neighbors:r.counts[j],local_alignment:r.values[j],endpoint_local_change:r.delta,contribution_to_global_change:r.contribution,endpoint_pattern:r.pattern,selected:state.selected.has(i)}));}));}
  function card(title,body) {const id=body.match(/id="([^"]+)" class="simple-plot"/)?.[1];if(!id)return `<article class="simple-card"><h3>${esc(title)}</h3>${body}</article>`;const guide=readingGuide(id,title);return `<article class="simple-card" data-reading-plot="${id}"><h3>${esc(guide.question)}</h3><div class="simple-plot-key" aria-label="Plot legend">${esc(guide.key)}</div>${body}<p class="simple-read"><strong>How to read this:</strong> ${esc(guide.caption)}</p></article>`;}
  function table(headers,rows,id='') {return `<div class="simple-table-wrap"><table class="simple-table"${id?` id="${id}"`:''}><thead><tr>${headers.map(h=>`<th>${esc(h)}</th>`).join('')}</tr></thead><tbody>${rows.map(r=>`<tr>${r.map(c=>`<td>${c}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;}
  async function plot(id,traces,layout) {
    if(!$(id)?.getClientRects().length)return;
    layout.xaxis={automargin:true,tickformat:'.3~g',...layout.xaxis};layout.yaxis={automargin:true,tickformat:'.3~g',...layout.yaxis};
    for(const trace of traces)if(trace.colorbar)trace.colorbar={tickformat:'.3~g',...trace.colorbar};
    if(layout.height)$(id).style.height=`${layout.height}px`;
    plots.add(id);await Plotly.newPlot($(id),traces,{paper_bgcolor:'#fff',plot_bgcolor:'#fff',font:{family:'Arial, sans-serif',size:12,color:'#193243'},margin:{l:62,r:22,t:30,b:65},showlegend:true,legend:{orientation:'h',y:1.15},...layout},{responsive:false,displaylogo:false,scrollZoom:false});
  }
  function familyLabel(i) {const a=D.panel.annotations,g=a.group?.[i],entry=D.panel.group_names?.[g];return entry?.description||entry?.name||a.group_description?.[i]||a.group_name?.[i]||'No recorded Pfam family';}
  function functionCell(i) {const value=D.panel.annotations.function?.[i];if(!value)return 'Not recorded';return value.length>130?`<details class="simple-function"><summary>${esc(value.slice(0,130))}…</summary><p>${esc(value)}</p></details>`:esc(value);}
  function tableContent() {
    const a=series(state.a),b=series(state.b),rows=filtered(),pages=Math.max(1,Math.ceil(rows.length/state.pageSize));state.page=Math.min(state.page,pages-1);
    const shown=rows.slice(state.page*state.pageSize,(state.page+1)*state.pageSize);
    return `<div class="table-toolbar"><span>${rows.length.toLocaleString()} matching proteins · ${N.toLocaleString()} in the fixed panel</span><div class="pagination"><button id="simple-prev" ${state.page?'':'disabled'}>Previous</button><span>${state.page+1} / ${pages}</span><button id="simple-next" ${state.page<pages-1?'':'disabled'}>Next</button></div></div>`+table(['Select','Gene / protein',...qs.map(q=>`${q.size_b}B`),'Local change','Endpoint pattern','Across all sizes','Pfam family','Length (amino acids)','Function'],shown.map(i=>{const r=a.rows[i];return [`<input type="checkbox" data-simple-select="${i}" aria-label="Select ${esc(label(i))}" ${state.selected.has(i)?'checked':''}>`,`<button class="protein-link" data-simple-focus="${i}">${esc(label(i))}</button>`,...r.values.map(v=>number(v,2)),`<span class="${r.delta>0?'simple-gain':r.delta<0?'simple-loss':''}">${signed(r.delta,2)}</span>`,esc(patternName(r.pattern)),esc(r.flat?'Unchanged throughout':r.steadyGain?'Never decreases':r.steadyLoss?'Never increases':'Varies across sizes'),esc(familyLabel(i)),D.panel.annotations.length?.[i]??'Not recorded',functionCell(i)];}),'simple-proteins')+`<p class="simple-explainer">Values are shared neighbors ÷ ${state.k}; 0.30 means ${number(.3*state.k,0)} of ${state.k} neighbors agree. Change is ${first.size_b}B → ${last.size_b}B. “Same endpoints” can include changes in between. Use Model differences to inspect the same protein with both PLMs.</p>`;
  }
  const shortNames={prot_t5_xl_uniref50:'ProtT5 XL',esmc_300m:'ESMC 300M',esmc_600m:'ESMC 600M',esm3_sm_open_v1:'ESM3 small',esm2_t6_8m:'ESM2 8M',esm2_t12_35m:'ESM2 35M',esm2_t30_150m:'ESM2 150M',esm2_t33_650m:'ESM2 650M',esmfold_v1:'ESMFold v1'};
  const shortModel=id=>shortNames[id]||model(id);
  const modelColors=['#0072b2','#d55e00','#009e73','#cc79a7','#7057a3','#7a6a31','#4b7f86','#b34b59','#696969'];
  const colorOf=p=>modelColors[plms.findIndex(m=>m.id===p)%modelColors.length];
  const modelDash=p=>['solid','dash','dot','dashdot'][plms.findIndex(m=>m.id===p)%4];
  const modelSymbol=p=>['circle','square','diamond','cross','x','triangle-up','triangle-down','pentagon','hexagram'][plms.findIndex(m=>m.id===p)%9];
  const checkpointIndex=()=>qs.findIndex(q=>q.id===state.qwen);
  const localAt=p=>series(p).rows.map(r=>r.values[checkpointIndex()]);
  const pairCounts=()=>ix.map(i=>series(state.a).rows[i].counts[checkpointIndex()]-series(state.b).rows[i].counts[checkpointIndex()]);
  function pairLegend(){return [{type:'scatter',mode:'markers',x:[null],y:[null],name:`Positive: ${shortModel(state.a)} higher`,marker:{color:A,size:8},hoverinfo:'skip'},{type:'scatter',mode:'markers',x:[null],y:[null],name:`Negative: ${shortModel(state.b)} higher`,marker:{color:B,size:8},hoverinfo:'skip'}];}
  function readingGuide(id,title){
    const chosen=`Qwen ${qs[checkpointIndex()].size_b}B · ${shortModel(state.a)} · k=${state.k}`;
    const pair=`${shortModel(state.a)} versus ${shortModel(state.b)}`;
    const guides={
      'guided-top-family-plot':['How many top-alignment proteins belong to each Pfam family?',`Each bar counts members of one Pfam family in the exact top-protein cohort. The numbers at bar ends give the counts.`,chosen],
      'guided-global-heatmap':['How does global kNN alignment compare across Qwen and PLM checkpoints?',`Each cell is the mean local alignment across all ${N.toLocaleString()} proteins for one Qwen–PLM pair. Printed values and the color scale show the same score, from 0 to 1.`,plms.map(p=>shortModel(p.id)).join(' · ')],
      'guided-global-lines':['How does global kNN alignment change across Qwen sizes?',`Each named line follows one PLM across the five Qwen checkpoints. Markers show measured global scores; line patterns and marker shapes also distinguish models.`,plms.map(p=>shortModel(p.id)).join(' · ')],
      'guided-ecdf':['How are local kNN scores distributed across proteins?',`Each model’s curve gives the percentage of proteins at or below a local score. A curve further right has higher scores across more of the panel.`,plms.map(p=>shortModel(p.id)).join(' · ')],
      'guided-gap-distribution':['What percentage of proteins has each signed local kNN difference?',`Bars count the percentage of proteins at each local-score difference. Positive values favor ${shortModel(state.a)}; negative values favor ${shortModel(state.b)}; zero means equal scores.`,pair],
      'guided-gap-concentration':['How concentrated are the absolute local kNN differences?',`The solid curve adds absolute ${shortModel(state.a)} − ${shortModel(state.b)} differences, largest first. The dotted diagonal represents equal contributions from every protein. A steep initial rise means a small subset accounts for much of the difference.`,`Solid: cumulative absolute difference · Dotted: equal contributions · ${pair}`],
      'guided-gap-proteins':['Which proteins have the largest local kNN differences in either direction?',`Each bar is one protein’s local ${shortModel(state.a)} − ${shortModel(state.b)} difference. Right of zero favors ${shortModel(state.a)}; left favors ${shortModel(state.b)}. Bar position and sign carry the meaning as well as color.`,pair],
      'simple-difference-heatmap':['How do selected proteins’ signed local kNN differences change across Qwen sizes?',`Each cell is local alignment with ${shortModel(state.a)} minus local alignment with ${shortModel(state.b)}, for the same protein at that Qwen size. Blue / positive favors ${shortModel(state.a)}; red / negative favors ${shortModel(state.b)}; the light midpoint at zero means equal alignment.`,pair],
      'simple-heat-0':[`How do selected proteins’ local kNN scores with ${shortModel(state.a)} change across Qwen sizes?`,`Each cell is that protein’s local shared-neighbor fraction between Qwen and ${shortModel(state.a)}. Numbers and shading both show scores from 0 to 1.`,`Protein model: ${model(state.a)}`],
      'simple-heat-1':[`How do selected proteins’ local kNN scores with ${shortModel(state.b)} change across Qwen sizes?`,`Each cell is that protein’s local shared-neighbor fraction between Qwen and ${shortModel(state.b)}. Numbers and shading both show scores from 0 to 1.`,`Protein model: ${model(state.b)}`],
      'simple-focus-plot':[`How does local kNN alignment for ${label(state.focus)} change across Qwen sizes?`,`The two named lines follow the same protein. Higher points mean more shared neighbors with Qwen. Marker shapes and line patterns distinguish the models.`,pair],
      'simple-trajectories':['How do tracked proteins’ local kNN scores change across Qwen sizes?',`Each line follows one tracked protein. Solid lines have an endpoint gain, dashed lines a loss, dotted lines equal endpoints. The thick line with diamond markers is ${label(state.focus)}.`,`PLM: ${model(state.a)} · Solid: gain · Dashed: loss · Dotted: same endpoints`],
      'simple-plm-change':[`Do local kNN changes differ between ${shortModel(state.a)} and ${shortModel(state.b)}?`,`Each point is one of the ${N.toLocaleString()} reference proteins, showing its endpoint change with ${shortModel(state.a)} and ${shortModel(state.b)}. The dotted diagonal marks equal changes. Above it means a larger gain with ${shortModel(state.a)}; below it a larger gain with ${shortModel(state.b)}. Diamonds mark tracked proteins; other up/down triangles mark ${shortModel(state.a)} gains/losses, circles unchanged endpoints. Zero axes separate gains from losses.`,pair+' · Dotted line: equal endpoint changes'],
      'simple-account-plot':['How do local kNN gains and losses sum to the global change?',`Each model has bars for positive changes, negative changes, and their net sum. Heights are percentage points of the global kNN score. Patterns and named legends distinguish the models.`,pair],
      'guided-biology-family':['What is the mean signed local kNN difference in each of the 12 largest Pfam groups?',`Each point is a family’s mean local alignment difference; its stem connects to zero (equal alignment). Left favors ${shortModel(state.b)}; right favors ${shortModel(state.a)}. Labels show the difference and family size n.`,`Qwen ${qs[checkpointIndex()].size_b}B · k=${state.k} · ${pair}`],
      'guided-biology-taxon':['How does the mean signed local kNN difference vary by collection group?',`Each bar is the mean local ${shortModel(state.a)} − ${shortModel(state.b)} difference for a named collection group. Zero means equal mean local scores.`,pair],
      'guided-redundancy-gap':['Does the global kNN gap remain after keeping one protein per name?',`Each line shows the global ${shortModel(state.a)} − ${shortModel(state.b)} gap under one named control. Positive favors ${shortModel(state.a)} and negative favors ${shortModel(state.b)}. Distinct line patterns mark the full panel, one-per-name, and matched random-deletion results.`,pair],
      'guided-redundancy-paired':['How do global kNN scores compare for name selection and matched deletion?',`Each model’s line shows one-per-name minus matched-deletion alignment. Positive means the name-selected panel scored higher. Whiskers span the 20 recorded draws; they are ranges, not confidence intervals.`,pair]
    };
    const v=guides[id];return v?{question:v[0],caption:v[1],key:v[2]}:{question:title.endsWith('?')?title:title+'?',caption:'Hover over a mark for its protein, model and measured value.',key:''};
  }

  function globalMarkup(){return `<section class="simple-section" id="guided-global"><h2>How does global kNN alignment compare across protein models?</h2>${card('All PLM checkpoints across Qwen sizes','<div class="simple-heatmap-scroll"><div id="guided-global-heatmap" class="simple-plot"></div></div><p class="simple-math-key">Global kNN alignment α: mean local alignment across all '+N.toLocaleString()+' proteins. Click a cell to inspect that Qwen checkpoint and PLM.</p>')}<details id="guided-global-lines-details"><summary>Global alignment trajectories for all models</summary>${card('Global alignment across Qwen sizes','<div id="guided-global-lines" class="simple-plot"></div>')}</details></section>`;}
  function distributionMarkup(){const counts=pairCounts(),nPos=counts.filter(v=>v>0).length,nNeg=counts.filter(v=>v<0).length;return `<section class="simple-section" id="guided-distribution"><h2>How are local kNN differences distributed across proteins?</h2><p class="simple-comparison-key"><span><strong>Qwen:</strong> ${qs[checkpointIndex()].size_b}B · k = ${state.k}</span><span><strong>${esc(model(state.a))}</strong></span><span>versus <strong>${esc(model(state.b))}</strong></span></p>${card('What percentage of proteins has each signed local kNN difference?','<div id="guided-gap-distribution" class="simple-plot"></div>')}<p class="simple-explainer">${nPos.toLocaleString()} proteins favor ${esc(shortModel(state.a))}; ${nNeg.toLocaleString()} favor ${esc(shortModel(state.b))}; ${(N-nPos-nNeg).toLocaleString()} are tied.</p><div class="simple-equation" role="math" aria-label="The paired local difference is A minus B">D<sub>i,k</sub>(s) = ℓ<sub>i,k</sub><sup>A</sup>(s) − ℓ<sub>i,k</sub><sup>B</sup>(s)</div><p class="simple-math-key"><strong>D</strong> is one protein’s local alignment difference, A minus B. <strong>ℓ</strong> is its local alignment score (shared neighbors divided by k). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters.<br><strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. Positive D favors A; negative D favors B; zero means equal local scores.</p><details id="guided-ecdf-details"><summary>Compare local-score distributions across all models</summary>${card('Local alignment distributions · all PLMs','<div id="guided-ecdf" class="simple-plot"></div>')}</details><details id="guided-concentration-details"><summary>See how much of the difference comes from a subset</summary>${card('How concentrated are the differences?','<div id="guided-gap-concentration" class="simple-plot"></div>')}<div class="simple-equation" role="math" aria-label="Share of absolute difference accounted for by the m largest absolute differences">C(m) = ∑<sub>j=1</sub><sup>m</sup> |D<sub>(j),k</sub>| / ∑<sub>i=1</sub><sup>N</sup> |D<sub>i,k</sub>|</div><p class="simple-math-key"><strong>C(m)</strong> is the fraction of all absolute differences accounted for by the first m ranked proteins; the plot shows 100 × C(m), a percentage. <strong>m</strong> is how many proteins are included; <strong>j</strong> is a rank; <strong>D<sub>(j),k</sub></strong> is the j-th protein’s local A−B difference after sorting by largest |D| first.<br><strong>D<sub>i,k</sub></strong> is protein i’s local score with A minus its score with B, at the selected Qwen size; <strong>k = ${state.k}</strong>. <strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. <strong>|D|</strong> is the magnitude, ignoring the sign. <strong>∑</strong> means add. <strong>N = ${N.toLocaleString()}</strong> is the full reference-panel protein count. Both directions count; C(m) is undefined if every difference is zero.</p></details></section>`;}
  function differentialMarkup(){return `<div id="guided-differential-current"><p class="simple-comparison-key"><span><strong>Qwen:</strong> ${qs[checkpointIndex()].size_b}B · k = ${state.k}</span><span><strong>${esc(model(state.a))}</strong></span><span>versus <strong>${esc(model(state.b))}</strong></span></p>${card('Proteins with the largest differences at this Qwen size','<div class="simple-heatmap-scroll"><div id="guided-gap-proteins" class="simple-plot"></div></div><p class="simple-math-key">Up to 10 strongest differences in each direction; ties are ordered by protein ID. Click a bar to inspect a protein.</p>')}<p id="guided-pair-gap" class="simple-explainer"></p><div class="simple-equation" role="math" aria-label="The global alignment difference is the mean paired local alignment difference">α<sub>k</sub><sup>A</sup>(s) − α<sub>k</sub><sup>B</sup>(s) = (1/N) ∑<sub>i=1</sub><sup>N</sup> D<sub>i,k</sub>(s)</div><p class="simple-math-key"><strong>α</strong> is the global kNN alignment score, the mean of all local scores. <strong>D<sub>i,k</sub>(s)</strong> is protein i’s local score with A minus its local score with B. <strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>.<br><strong>N = ${N.toLocaleString()}</strong> is the total protein count in the fixed reference panel; <strong>∑</strong> adds all N differences before dividing by N. <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters.</p></div>`;}
  async function drawGlobal(){
    const values=qs.map((q,j)=>plms.map(p=>sum(series(p.id).rows.map(r=>r.counts[j]))/(N*state.k)));
    await plot('guided-global-heatmap',[{type:'heatmap',x:plms.map(p=>shortModel(p.id)),y:qs.map(q=>q.size_b+'B'),z:values,text:values.map(row=>row.map(v=>number(v,3))),texttemplate:'%{text}',customdata:qs.map(q=>plms.map(p=>`${q.id}|${p.id}`)),hovertemplate:'Qwen %{y}<br>%{x}<br>Global kNN alignment: %{z:.3g}<extra></extra>',zmin:0,colorscale:'Viridis',textfont:{color:'auto'},colorbar:{title:{text:'Global<br>alignment α'},thickness:12}}],{height:410,showlegend:false,xaxis:{title:{text:'Protein model'},tickangle:-35},yaxis:{title:{text:'Qwen model size'},autorange:'reversed'},margin:{l:95,r:90,t:16,b:125}});
    $('guided-global-heatmap').on?.('plotly_click',e=>{const point=e.points?.[0];if(!point?.customdata)return;const [q,p]=point.customdata.split('|');state.qwen=q;state.a=p;state.proteinPage=0;window.DESCRIPTIVE.ready=false;setTimeout(render,0);});
    await plot('guided-global-lines',plms.map((p,j)=>({type:'scatter',mode:'lines+markers',x:qs.map(q=>q.size_b+'B'),y:values.map(row=>row[j]),name:shortModel(p.id),line:{color:modelColors[j],dash:modelDash(p.id)},marker:{symbol:modelSymbol(p.id),size:6},hovertemplate:`${esc(model(p.id))}<br>Qwen %{x}<br>Global alignment: %{y:.3g}<extra></extra>`})),{height:450,xaxis:{title:{text:'Qwen model size'}},yaxis:{title:{text:`Global kNN alignment<br>(mean fraction, k=${state.k})`},rangemode:'tozero'},legend:{orientation:'h',y:1.25},margin:{l:70,r:20,t:100,b:60}});
  }
  async function drawDistribution(){
    await plot('guided-ecdf',plms.map((p,j)=>{const vals=localAt(p.id).sort((a,b)=>a-b),x=[],y=[];for(let i=0;i<N;i++)if(i===N-1||vals[i]!==vals[i+1]){x.push(vals[i]);y.push(100*(i+1)/N);}return {type:'scatter',mode:'lines+markers',x,y,name:shortModel(p.id),line:{color:modelColors[j],width:p.id===state.a||p.id===state.b?2.8:1.3,shape:'hv',dash:modelDash(p.id)},marker:{symbol:modelSymbol(p.id),size:4},hovertemplate:`${esc(model(p.id))}<br>Local alignment: %{x:.3g}<br>Proteins at or below: %{y:.3g}%<extra></extra>`};}),{height:460,xaxis:{title:{text:'Local alignment (shared-neighbor fraction)'},range:[0,1.02]},yaxis:{title:{text:'Proteins at or below score (%)'},range:[0,102]},legend:{orientation:'h',y:1.3},margin:{l:70,r:20,t:110,b:60}});
    const delta=pairCounts(),bins=Array.from({length:2*state.k+1},(_,j)=>j-state.k),frequencies=bins.map(b=>delta.filter(v=>v===b).length),ordered=delta.map((v,i)=>({v:Math.abs(v),i})).sort((a,b)=>b.v-a.v||ids[a.i].localeCompare(ids[b.i])),total=sum(ordered.map(r=>r.v));
    await plot('guided-gap-distribution',[{type:'bar',name:'Protein frequencies',showlegend:false,x:bins.map(c=>c/state.k),y:frequencies.map(n=>100*n/N),customdata:frequencies,marker:{color:bins.map(v=>v>0?A:v<0?B:'#9ba8b0')},hovertemplate:`${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: %{x:.3g}<br>%{customdata} proteins (%{y:.3g}%)<extra></extra>`},...pairLegend()],{height:430,showlegend:true,legend:{orientation:'h',y:1.25,font:{size:10}},xaxis:{title:{text:`Local score difference<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}`},range:[-1.05,1.05]},yaxis:{title:{text:'Proteins (%)'}},margin:{l:60,r:20,t:82,b:64}});
    if(total){let cumulative=0;const ys=[0,...ordered.map(r=>100*(cumulative+=r.v)/total)],xs=Array.from({length:N+1},(_,j)=>100*j/N);await plot('guided-gap-concentration',[{type:'scatter',mode:'lines',x:xs,y:ys,line:{color:A,width:2.5},hovertemplate:'Proteins included: %{x:.3g}%<br>Absolute difference accounted for: %{y:.3g}%<extra></extra>'},{type:'scatter',mode:'lines',x:[0,100],y:[0,100],line:{color:'#8998a2',dash:'dot'},hoverinfo:'skip'}],{height:380,showlegend:false,xaxis:{title:{text:'Proteins included (%)'},range:[0,100]},yaxis:{title:{text:'Total absolute difference (%)'},range:[0,102]},margin:{l:68,r:20,t:16,b:64}});}else $('guided-gap-concentration').innerHTML='<p class="empty">All paired local scores are equal; total absolute difference is zero.</p>';
  }
  async function drawDifferential(){const delta=pairCounts(),positive=ix.filter(i=>delta[i]>0).sort((i,j)=>delta[j]-delta[i]||ids[i].localeCompare(ids[j])).slice(0,10),negative=ix.filter(i=>delta[i]<0).sort((i,j)=>delta[i]-delta[j]||ids[i].localeCompare(ids[j])).slice(0,10),chosen=[...positive,...negative].sort((i,j)=>delta[j]-delta[i]||ids[i].localeCompare(ids[j]));
    if(chosen.length){await plot('guided-gap-proteins',[{type:'bar',name:'Protein differences',showlegend:false,orientation:'h',x:chosen.map(i=>delta[i]/state.k),y:chosen.map(i=>ids[i]),customdata:chosen,text:chosen.map(i=>label(i)),marker:{color:chosen.map(i=>delta[i]>0?A:B)},hovertemplate:`%{text}<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: %{x:.3g}<extra></extra>`},...pairLegend()],{height:Math.max(330,chosen.length*25+160),showlegend:true,legend:{orientation:'h',y:1.12,font:{size:10}},xaxis:{title:{text:`Local score difference<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}`},range:[-1.03,1.03],zeroline:true},yaxis:{title:{text:'Gene / protein name'},tickvals:chosen.map(i=>ids[i]),ticktext:chosen.map(i=>esc(label(i))),autorange:'reversed'},margin:{l:190,r:30,t:85,b:60}});$('guided-gap-proteins').on?.('plotly_click',e=>{const i=e.points?.[0]?.customdata;if(Number.isInteger(i)){state.focus=i;window.DESCRIPTIVE.ready=false;setTimeout(render,0);}});}else $('guided-gap-proteins').innerHTML='<p class="empty">All paired local scores are equal at this checkpoint.</p>';
    const positiveGap=sum(delta.filter(v=>v>0))/(N*state.k),negativeGap=sum(delta.filter(v=>v<0))/(N*state.k);$('guided-pair-gap').innerHTML=`Global ${esc(shortModel(state.a))} − ${esc(shortModel(state.b))} gap: <strong>${signed(100*(positiveGap+negativeGap),3)} pp</strong> = ${signed(100*positiveGap,3)} pp from proteins favoring ${esc(shortModel(state.a))} + ${signed(100*negativeGap,3)} pp from proteins favoring ${esc(shortModel(state.b))}. All ${N.toLocaleString()} proteins enter this sum.`;
  }

/* Insert inside dashboard_descriptive.js's closure. Original counts; no metric refits. */
function topProteinData() {
  const qIndex = qs.findIndex(q => q.id === state.qwen);
  if (qIndex < 0) throw Error('Choose a recorded Qwen checkpoint for the top-protein view.');
  const requested = Number(state.topN ?? 20);
  if (![10, 20, 50, 100].includes(requested)) throw Error('Choose 10, 20, 50, or 100 top proteins.');
  const annotations = D.panel.annotations;
  const source = series(state.a).rows;
  const ordered = ix.map(i => {
    const count = source[i].counts[qIndex];
    const raw = String(annotations.group?.[i] ?? '').trim();
    const familyId = !raw || ['unknown', 'none', 'null', 'nan'].includes(raw.toLowerCase()) ? 'unknown' : raw;
    const catalog = D.panel.group_names?.[familyId];
    const familyShortName = familyId === 'unknown' ? '' : String(catalog?.name || annotations.group_name?.[i] || '');
    const familyName = familyId === 'unknown' ? 'No recorded Pfam family'
      : String(catalog?.description || catalog?.name || annotations.group_description?.[i] || annotations.group_name?.[i] || 'Unlabeled Pfam family');
    return {i, count, local: count / state.k, familyId, familyName, familyShortName, familyLabel: familyName};
  }).sort((a, b) => b.count - a.count || ids[a.i].localeCompare(ids[b.i]));
  const cutoffCount = ordered[Math.min(requested, N) - 1].count;
  const rows = ordered.filter(row => row.count >= cutoffCount);
  let previousCount = null, rank = 0;
  rows.forEach((row, index) => {if (row.count !== previousCount) rank = index + 1; row.rank = rank; previousCount = row.count;});
  const groups = new Map();
  for (const row of rows) {
    if (!groups.has(row.familyId)) groups.set(row.familyId, {id: row.familyId, name: row.familyName, shortName: row.familyShortName, label: row.familyLabel, n: 0});
    groups.get(row.familyId).n++;
  }
  const labelCounts = new Map();
  for (const group of groups.values()) labelCounts.set(group.label, (labelCounts.get(group.label) || 0) + 1);
  for (const group of groups.values()) {
    if (labelCounts.get(group.label) > 1 && group.shortName) group.label += ` — ${group.shortName}`;
  }
  for (const row of rows) row.familyLabel = groups.get(row.familyId).label;
  const families = [...groups.values()].sort((a, b) => b.n - a.n || a.id.localeCompare(b.id));
  let displayFamilies = families;
  let otherGroupCount = 0;
  let familyDisplayNote = 'Every recorded group in the top-protein set is shown.';
  if (families.length > 20) {
    const unknown = families.find(group => group.id === 'unknown');
    const named = families.filter(group => group.id !== 'unknown');
    const namedLimit = unknown ? 19 : 20;
    displayFamilies = named.slice(0, namedLimit);
    if (unknown) displayFamilies.push(unknown);
    const remainder = named.slice(namedLimit);
    otherGroupCount = remainder.length;
    displayFamilies.push({id: '__other_recorded_groups__', name: 'Other recorded groups',
      label: `Other recorded families (${otherGroupCount} groups)`, n: remainder.reduce((total, group) => total + group.n, 0)});
    familyDisplayNote = `The ${namedLimit} largest annotated groups${unknown ? ' and unknown annotations' : ''} are shown separately; “Other” combines the remaining ${otherGroupCount} groups.`;
  }
  if (families.reduce((total, group) => total + group.n, 0) !== rows.length ||
      displayFamilies.reduce((total, group) => total + group.n, 0) !== rows.length) {
    throw Error('Family counts must cover the exact top-protein set.');
  }
  return {q: qs[qIndex], requested, cutoffCount, rows, families,
    displayFamilies, otherGroupCount, familyDisplayNote, n: rows.length, tieExtra: rows.length - Math.min(requested, N)};
}

function topMarkup() {
  const data = topProteinData();
  const comparison = `Qwen ${data.q.size_b}B × ${shortModel(state.a)} · k = ${state.k}`;
  const cutoff = `${data.cutoffCount}/${state.k} shared neighbors`;
  const rows = data.rows.map(row => [row.rank, `<button class="protein-link" data-top-focus="${row.i}">${esc(label(row.i))}</button>`,
    esc(row.familyLabel), `${row.count}/${state.k}`, row.local.toFixed(2), D.panel.annotations.length?.[row.i] ?? 'Not recorded', functionCell(row.i)]);
  const scoreGroups = new Map();
  for (const row of data.rows) {
    if (!scoreGroups.has(row.count)) scoreGroups.set(row.count, []);
    scoreGroups.get(row.count).push(row);
  }
  const geneList = `<div id="guided-top-protein-list" class="top-score-list">${[...scoreGroups].map(([count, members]) =>
    `<section class="top-score-group" data-shared-count="${count}"><h4 class="top-score-header"><span class="top-score-value">${number(count/state.k,2)}</span><span>local alignment <span class="top-score-detail">${count}/${state.k} shared neighbors · ${members.length} protein${members.length===1?'':'s'}</span></span></h4><ul class="top-gene-list">${members.map(row=>
      `<li><button class="top-gene-name" data-top-focus="${row.i}" data-protein-id="${esc(ids[row.i])}" aria-pressed="${row.i===state.focus}" title="${esc(ids[row.i]+' · '+row.familyLabel)}">${esc(label(row.i))}</button></li>`).join('')}</ul></section>`).join('')}</div>`;
  return `<section class="simple-section" id="guided-top-proteins">
    <h2>Which proteins have the highest local kNN alignment, and which families contain them?</h2>
    <p class="simple-explainer"><strong>${esc(comparison)}</strong></p>
    <div class="simple-actions"><label>Top proteins to display <select id="guided-top-n" aria-label="Requested number of top proteins">${[10,20,50,100].map(n => `<option value="${n}" ${n === data.requested ? 'selected' : ''}>${n}</option>`).join('')}</select></label><button id="guided-top-track">Track these ${data.n} proteins</button><button id="guided-top-export">All ${data.n} top proteins · CSV</button></div>
    <p class="simple-selection-note" id="guided-top-cohort"><strong>Top ${data.requested} + ties = ${data.n} proteins.</strong> Include everyone at the cutoff: ${cutoff}. Names within each score group are ordered by protein ID.</p><p class="simple-selection-note">Labels show gene names. Parenthesized IDs distinguish repeated gene labels; they do not imply identical sequences. A bare ID means no gene name was recorded.</p>
    <div class="simple-equation" role="math" aria-label="Local alignment equals shared Qwen and PLM neighbors divided by k">ℓ<sub>i,k</sub><sup>P</sup>(s) = |𝒩<sub>k</sub><sup>Q(s)</sup>(i) ∩ 𝒩<sub>k</sub><sup>P</sup>(i)| / k</div><p class="simple-math-key"><strong>ℓ</strong> is this protein’s local alignment score: the fraction of its neighbors shared by Qwen and the PLM (0 to 1). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the number of neighbors; <strong>s</strong> is Qwen size in billions of parameters; <strong>P</strong> is the PLM checkpoint.<br><strong>Q(s)</strong> is Qwen at size s. <strong>𝒩</strong> is a set of nearest neighbors of i in the indicated model. <strong>∩</strong> means the shared members; <strong>|·|</strong> counts them. The protein itself is excluded. Superscripts identify models, not powers.</p>
    ${card(`Which proteins have the highest local kNN alignment with ${shortModel(state.a)}?`, geneList)}
    ${card('2. How many of these proteins share a recorded family?', `<p class="simple-explainer">The same <strong>${data.n} proteins</strong>, grouped by their original single recorded Pfam assignment. Counts add to ${data.n}. ${esc(data.familyDisplayNote)}</p><div class="simple-heatmap-scroll" style="max-height:650px;overflow:auto"><div id="guided-top-family-plot" class="simple-plot" aria-label="Recorded Pfam group counts for the exact same top-protein set"></div></div><p class="simple-explainer">Counts use the first Pfam accession in sorted order for each protein; additional Pfam annotations are not represented here. Missing annotations are grouped separately.</p>`)}
    <details id="guided-top-table-details"><summary>All ${data.n} proteins: family, length and function</summary><div style="max-height:460px;overflow:auto">${table(['Rank', 'Gene / protein name', 'Pfam family', 'Shared neighbors', 'Local alignment', 'Length (amino acids)', 'Function'], rows, 'guided-top-table')}</div></details>
    <p class="simple-explainer">“Track these proteins” carries this exact set into the trajectory and comparison views. Changing the displayed checkpoint recalculates this ranking, while the tracked IDs stay fixed until you choose a new set.</p>
  </section>`;
}

async function drawTop() {
  const data = topProteinData();
  const wrapFamily = value => {
    const lines = [''];
    for (const word of value.split(/\s+/)) {
      const last = lines.length - 1;
      if (lines[last] && lines[last].length + word.length + 1 > 38) lines.push(word);
      else lines[last] += (lines[last] ? ' ' : '') + word;
    }
    return lines.map(esc).join('<br>');
  };
  const familyTicks = data.displayFamilies.map(group => wrapFamily(group.label));
  const familyRowHeight = Math.max(29, ...familyTicks.map(value => (value.split('<br>').length * 13) + 6));
  await plot('guided-top-family-plot', [{type: 'bar', orientation: 'h',
    x: data.displayFamilies.map(group => group.n), y: data.displayFamilies.map(group => group.id),
    text: data.displayFamilies.map(group => `${esc(group.label)}<br>${group.n} of ${data.n} top proteins${group.id.startsWith('__') ? '' : `<br>Pfam accession: ${esc(group.id)}`}`),
    hovertemplate: '%{text}<extra></extra>', textposition: 'outside', cliponaxis: false,
    texttemplate: '%{x}', marker: {color: data.displayFamilies.map(group => group.id === '__other_recorded_groups__' ? '#a2aeb7' : group.id === 'unknown' ? '#667985' : colorOf(state.a))}}],
    {height: Math.max(280, data.displayFamilies.length * familyRowHeight + 90), showlegend: false,
      xaxis: {title: {text: 'Number of top proteins'}, side: 'top', rangemode: 'tozero', tickformat: 'd', dtick: Math.max(1, Math.ceil(Math.max(...data.displayFamilies.map(group => group.n)) / 5))},
      yaxis: {title: {text: 'Pfam family'}, type: 'category', categoryorder: 'array', categoryarray: data.displayFamilies.map(group => group.id), autorange: 'reversed',
        tickvals: data.displayFamilies.map(group => group.id), ticktext: familyTicks},
      margin: {l: 270, r: 42, t: 68, b: 20}});
}

function bindTop() {
  const data = topProteinData();
  const select = document.getElementById('guided-top-n');
  if (select) select.onchange = event => {state.topN = Number(event.target.value); render();};
  const track = document.getElementById('guided-top-track');
  if (track) track.onclick = () => {
    state.selected = new Set(data.rows.map(row => row.i));
    state.focus = data.rows[0].i;
    state.selectionOrigin = `Top ${data.requested} plus boundary ties: ${data.n} proteins for Qwen ${data.q.size_b}B and ${model(state.a)}, k=${state.k}; cutoff ${data.cutoffCount}/${state.k}. These exact IDs stay fixed across comparison changes.`;
    render();
  };
  const exportButton = document.getElementById('guided-top-export');
  if (exportButton) exportButton.onclick = () => download('top_proteins_and_families.csv', data.rows.map(row => ({
    panel_id: D.panel.id, panel_n: N, qwen_checkpoint: data.q.id, plm_checkpoint: state.a, k: state.k,
    requested_top_n: data.requested, top_cohort_n: data.n, cutoff_shared_neighbors: data.cutoffCount,
    rank: row.rank, protein_id: ids[row.i], gene_name: D.panel.annotations.gene_name?.[row.i] || '',
    protein_label: label(row.i), protein_name: D.panel.annotations.protein_name?.[row.i] || '',
    recorded_pfam_id: row.familyId, recorded_pfam_name: row.familyShortName,
    recorded_pfam_description: row.familyName, recorded_family_label: row.familyLabel,
    sequence_length: D.panel.annotations.length?.[row.i] ?? '', function: D.panel.annotations.function?.[row.i] ?? '',
    shared_neighbors: row.count, local_alignment: row.local
  })));
  document.querySelectorAll('[data-top-focus]').forEach(button => {
    button.onclick = () => {state.focus = Number(button.dataset.topFocus); render();};
  });

}

/* Insert inside dashboard_descriptive.js's closure. Uses the existing counts and stored controls only. */
function biologyData() {
  const qIndex = qs.findIndex(q => q.id === state.qwen);
  if (qIndex < 0) throw Error('Choose a recorded Qwen checkpoint for the biology view.');
  const a = series(state.a).rows, b = series(state.b).rows;
  const annotations = D.panel.annotations;
  const rows = ix.map(i => ({
    i, length: annotations.length?.[i],
    a: a[i].counts[qIndex] / state.k, b: b[i].counts[qIndex] / state.k,
    gap: (a[i].counts[qIndex] - b[i].counts[qIndex]) / state.k,
    taxon: String(annotations.taxon?.[i] || 'Unknown collection group'),
    family: String(annotations.group?.[i] || 'unknown')
  }));
  const group = field => {
    const groups = new Map();
    for (const row of rows) {
      const key = row[field];
      if (!groups.has(key)) {
        const missing = ['unknown', 'none', 'null', 'nan', ''].includes(key.toLowerCase());
        const catalog = D.panel.group_names?.[key];
        const shortName = field === 'family' ? String(catalog?.name || annotations.group_name?.[row.i] || '') : '';
        const name = field === 'family' ? (missing ? 'No recorded Pfam family' :
          String(catalog?.description || catalog?.name || annotations.group_description?.[row.i] || annotations.group_name?.[row.i] || 'Unlabeled Pfam family')) : key;
        groups.set(key, {id: key, label: name, shortName, n: 0, a: 0, b: 0, gap: 0});
      }
      const g = groups.get(key);
      g.n++; g.a += row.a; g.b += row.b; g.gap += row.gap;
    }
    const labelCounts = new Map();
    for (const group of groups.values()) labelCounts.set(group.label, (labelCounts.get(group.label) || 0) + 1);
    for (const group of groups.values()) {
      if (labelCounts.get(group.label) > 1 && group.shortName) group.label += ` — ${group.shortName}`;
    }
    return [...groups.values()].map(g => ({...g, a: g.a / g.n, b: g.b / g.n, gap: g.gap / g.n}))
      .sort((a, b) => a.label.localeCompare(b.label));
  };
  return {q: qs[qIndex], rows, families: group('family'), taxa: group('taxon')};
}

function biologyMarkup() {
  const data = biologyData();
  const fmt = value => Number.isFinite(value) ? number(value) : 'Unavailable';
  const groupRows = groups => groups.map(g => [esc(g.label), g.n, fmt(g.a), fmt(g.b), fmt(g.gap)]);
  return `<section class="simple-section" id="guided-question-4">
    <h2>How do mean local kNN differences vary by Pfam family and collection group?</h2>
    <p class="simple-explainer">The 12 largest recorded families, using the same ${N.toLocaleString()}-protein panel.</p>
    <div class="simple-equation" role="math" aria-label="The local alignment difference is the score with A minus the score with B for the same protein">D<sub>i,k</sub>(s) = ℓ<sub>i,k</sub><sup>A</sup>(s) − ℓ<sub>i,k</sub><sup>B</sup>(s)</div><p class="simple-math-key"><strong>D</strong> is one protein’s local alignment difference, A minus B. <strong>ℓ</strong> is its local alignment score (shared neighbors divided by k). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters.<br><strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. Positive D favors A; negative D favors B; zero means equal local scores.</p>
    ${card('The 12 largest recorded Pfam groups', '<div class="simple-heatmap-scroll" tabindex="0"><div id="guided-biology-family" class="simple-plot" aria-label="Family mean local alignment differences shown as signed points from zero"></div></div><p class="simple-explainer">Families are chosen by protein count, not by their score. Each protein uses its first Pfam accession in sorted order; additional annotations are not represented here.</p>')}
    <details id="guided-biology-taxon-details"><summary>Collection groups and how group means are calculated</summary>
      ${card('Differences by collection group', '<div id="guided-biology-taxon" class="simple-plot" aria-label="Mean local alignment difference by recorded collection group"></div><p class="simple-explainer">Colors identify collection groups, not the sign of the difference.</p>')}
      <div class="simple-equation" role="math" aria-label="A group mean difference is the sum of the protein differences in that group divided by the group's protein count">D̄<sub>G,k</sub>(s) = (1 / n<sub>G</sub>) ∑<sub>i∈G</sub> D<sub>i,k</sub>(s)</div>
      <p class="simple-math-key"><strong>D̄</strong> (D-bar) is the group’s mean local alignment difference. <strong>G</strong> is the annotation group; <strong>n<sub>G</sub></strong> is the number of proteins in that group. <strong>∑<sub>i∈G</sub></strong> means add over proteins i belonging to G.<br><strong>D<sub>i,k</sub>(s)</strong> is protein i’s local score with A minus its local score with B. <strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters. A group mean divides by n<sub>G</sub>, not the full-panel count N = ${N.toLocaleString()}.</p>
      ${table(['Collection group', 'Proteins', `Mean ${shortModel(state.a)} score`, `Mean ${shortModel(state.b)} score`, `Mean ${shortModel(state.a)} − ${shortModel(state.b)}`], groupRows(data.taxa), 'guided-biology-taxa-table')}
    </details>
    <details id="guided-biology-family-details"><summary>Browse all ${data.families.length.toLocaleString()} named Pfam groups</summary><div style="max-height:360px;overflow:auto">${table(['Recorded Pfam group', 'Proteins', `Mean ${shortModel(state.a)} score`, `Mean ${shortModel(state.b)} score`, `Mean ${shortModel(state.a)} − ${shortModel(state.b)}`], groupRows(data.families), 'guided-biology-families-table')}</div></details>
    <p class="simple-explainer"><strong>Available next:</strong> recorded GO annotations and PDB cross-references remain in the source manifest. Structural measurements and a fitted functional-property predictor are not part of this report.</p>
  </section>`;
}

async function drawBiology() {
  const data = biologyData();
  const palette = ['#0072b2', '#d55e00', '#009e73', '#8d3d91', '#778896'];
  const fmt = value => Number.isFinite(value) ? number(value) : 'Unavailable';
  await plot('guided-biology-taxon', [{type: 'bar', x: data.taxa.map(g => g.label), y: data.taxa.map(g => g.gap),
    text: data.taxa.map(g => `${esc(g.label)}<br>Proteins: ${g.n}<br>${esc(shortModel(state.a))} mean: ${fmt(g.a)} · ${esc(shortModel(state.b))} mean: ${fmt(g.b)}<br>Mean ${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: ${fmt(g.gap)}`),
    hovertemplate: '%{text}<extra></extra>', textposition: 'none', marker: {color: data.taxa.map((_, j) => palette[j % palette.length])}}],
    {height: 330, showlegend: false, xaxis: {title: {text: 'Recorded collection group'}},
      yaxis: {title: {text: `Mean local score difference<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}`}, nticks: 5, tickformat: '.2f', zeroline: true, zerolinecolor: '#193243'},
      margin: {l: 78, r: 20, t: 12, b: 70}});
  const largest = [...data.families].sort((a, b) => b.n - a.n || a.id.localeCompare(b.id)).slice(0, 12);
  const familyTicks = largest.map(group => {
    const lines = [''];
    for (const word of group.label.split(/\s+/)) {
      const last = lines.length - 1;
      if (lines[last] && lines[last].length + word.length + 1 > 38) lines.push(word);
      else lines[last] += (lines[last] ? ' ' : '') + word;
    }
    return lines.map(esc).join('<br>') + ` · n=${group.n}`;
  });
  const extent = Math.max(.05, ...largest.map(g => Math.abs(g.gap))) * 1.35;
  const tickBase = 10 ** Math.floor(Math.log10(extent / 3));
  const tickStep = [1,2,5,10].map(v=>v*tickBase).find(v=>v>=extent/3);
  const tickCount = Math.floor(extent/tickStep);
  const familyTickValues = Array.from({length:2*tickCount+1},(_,j)=>(j-tickCount)*tickStep);
  const familyColors = largest.map(g => g.gap > 0 ? colorOf(state.a) : g.gap < 0 ? colorOf(state.b) : '#667985');
  const familyShapes = largest.map(g => g.gap > 0 ? 'triangle-right' : g.gap < 0 ? 'triangle-left' : 'diamond');
  const rowShapes = largest.flatMap((g, j) => [
    {type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:j-.5,y1:j+.5,fillcolor:j%2?'#ffffff':'#f4f7fa',line:{width:0},layer:'below'},
    {type:'line',xref:'x',yref:'y',x0:0,x1:g.gap,y0:j,y1:j,line:{color:familyColors[j],width:2},opacity:.5,layer:'below'}
  ]);
  rowShapes.push({type:'line',xref:'x',yref:'paper',x0:0,x1:0,y0:0,y1:1,line:{color:'#566e7b',width:1.5},layer:'below'});
  await plot('guided-biology-family', [{type:'scatter',mode:'markers+text',showlegend:false,
    x:largest.map(g=>g.gap),y:largest.map(g=>g.id),
    text:largest.map(g=>signed(g.gap,3)),textposition:largest.map(g=>g.gap<0?'middle left':'middle right'),
    textfont:{size:12,color:'#193243'},cliponaxis:false,
    hovertext:largest.map(g=>`${esc(g.label)}<br>Proteins: ${g.n}<br>${esc(shortModel(state.a))}: ${fmt(g.a)}<br>${esc(shortModel(state.b))}: ${fmt(g.b)}<br>Mean ${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: ${signed(g.gap,3)}<br>Pfam accession: ${esc(g.id)}`),
    hovertemplate:'%{hovertext}<extra></extra>',
    marker:{color:familyColors,symbol:familyShapes,size:12,line:{color:'#ffffff',width:1.4}}},
    {type:'scatter',mode:'markers',x:[null],y:[null],name:`Negative: ${shortModel(state.b)} higher`,marker:{color:colorOf(state.b),symbol:'triangle-left',size:10},hoverinfo:'skip'},
    {type:'scatter',mode:'markers',x:[null],y:[null],name:`Positive: ${shortModel(state.a)} higher`,marker:{color:colorOf(state.a),symbol:'triangle-right',size:10},hoverinfo:'skip'}],
    {height:Math.max(570,largest.length*Math.max(37,...familyTicks.map(value=>value.split('<br>').length*14+9))+145),
      showlegend:true,legend:{orientation:'h',x:0,y:1.15,font:{size:12}},
      xaxis:{title:{text:`Mean local score difference<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}`},range:[-extent,extent],tickmode:'array',tickvals:familyTickValues,ticktext:familyTickValues.map(v=>number(v,3)),zeroline:false,showgrid:true,gridcolor:'#e2eaf0',ticks:'outside',tickfont:{size:12}},
      yaxis:{title:{text:'Pfam family'},type:'category',categoryorder:'array',categoryarray:largest.map(g=>g.id),autorange:'reversed',
        tickvals:largest.map(g=>g.id),ticktext:familyTicks,tickfont:{size:12},showgrid:false,ticks:''},
      shapes:rowShapes,annotations:[{text:'Equal alignment',x:0,y:1.025,xref:'x',yref:'paper',showarrow:false,font:{size:11,color:'#536976'},xanchor:'center'}],
      margin:{l:290,r:60,t:92,b:75}});
}

function redundancyData() {
  const metric = `knn_k${state.k}`;
  const stored = D.review?.name_controls;
  if (!Array.isArray(stored)) return null;
  const rows = [state.a, state.b].map(p => ({p,
    values: qs.map(q => stored.find(row => row.plm_model === p && row.qwen_model === q.id && row.metric === metric))}));
  return rows.some(row => row.values.some(Boolean)) ? rows : null;
}

function redundancyMarkup() {
  const data = redundancyData();
  const fmt = value => Number.isFinite(value) ? number(value) : 'Unavailable';
  const present = data ? data.flatMap(group => group.values.filter(Boolean)) : [];
  const counts = [...new Set(present.map(row => row.retained_n))].join(', ');
  const draws = [...new Set(present.map(row => row.draws))].join(', ');
  const result = !data ? card('Name controls are unavailable in this build', '<p>No stored results match the selected checkpoints and neighbor count.</p>') :
    `<p class="simple-explainer">Compare the full panel with one-per-protein-name panels and organism-matched random deletion. Control lines average ${esc(draws)} recorded draws.</p>
    ${card(`Does the ${shortModel(state.a)} − ${shortModel(state.b)} difference persist in the name controls?`, '<div id="guided-redundancy-gap" class="simple-plot" aria-label="Global alignment difference between the selected PLMs under the stored name controls"></div><p class="simple-explainer">Line colors identify controls, not PLMs. Compare one-per-name with matched deletion to distinguish name selection from simply having fewer proteins.</p>')}
    <div class="simple-equation" role="math" aria-label="The model gap for a control is the average across the matched draws of the score with A minus the score with B">Gap<sub>control,k</sub>(s) = mean<sub>draws</sub>[α<sub>control,k</sub><sup>A</sup>(s) − α<sub>control,k</sub><sup>B</sup>(s)]</div><p class="simple-math-key"><strong>Gap</strong> is the mean A−B difference across recorded random draws. <strong>α</strong> is the global kNN alignment score within a control’s retained protein panel. <strong>control</strong> identifies the one-per-name or matched-deletion procedure; a <strong>draw</strong> is one random selection of proteins; <strong>mean</strong> averages the draw values.<br><strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters. Each control score averages its retained proteins; the original full panel contains N = ${N.toLocaleString()}.</p>
    <details id="guided-redundancy-paired-details"><summary>One-per-name versus matched deletion for each PLM</summary>
      ${card('One-per-name selection compared with matched deletion', '<div id="guided-redundancy-paired" class="simple-plot" aria-label="Stored paired differences between one-per-name and matched deletion for each selected PLM"></div>')}
      <div class="simple-equation" role="math" aria-label="Within each draw, the name-control difference is the one-per-name score minus the matched-deletion score">δ<sub>draw,k</sub><sup>P</sup>(s) = α<sub>one-per-name,k</sub><sup>P</sup>(s) − α<sub>matched,k</sub><sup>P</sup>(s)</div><p class="simple-math-key"><strong>δ (delta)</strong> is the control-score difference within one paired random draw. <strong>α</strong> is global kNN alignment in the indicated retained panel. <strong>one-per-name</strong> retains one random protein per full protein name; <strong>matched</strong> retains the same number per organism taxon by random deletion.<br><strong>P</strong> is the PLM checkpoint; <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters. Positive δ means higher alignment in the one-per-name panel.</p>
    </details>
    <details id="guided-redundancy-method-details"><summary>How the name controls were run</summary><p class="simple-explainer">Each of ${esc(draws)} draws retains one random protein per exact recorded full protein name (${esc(counts)} proteins per panel). Its random-deletion comparator retains exactly the same count in every organism taxon ID. These groups use full protein names, not the gene labels shown in the plots. The same selections are applied to all checkpoints. Each smaller candidate pool is refitted using submatrices of the original fixed similarities. Full-panel scores still use all ${N.toLocaleString()} original proteins.</p></details>
    <details id="guided-redundancy-table-details"><summary>Recorded scores for both selected PLMs</summary>${table(['PLM', 'Qwen size', 'Full panel', 'One/name mean', 'Matched mean', 'Paired difference mean', 'Draw minimum', 'Draw maximum'], data.flatMap((group, j) => group.values.map((row, index) => [esc(model(group.p)), esc(qs[index].size_b + 'B'), fmt(row?.full_score), fmt(row?.one_per_name_mean), fmt(row?.taxon_matched_mean), fmt(row?.paired_difference_mean), fmt(row?.paired_difference_min), fmt(row?.paired_difference_max)])), 'guided-redundancy-table')}</details>`;
  return `<section class="simple-section" id="guided-question-5"><h2>How does global kNN alignment compare across protein-name controls?</h2>
    <p class="simple-explainer"><strong>Completed:</strong> exact-name controls shown below. <strong>Not run:</strong> family-balanced scores and sequence-identity-filtered comparisons.</p>
    <p class="simple-comparison-key"><span><strong>${esc(model(state.a))}</strong></span><span>versus <strong>${esc(model(state.b))}</strong></span></p>
    ${result}
    <p class="simple-explainer">The existing component bootstrap resamples related proteins together; it does not equalize family weights or remove related neighbors. The homology summaries describe existing pair contributions. Testing family balancing or identity filtering still requires a family/identity rule and a decision about whether to change the candidate panel.</p>
  </section>`;
}

async function drawRedundancy() {
  const data = redundancyData();
  if (!data) return;
  const x = qs.map(q => q.size_b + 'B');
  const difference = field => qs.map((_, j) => {
    const a = data[0].values[j]?.[field], b = data[1].values[j]?.[field];
    return Number.isFinite(a) && Number.isFinite(b) ? a - b : null;
  });
  const settings = [
    ['Full panel', 'full_score', '#193243', 'solid'],
    ['One per exact protein name', 'one_per_name_mean', '#0072b2', 'dash'],
    ['Organism-matched deletion', 'taxon_matched_mean', '#d55e00', 'dot']
  ];
  await plot('guided-redundancy-gap', settings.map(([name, field, color, dash]) => ({type: 'scatter', mode: 'lines+markers',
    name, x, y: difference(field), connectgaps: false, line: {color, dash, width: 2}, marker: {size: 7},
    hovertemplate: `${esc(name)}<br>Qwen %{x}<br>Global ${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: %{y:.3g}<extra></extra>`})),
    {height: 340, xaxis: {title: {text: 'Qwen model size'}, type: 'category', categoryorder: 'array', categoryarray: x},
      yaxis: {title: {text: `Global score difference<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}`}, nticks: 5, tickformat: '.3~g', zeroline: true, zerolinecolor: '#193243'},
      margin: {l: 80, r: 20, t: 48, b: 68}});
  await plot('guided-redundancy-paired', data.map((group, j) => ({type: 'scatter', mode: 'lines+markers',
    name: shortModel(group.p), x,
    y: group.values.map(row => row?.paired_difference_mean ?? null), connectgaps: false,
    error_y: {type: 'data', symmetric: false,
      array: group.values.map(row => row ? Math.max(0, row.paired_difference_max - row.paired_difference_mean) : null),
      arrayminus: group.values.map(row => row ? Math.max(0, row.paired_difference_mean - row.paired_difference_min) : null),
      visible: true, thickness: 1.3, width: 4},
    text: group.values.map(row => row ? `${esc(model(group.p))}<br>${row.draws} recorded paired draws<br>Mean: ${number(row.paired_difference_mean)}<br>Minimum: ${number(row.paired_difference_min)}<br>Maximum: ${number(row.paired_difference_max)}` : 'Unavailable'),
    hovertemplate: '%{text}<extra></extra>', line: {color: colorOf(group.p), width: 2, dash:modelDash(group.p)}, marker: {size: 7,symbol:modelSymbol(group.p)}})),
    {height: 340, xaxis: {title: {text: 'Qwen model size'}, type: 'category', categoryorder: 'array', categoryarray: x},
      yaxis: {title: {text: 'Score difference<br>(one/name − matched deletion)'}, nticks: 5, tickformat: '.3~g', zeroline: true, zerolinecolor: '#193243'},
      margin: {l: 88, r: 20, t: 48, b: 68}});
}

  async function renderNow() {
    const openDetails=[...document.querySelectorAll('#simple-view details[open][id]')].map(el=>el.id);
    for(const id of plots){if($(id))Plotly.purge($(id));}plots.clear();
    for(const [id,key] of [['simple-a','a'],['simple-b','b'],['simple-k','k'],['simple-search','search'],['simple-pattern','pattern'],['simple-sort','sort']])$(id).value=String(state[key]);
    const a=series(state.a),b=series(state.b),selected=[...state.selected],selectedRows=selected.slice(0,30),differenceRows=selected.slice(0,40),focus=state.focus;
    const trajectoryRows=[...differenceRows.filter(i=>i!==focus),focus];
    $('simple-summary').innerHTML=`<p class="simple-explainer"><strong>The same ${N.toLocaleString()} proteins · k=${state.k} · Qwen ${qs.map(q=>q.size_b+'B').join(' → ')}</strong><br>Local alignment is the fraction of nearest neighbors shared by Qwen and the selected PLM.</p><div class="simple-equation" role="math" aria-label="Local alignment equals the number of neighbors shared by Qwen and PLM P, divided by k">ℓ<sub>i,k</sub><sup>P</sup>(s) = |𝒩<sub>k</sub><sup>Q(s)</sup>(i) ∩ 𝒩<sub>k</sub><sup>P</sup>(i)| / k</div><p class="simple-math-key"><strong>ℓ</strong> is this protein’s local alignment score: the fraction of its neighbors shared by Qwen and the PLM (0 to 1). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the number of neighbors; <strong>s</strong> is Qwen size in billions of parameters; <strong>P</strong> is the PLM checkpoint.<br><strong>Q(s)</strong> is Qwen at size s. <strong>𝒩</strong> is a set of nearest neighbors of i in the indicated model. <strong>∩</strong> means the shared members; <strong>|·|</strong> counts them. The protein itself is excluded. Superscripts identify models, not powers.</p><p class="simple-explainer">All comparisons retain the original candidates, layers, pooling, cosine similarity, and self-exclusion: i is excluded from its own neighbor sets. The global kNN alignment score (k = ${state.k}) is the mean over all ${N.toLocaleString()} proteins:</p><div class="simple-equation" role="math" aria-label="The global kNN alignment score alpha is the sum of all protein local alignment scores divided by N">α<sub>k</sub><sup>P</sup>(s) = (1/N) ∑<sub>i=1</sub><sup>N</sup> ℓ<sub>i,k</sub><sup>P</sup>(s)</div><p class="simple-math-key"><strong>α</strong> is the global kNN alignment score. <strong>ℓ</strong> is one protein’s local alignment score (shared neighbors divided by k). <strong>N = ${N.toLocaleString()}</strong> is the total number of proteins in the fixed reference panel, including proteins hidden by display filters or outside the selection.<br><strong>i</strong> identifies a protein; <strong>∑</strong> means add the scores for all N proteins, then divide by N. <strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint; <strong>s</strong> is Qwen size in billions of parameters.</p><div class="simple-summary-grid">${[state.a,state.b].map(p=>{const s=series(p),current=sum(s.rows.map(r=>r.counts[checkpointIndex()]))/(N*state.k);return `<div class="summary-box"><span class="summary-label">${esc(model(p))} · global kNN alignment (k = ${state.k})</span><strong class="summary-value">${number(current,4)}</strong><span class="summary-detail">Qwen ${qs[checkpointIndex()].size_b}B · ${N.toLocaleString()} proteins</span></div>`;}).join('')}</div>`;
    const summaries=[state.a,state.b].map(p=>{const r=series(p).rows;return [esc(model(p)),r.filter(x=>x.pattern==='gain').length,r.filter(x=>x.pattern==='retained').length,r.filter(x=>x.pattern==='loss').length];});
    const focusRows=[state.a,state.b].map(p=>{const r=series(p).rows[focus];return [esc(model(p)),...r.values.map(v=>number(v,2)),signed(r.delta,2),esc(patternName(r.pattern))];});
    const accounts=[state.a,state.b].map(p=>({p,...accounting(p)}));
    const accountRows=accounts.map(s=>[esc(model(s.p)),signed(s.positive*100,3),signed(s.negative*100,3),signed(s.net*100,3),signed(s.selectedNet*100,3),signed(s.remaining*100,3),s.selectedShare===null?'Undefined':signed(s.selectedShare*100,1)+'%']);
    $('simple-view').innerHTML=`<section class="simple-section" id="simple-question-1"><h2>Which proteins gain or lose local kNN alignment across Qwen sizes?</h2><p class="simple-explainer">Start with each protein’s five measured scores. Gains and losses describe the endpoints; the full row shows any intermediate dips or recoveries.</p><div class="simple-equation" role="math" aria-label="Local change equals local alignment at the largest Qwen minus local alignment at the smallest Qwen">Δℓ<sub>i,k</sub><sup>P</sup> = ℓ<sub>i,k</sub><sup>P</sup>(${last.size_b}B) − ℓ<sub>i,k</sub><sup>P</sup>(${first.size_b}B)</div><p class="simple-math-key"><strong>ℓ</strong> is one protein’s local alignment score (shared neighbors divided by k). <strong>Δℓ</strong> is its change from Qwen ${first.size_b}B to ${last.size_b}B: largest minus smallest. <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint; <strong>B</strong> after a size means billion parameters.<br>Δℓ &gt; 0 is a gain, Δℓ &lt; 0 is a loss, and Δℓ = 0 means the same endpoints.</p>${table(['PLM checkpoint','Gain','Same endpoints','Loss'],summaries,'simple-pattern-summary')}<p class="small-note">These counts include all ${N} proteins. The plots below use the same fixed panel.</p>${card(`Protein paths across Qwen sizes · ${model(state.a)}`,`<p class="simple-selection-note">Showing ${trajectoryRows.length} proteins: up to 40 from the fixed selection, plus the inspected protein if needed. <strong>Highlighted: ${esc(label(focus))}.</strong></p><div id="simple-trajectories" class="simple-plot" aria-label="Each protein's local alignment across all five Qwen sizes"></div><p class="simple-explainer">Click a line to inspect it; the tracked protein set stays fixed.</p><label class="simple-trajectory-picker">Highlight a gene / protein<select id="simple-trajectory-focus">${trajectoryRows.map(i=>`<option value="${i}" ${i===focus?'selected':''}>${esc(label(i))}</option>`).join('')}</select></label>`)}<details id="simple-table-details"><summary>Browse every protein’s scores across Qwen sizes</summary>${card(`What are each protein’s local kNN scores with ${shortModel(state.a)} across Qwen sizes?`,tableContent())}</details><div class="simple-actions"><button id="simple-examples">Choose 6 largest gains + 6 largest losses</button><button id="simple-select-page">Select this table page</button><button id="simple-select-matching">Select all matching proteins</button><button id="simple-export-all">All trajectories · CSV</button></div></section>
    <section class="simple-section" id="simple-question-2"><h2>Which proteins have higher local kNN alignment with ${esc(shortModel(state.a))} or ${esc(shortModel(state.b))}?</h2><p class="simple-selection-note"><strong>${selected.length} proteins selected.</strong> ${esc(state.selectionOrigin)} Changing a PLM or k keeps these exact protein IDs.</p><div class="simple-actions"><button id="simple-contrast-examples">Show 40 strongest PLM differences</button></div>${card('Difference heatmap across Qwen sizes','<div class="simple-heatmap-scroll" tabindex="0"><div id="simple-difference-heatmap" class="simple-plot" aria-label="Gene-labeled PLM difference heatmap across Qwen sizes"></div></div>')}<p class="simple-explainer">${selected.length>40?`Showing 40 of ${selected.length} selected proteins; all selected proteins enter global-change accounting.`:''} Click a row to inspect both underlying trajectories. On a narrow screen, scroll the heatmap horizontally.</p><div class="simple-equation" role="math" aria-label="Each difference heatmap cell equals local alignment with PLM A minus local alignment with PLM B, at the same Qwen size">D<sub>i,k</sub>(s) = ℓ<sub>i,k</sub><sup>A</sup>(s) − ℓ<sub>i,k</sub><sup>B</sup>(s)</div><p class="simple-math-key"><strong>D</strong> is one protein’s local alignment difference, A minus B. <strong>ℓ</strong> is its local alignment score (shared neighbors divided by k). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the neighbor count; <strong>s</strong> is Qwen size in billions of parameters.<br><strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. Positive D favors A; negative D favors B; zero means equal local scores.</p><details id="simple-delta-details"><summary>Compare endpoint gains and losses between PLMs</summary>${card('Change with one PLM versus change with the other',`<div id="simple-plm-change" class="simple-plot" aria-label="Same-protein changes compared between PLMs"></div><div class="simple-equation" role="math" aria-label="The difference between endpoint changes equals the change with A minus the change with B">Δℓ<sub>i,k</sub><sup>A</sup> − Δℓ<sub>i,k</sub><sup>B</sup></div><p class="simple-math-key"><strong>ℓ</strong> is one protein’s local alignment score (shared neighbors divided by k). <strong>Δℓ</strong> is its change from the Qwen ${first.size_b}B to ${last.size_b}B. <strong>i</strong> identifies the same protein in both terms; <strong>k = ${state.k}</strong> is the neighbor count; <strong>A = ${esc(shortModel(state.a))}</strong>; <strong>B = ${esc(shortModel(state.b))}</strong>. Positive means a larger gain (or smaller loss) with A.</p>`)}</details><details id="simple-heatmap-details"><summary>Selected proteins side by side across all Qwen sizes</summary><div class="simple-equation" role="math" aria-label="Each heatmap cell is the local alignment score: shared Qwen and PLM neighbors divided by k">ℓ<sub>i,k</sub><sup>P</sup>(s) = |𝒩<sub>k</sub><sup>Q(s)</sup>(i) ∩ 𝒩<sub>k</sub><sup>P</sup>(i)| / k</div><p class="simple-math-key"><strong>ℓ</strong> is this protein’s local alignment score: the fraction of its neighbors shared by Qwen and the PLM (0 to 1). <strong>i</strong> identifies the protein; <strong>k = ${state.k}</strong> is the number of neighbors; <strong>s</strong> is Qwen size in billions of parameters; <strong>P</strong> is the PLM checkpoint.<br><strong>Q(s)</strong> is Qwen at size s. <strong>𝒩</strong> is a set of nearest neighbors of i in the indicated model. <strong>∩</strong> means the shared members; <strong>|·|</strong> counts them. The protein itself is excluded. Superscripts identify models, not powers.</p><div class="simple-grid">${[state.a,state.b].map((p,j)=>card(model(p),`<div class="simple-heatmap-scroll" tabindex="0"><div id="simple-heat-${j}" class="simple-plot" aria-label="${esc(model(p))} selected protein trajectories"></div></div>`)).join('')}</div><p class="simple-explainer">Same proteins, same row order, same 0–1 color scale. ${selected.length>30?`Showing the first 30 of ${selected.length} selected proteins; all ${selected.length} enter global-change accounting.`:selected.length?'Each row follows one selected protein across all Qwen sizes.':'Choose proteins in the table to compare their rows.'} Click a heatmap cell to inspect that protein.</p></details>${card(`Inspect ${label(focus)}`,`<div id="simple-focus-plot" class="simple-plot" aria-label="Same protein across PLMs"></div>${table(['PLM checkpoint',...qs.map(q=>`${q.size_b}B`),'Local change','Endpoint pattern'],focusRows,'simple-focus-scores')}<p class="simple-explainer" id="simple-observation"><strong>Observed:</strong> ${esc(label(focus))} changes by ${signed(a.rows[focus].delta,2)} with ${esc(model(state.a))} and ${signed(b.rows[focus].delta,2)} with ${esc(model(state.b))}. The difference between those changes is ${signed(a.rows[focus].delta-b.rows[focus].delta,2)}.</p>`)}</section>
    <section class="simple-section" id="simple-question-3"><h2>How do protein changes add to the global kNN alignment change?</h2><p class="simple-explainer">For k = ${state.k}, the <strong>global kNN alignment score α<sub>k</sub><sup>P</sup>(s)</strong> is the mean of all ${N.toLocaleString()} local alignment scores between Qwen and PLM P. Each protein’s local-score change contributes <strong>Δℓ<sub>i,k</sub><sup>P</sup> / ${N.toLocaleString()}</strong> to the change in this global score.</p><div class="simple-equation" role="math" aria-label="The global kNN alignment score is the mean of all N local alignment scores">α<sub>k</sub><sup>P</sup>(s) = (1/N) ∑<sub>i=1</sub><sup>N</sup> ℓ<sub>i,k</sub><sup>P</sup>(s)</div><p class="simple-math-key"><strong>α</strong> is the global kNN alignment score. <strong>ℓ</strong> is one protein’s local alignment score (shared neighbors divided by k). <strong>N = ${N.toLocaleString()}</strong> is the total number of proteins in the fixed reference panel, including proteins hidden by display filters or outside the selection.<br><strong>i</strong> identifies a protein; <strong>∑</strong> means add the scores for all N proteins, then divide by N. <strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint; <strong>s</strong> is Qwen size in billions of parameters.</p><div class="simple-equation" role="math" aria-label="Protein i contributes its local change divided by N to the global alignment change">g<sub>i,k</sub><sup>P</sup> = Δℓ<sub>i,k</sub><sup>P</sup> / N</div><p class="simple-math-key"><strong>g</strong> is protein i’s additive contribution to the global alignment change. <strong>ℓ</strong> is its local score (shared neighbors divided by k); <strong>Δℓ</strong> is the local change from Qwen ${first.size_b}B to ${last.size_b}B.<br><strong>N = ${N.toLocaleString()}</strong> is the total protein count in the fixed reference panel, not the selected count. <strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint.</p><div class="simple-equation" role="math" aria-label="The global alignment endpoint change is the sum of all individual contributions">Δα<sub>k</sub><sup>P</sup> = α<sub>k</sub><sup>P</sup>(${last.size_b}B) − α<sub>k</sub><sup>P</sup>(${first.size_b}B) = ∑<sub>i=1</sub><sup>N</sup> g<sub>i,k</sub><sup>P</sup></div><p class="simple-math-key"><strong>Δα</strong> is the change in global kNN alignment from Qwen ${first.size_b}B to ${last.size_b}B; <strong>α</strong> is the mean local score over the full panel. <strong>g</strong> is one protein’s contribution to that change: its local-score change divided by <strong>N = ${N.toLocaleString()}</strong> total reference proteins.<br><strong>∑</strong> adds g over every protein i. <strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint; <strong>B</strong> after a size means billion parameters. Multiply g or Δα by 100 for percentage points. These are additive contributions, not deletion effects.</p>${card('Gains and losses across the whole panel','<div id="simple-account-plot" class="simple-plot" aria-label="Positive and negative contributions to the global alignment change"></div>')}${table(['PLM checkpoint','All gains (pp)','All losses (pp)','Net change (pp)','Selected net (pp)','Other proteins (pp)','Selected / net change'],accountRows,'simple-account-summary')}<p class="simple-explainer">Units: percentage points (pp) of the global kNN alignment score, k = ${state.k}. <strong>All gains + all losses = net change; selected net + other proteins = net change.</strong> A selection’s share of the net can be negative or exceed 100% when gains and losses offset each other.</p><div class="simple-equation" role="math" aria-label="The selected proteins' contribution is the sum of their individual contributions">G<sub>S,k</sub><sup>P</sup> = ∑<sub>i∈S</sub> g<sub>i,k</sub><sup>P</sup></div><p class="simple-math-key"><strong>S</strong> is the selected protein set; <strong>G<sub>S,k</sub><sup>P</sup></strong> is its combined contribution to the global score change. <strong>∑<sub>i∈S</sub></strong> means add over proteins i in S. Each <strong>g</strong> is that protein’s local-score change from Qwen ${first.size_b}B to ${last.size_b}B, divided by <strong>N = ${N.toLocaleString()}</strong> total reference proteins, not the number selected.<br><strong>k = ${state.k}</strong> is the neighbor count; <strong>P</strong> is the PLM checkpoint. <strong>Δα</strong> is the global score change across those endpoints. “Selected / net change” is 100 × G / Δα, a percentage shown only when Δα ≠ 0.</p><details id="simple-contribution-details"><summary>Individual contributions from the selected proteins</summary><p class="simple-comparison-key"><span><strong>${esc(model(state.a))}</strong></span><span>versus <strong>${esc(model(state.b))}</strong></span></p>${table(['Gene / protein name',`${shortModel(state.a)}: local change`,`${shortModel(state.a)}: contribution (pp)`,`${shortModel(state.b)}: local change`,`${shortModel(state.b)}: contribution (pp)`],selected.map(i=>[esc(label(i)),signed(a.rows[i].delta,2),signed(100*a.rows[i].contribution,5),signed(b.rows[i].delta,2),signed(100*b.rows[i].contribution,5)]),'simple-contributions')}</details><div class="simple-actions"><button id="simple-export-selected" ${selected.length?'':'disabled'}>Selected trajectories · CSV</button><button id="simple-export-account">Change accounting · CSV</button></div></section>
    <section class="simple-section simple-caveats" id="simple-followup"><h2>Interpretation and next questions</h2><p><strong>Observed patterns:</strong> these are recorded neighbor-overlap trajectories on the same proteins. Endpoint gains can include intermediate losses. They do not establish biological correctness or a causal effect of model size.</p><p><strong>Possible explanations to investigate:</strong> inspect Pfam families, homology, repeated names, and which neighbors change for selected proteins. Qwen’s generated text and thinking-token content also differ across sizes. These are hypotheses, not conclusions from the plots.</p><details><summary>Comparison settings and numerical accounting</summary><p>Counts are the unchanged evaluator’s shared-neighbor counts. Saved float32 global values differ by at most 0.00000005 in the validated endpoint comparison. Lines and points can overlap because neighbor counts are discrete; search the table to find a specific protein. Qwen sizes on trajectory plots are ordered categories; spacing is not proportional to parameter count.</p>${table(['Model','Recorded layer','Recorded pooling'],[...qs,...[state.a,state.b].map(p=>byPLM[p])].map(m=>[esc(m.label),esc(m.layer??m.selected_layer??'Not recorded'),esc(m.pooling??'Not recorded')]))}</details><p class="simple-explainer">Bootstrap, regression, deletion influence, family annotations, and the earlier metric views remain available as supporting analyses.</p><button data-report-mode="supporting">Open supporting analyses</button></section>`;
    const host=$('simple-view');
    host.insertAdjacentHTML('afterbegin',globalMarkup()+distributionMarkup()+topMarkup()+biologyMarkup()+redundancyMarkup());
    $('simple-question-2').querySelector('h2').insertAdjacentHTML('afterend',differentialMarkup());
    const scaling=document.createElement('div');scaling.id='guided-scaling';host.appendChild(scaling);scaling.append($('simple-question-1'));
    const accountingDetails=document.createElement('details');accountingDetails.id='guided-scaling-account-details';accountingDetails.innerHTML='<summary>Contributions to the global alignment change</summary>';accountingDetails.append($('simple-question-3'));scaling.append(accountingDetails);
    const routes={global:'guided-global',distribution:'guided-distribution',proteins:'guided-top-proteins',differential:'simple-question-2',biology:'guided-question-4',redundancy:'guided-question-5',scaling:'guided-scaling'};
    for(const [route,id] of Object.entries(routes))$(id).hidden=state.section!==route;
    document.querySelectorAll('[data-guided-section]').forEach(button=>button.setAttribute('aria-pressed',String(button.dataset.guidedSection===state.section)));
    $('simple-summary').hidden=state.section!=='global';$('simple-followup').hidden=true;
    document.querySelector('.simple-filter-controls').hidden=state.section!=='scaling';
    $('simple-b').closest('label').hidden=state.section==='proteins';
    $('simple-pattern').closest('label').hidden=state.section!=='scaling';$('simple-sort').closest('label').hidden=state.section!=='scaling';
    $('simple-qwen').value=state.qwen;$('simple-qwen').closest('label').hidden=['scaling','redundancy'].includes(state.section);
    document.querySelectorAll('#descriptive-report .simple-equation').forEach((equation,index)=>{const details=document.createElement('details');details.className='simple-formula-details';details.id='simple-formula-'+index;details.innerHTML='<summary>Formula and symbol definitions</summary>';equation.before(details);details.append(equation);const key=details.nextElementSibling;if(key?.classList.contains('simple-math-key'))details.append(key);});
    for(const id of openDetails)if($(id))$(id).open=true;
    bindTop();

    $('simple-trajectory-focus').onchange=e=>{state.focus=Number(e.target.value);render();};
    $('simple-prev').onclick=()=>{state.page--;render();};$('simple-next').onclick=()=>{state.page++;render();};
    $('simple-contrast-examples').onclick=()=>{chooseDifferenceExamples();render();};
    $('simple-examples').onclick=()=>{chooseExamples();render();};
    const manual=()=>{state.selectionOrigin='Proteins selected from the table; all remaining proteins stay in the comparison panel.';};
    $('simple-select-page').onclick=()=>{state.selected=new Set(filtered().slice(state.page*state.pageSize,(state.page+1)*state.pageSize));manual();render();};
    $('simple-select-matching').onclick=()=>{state.selected=new Set(filtered());manual();render();};
    document.querySelectorAll('[data-simple-select]').forEach(el=>el.onchange=()=>{const i=Number(el.dataset.simpleSelect);el.checked?state.selected.add(i):state.selected.delete(i);manual();render();});
    document.querySelectorAll('[data-simple-focus]').forEach(el=>el.onclick=()=>{state.focus=Number(el.dataset.simpleFocus);render();});
    $('simple-export-all').onclick=()=>download('descriptive_all_trajectories.csv',exportRows());$('simple-export-selected').onclick=()=>download('descriptive_selected_trajectories.csv',exportRows().filter(r=>r.selected));
    $('simple-export-account').onclick=()=>download('descriptive_change_accounting.csv',accounts.map(s=>({panel_id:D.panel.id,panel_n:N,k:state.k,plm_checkpoint:s.p,qwen_start:first.id,qwen_end:last.id,start_score:s.start,end_score:s.end,all_positive_change:s.positive,all_negative_change:s.negative,net_change:s.net,selected_n:selected.length,selected_positive_change:s.selectedPositive,selected_negative_change:s.selectedNegative,selected_net_change:s.selectedNet,other_proteins_change:s.remaining,selected_share_of_net:s.selectedShare})));
    document.querySelectorAll('[data-report-mode="supporting"]').forEach(el=>el.onclick=()=>showMode('supporting'));
    const shown=ix,marker={size:8,opacity:.6,symbol:shown.map(i=>state.selected.has(i)?'diamond':a.rows[i].deltaCount>0?'triangle-up':a.rows[i].deltaCount<0?'triangle-down':'circle'),color:shown.map(i=>state.selected.has(i)?'#8d3d91':a.rows[i].deltaCount>0?'#087b65':a.rows[i].deltaCount<0?'#bc512d':'#7c8e98')};
    const pointText=i=>`${esc(label(i))}<br>${esc(model(state.a))}: ${number(a.rows[i].values[0],2)} → ${number(a.rows[i].values.at(-1),2)} (${signed(a.rows[i].delta,2)})<br>${esc(model(state.b))}: ${number(b.rows[i].values[0],2)} → ${number(b.rows[i].values.at(-1),2)} (${signed(b.rows[i].delta,2)})`;
    const linkPlot=id=>{const el=$(id);el.on?.('plotly_click',e=>{const i=e.points?.[0]?.customdata;if(Number.isInteger(i)){state.focus=i;state.selected=new Set([i]);state.selectionOrigin='One protein selected by clicking a plot.';render();}});el.on?.('plotly_selected',e=>{const chosen=[...new Set((e?.points||[]).map(p=>p.customdata).filter(Number.isInteger))];if(chosen.length){state.selected=new Set(chosen);state.focus=chosen[0];state.selectionOrigin='Fixed proteins selected by a plot box/lasso.';window.DESCRIPTIVE.ready=false;setTimeout(render,0);}});};
    await plot('simple-trajectories',trajectoryRows.map(i=>({type:'scatter',mode:'lines+markers',x:qs.map(q=>q.size_b+'B'),y:a.rows[i].values,customdata:qs.map(()=>i),name:label(i),text:qs.map((q,j)=>`${esc(label(i))}<br>Qwen ${q.size_b}B · ${esc(model(state.a))}<br>Shared neighbors: ${a.rows[i].counts[j]}/${state.k}<br>Local alignment: ${number(a.rows[i].values[j],2)}<br>Endpoint change: ${signed(a.rows[i].delta,2)}`),hovertemplate:'%{text}<extra></extra>',line:{color:colorOf(state.a),width:i===focus?3.5:1.3,dash:a.rows[i].deltaCount>0?'solid':a.rows[i].deltaCount<0?'dash':'dot'},marker:{size:i===focus?8:4,symbol:i===focus?'diamond':'circle'},opacity:i===focus?1:.3})),{height:400,showlegend:false,hovermode:'closest',xaxis:{title:{text:'Qwen model size'},type:'category',categoryorder:'array',categoryarray:qs.map(q=>q.size_b+'B')},yaxis:{title:{text:'Local alignment<br>(shared neighbors ÷ k)'},range:[-.03,1.03],dtick:.2},margin:{l:78,r:20,t:12,b:74}});
    $('simple-trajectories').on?.('plotly_click',e=>{const i=e.points?.[0]?.customdata;if(Number.isInteger(i)){state.focus=i;window.DESCRIPTIVE.ready=false;setTimeout(render,0);}});
    await plot('simple-plm-change',[{type:'scatter',mode:'markers',x:shown.map(i=>b.rows[i].delta),y:shown.map(i=>a.rows[i].delta),customdata:shown,text:shown.map(pointText),hovertemplate:'%{text}<extra></extra>',marker,name:'Proteins'},{type:'scatter',mode:'lines',x:[-1,1],y:[-1,1],line:{color:'#778896',dash:'dot'},hoverinfo:'skip',name:'Same change with both PLMs'}],{height:380,showlegend:false,dragmode:'lasso',xaxis:{title:{text:`${esc(shortModel(state.b))}: local score change<br>(${last.size_b}B − ${first.size_b}B Qwen)`},range:[-1.03,1.03],zeroline:true,zerolinecolor:'#9babb5'},yaxis:{title:{text:`${esc(shortModel(state.a))}: local score change<br>(${last.size_b}B − ${first.size_b}B Qwen)`},range:[-1.03,1.03],zeroline:true,zerolinecolor:'#9babb5'},margin:{l:72,r:20,t:12,b:64}});linkPlot('simple-plm-change');
    for(const id of [...document.querySelectorAll('#simple-view details[id]')].map(node=>node.id)) {
      const node=$(id);let wasOpen=node.open;
      node.ontoggle=()=>{const open=node.open;if(open===wasOpen)return;wasOpen=open;if(open&&!active&&node.isConnected)render();};
    }
    if(differenceRows.length) {
      const height=Math.max(250,differenceRows.length*22+100);
      $('simple-difference-heatmap').style.height=height+'px';
      await plot('simple-difference-heatmap',[{type:'heatmap',x:qs.map(q=>q.size_b+'B'),y:differenceRows.map(i=>ids[i]),z:differenceRows.map(i=>a.rows[i].values.map((v,j)=>v-b.rows[i].values[j])),customdata:differenceRows.map(i=>qs.map(()=>i)),text:differenceRows.map(i=>qs.map((q,j)=>signed(a.rows[i].values[j]-b.rows[i].values[j],2))),hovertext:differenceRows.map(i=>qs.map((q,j)=>`${esc(label(i))}<br>Qwen ${q.size_b}B<br>${esc(model(state.a))}: ${number(a.rows[i].values[j],2)}<br>${esc(model(state.b))}: ${number(b.rows[i].values[j],2)}<br>${esc(shortModel(state.a))} − ${esc(shortModel(state.b))}: ${signed(a.rows[i].values[j]-b.rows[i].values[j],2)}`)),hovertemplate:'%{hovertext}<extra></extra>',zmin:-1,zmax:1,colorscale:[[0,'#b2182b'],[.25,'#ef8a62'],[.5,'#f7f7f7'],[.75,'#67a9cf'],[1,'#2166ac']],texttemplate:'%{text}',textfont:{size:10,color:'auto'},colorbar:{title:{text:'Local score<br>difference'},thickness:12}},...pairLegend().map((trace,j)=>({...trace,marker:{color:j?'#b2182b':'#2166ac',size:8}}))],{height:height+65,showlegend:true,legend:{orientation:'h',y:1.10,font:{size:10}},margin:{l:170,r:60,t:80,b:55},xaxis:{title:{text:'Qwen model size'}},yaxis:{title:{text:'Gene / protein name'},autorange:'reversed',tickvals:differenceRows.map(i=>ids[i]),ticktext:differenceRows.map(i=>esc(label(i)))}});
      $('simple-difference-heatmap').on?.('plotly_click',e=>{const i=e.points?.[0]?.customdata;if(Number.isInteger(i)){state.focus=i;window.DESCRIPTIVE.ready=false;setTimeout(render,0);}});
    } else $('simple-difference-heatmap').innerHTML='<p class="empty">Select proteins or choose the 40 strongest PLM differences.</p>';
    for(const [j,p] of [state.a,state.b].entries()) {
      if(!selectedRows.length){$('simple-heat-'+j).innerHTML='<p class="empty">Track proteins in Top proteins & families or select them in Scaling.</p>';continue;}
      const s=series(p),height=Math.max(240,selectedRows.length*23+100);$('simple-heat-'+j).style.height=height+'px';
      await plot('simple-heat-'+j,[{type:'heatmap',x:qs.map(q=>q.size_b+'B'),y:selectedRows.map(i=>ids[i]),z:selectedRows.map(i=>s.rows[i].values),customdata:selectedRows.map(i=>qs.map(()=>i)),text:selectedRows.map(i=>s.rows[i].values.map(v=>number(v,2))),hovertext:selectedRows.map(i=>qs.map((q,k)=>`${esc(label(i))}<br>${esc(model(p))}<br>Qwen ${q.size_b}B<br>Local alignment: ${s.rows[i].counts[k]}/${state.k} = ${number(s.rows[i].values[k],2)}`)),hovertemplate:'%{hovertext}<extra></extra>',zmin:0,zmax:1,colorscale:'Viridis',texttemplate:'%{text}',textfont:{size:10,color:'auto'},colorbar:{title:{text:'Local alignment<br>score ℓ'},thickness:10},showscale:true}],{height,margin:{l:155,r:50,t:15,b:52},xaxis:{title:{text:'Qwen model size'}},yaxis:{title:{text:'Gene / protein name'},tickvals:selectedRows.map(i=>ids[i]),ticktext:selectedRows.map(i=>esc(label(i))),autorange:'reversed',automargin:true},showlegend:false});
      $('simple-heat-'+j).on?.('plotly_click',e=>{const i=e.points?.[0]?.customdata;if(Number.isInteger(i)){state.focus=i;window.DESCRIPTIVE.ready=false;setTimeout(render,0);}});
    }
    await plot('simple-focus-plot',[state.a,state.b].map((p,j)=>({type:'scatter',mode:'lines+markers',x:qs.map(q=>q.size_b+'B'),y:series(p).rows[focus].values,name:model(p),line:{color:colorOf(p),width:2.5,dash:modelDash(p)},marker:{size:8,symbol:modelSymbol(p)},hovertemplate:`${esc(label(focus))}<br>${esc(model(p))}<br>Qwen %{x}: %{y:.3g}<extra></extra>`})),{height:320,xaxis:{title:{text:'Qwen model size'}},yaxis:{title:{text:'Local alignment<br>(shared neighbors ÷ k)'},range:[-.03,1.03],dtick:.2},margin:{l:78,r:20,t:48,b:74}});
    await plot('simple-account-plot',accounts.map((s,j)=>({type:'bar',x:['All gains','All losses','Net change'],y:[s.positive*100,s.negative*100,s.net*100],name:model(s.p),marker:{color:colorOf(s.p),pattern:{shape:j?'/':''}},hovertemplate:`${esc(model(s.p))}<br>%{x}: %{y:.3g} pp<extra></extra>`})),{height:320,barmode:'group',xaxis:{title:{text:'Contribution group'}},yaxis:{title:{text:`Global alignment change<br>(percentage points; k = ${state.k})`},zeroline:true,zerolinecolor:'#193243'},margin:{l:84,r:20,t:44,b:66}});
    if(state.section==='global')await drawGlobal();
    if(state.section==='distribution')await drawDistribution();
    if(state.section==='proteins'){await drawTop();bindTop();}
    if(state.section==='differential')await drawDifferential();
    if(state.section==='biology')await drawBiology();
    if(state.section==='redundancy')await drawRedundancy();

  }
  let active=false,again=false,promise=Promise.resolve();
  function render() {if(active){again=true;return promise;}active=true;window.DESCRIPTIVE.ready=false;promise=(async()=>{try{do{again=false;await renderNow();}while(again);}catch(e){window.DESCRIPTIVE.error=e.message;$('simple-view').innerHTML=`<p class="notice error">${esc(e.message)}</p>`;console.error(e);}finally{active=false;window.DESCRIPTIVE.ready=true;}})();return promise;}
  async function showMode(mode) {
    const fromCKA=!$('cka-report').hidden,fromSupporting=!$('supporting-report').hidden;
    if(fromSupporting&&mode!=='supporting') {
      const shared=window.DASHBOARD.getState();
      for(const key of ['qwen','a','b','k'])state[key]=shared[key];
      state.selected=new Set(shared.selected.map(id=>idIndex.get(id)).filter(i=>i!==undefined));
      state.focus=idIndex.get(shared.focus)??state.focus;
      state.selectionOrigin='Exact protein IDs carried from All analyses.';
    }
    for(const name of ['descriptive','cka','supporting']) {
      $(name+'-report').hidden=name!==mode;
      $('show-'+name).setAttribute('aria-pressed',String(name===mode));
    }
    document.body.classList.toggle('descriptive-mode',mode==='descriptive');
    document.body.classList.toggle('cka-mode',mode==='cka');
    if(mode==='descriptive')await render();
    else if(mode==='cka')await window.CKA_REPORT.render(window.DESCRIPTIVE.getState());
    else if(!fromSupporting)await window.DASHBOARD.update({qwen:state.qwen,a:state.a,b:state.b,k:state.k,metric:fromCKA?'cka':'knn',value:fromCKA?'contribution':'local',selected:[...state.selected].map(i=>ids[i]),focus:ids[state.focus]});
  }
  window.DESCRIPTIVE={ready:false,getState:()=>({...state,selected:[...state.selected].map(i=>ids[i]),focus:ids[state.focus]}),getSeries:p=>series(p),getAccounting:p=>accounting(p),getRows:exportRows,showMode,update:async patch=>{for(const [key,value] of Object.entries(patch)){if(key==='selected')state.selected=new Set(value.map(id=>idIndex.get(id)).filter(i=>i!==undefined));else if(key==='focus')state.focus=idIndex.get(value)??0;else if(Object.hasOwn(state,key))state[key]=value;}return render();}};
  for(const id of ['simple-a','simple-b'])$(id).innerHTML=plms.map(p=>`<option value="${esc(p.id)}">${esc(p.label)}</option>`).join('');
  $('simple-qwen').innerHTML=qs.map(q=>`<option value="${q.id}">Qwen ${q.size_b}B</option>`).join('');
  $('simple-qwen').onchange=e=>{state.qwen=e.target.value;state.proteinPage=0;render();};
  document.querySelectorAll('[data-guided-section]').forEach(button=>button.onclick=()=>{state.section=button.dataset.guidedSection;render();});
  $('simple-k').innerHTML=D.settings.knn_ks.map(k=>`<option value="${k}">${k}</option>`).join('');
  for(const [id,key] of [['simple-a','a'],['simple-b','b'],['simple-k','k'],['simple-pattern','pattern'],['simple-sort','sort']])$(id).onchange=e=>{state[key]=key==='k'?Number(e.target.value):e.target.value;state.page=0;state.proteinPage=0;render();};
  let searchTimer;$('simple-search').oninput=e=>{state.search=e.target.value;state.page=0;state.proteinPage=0;clearTimeout(searchTimer);searchTimer=setTimeout(render,180);};
  $('simple-clear').onclick=()=>{state.selected.clear();state.selectionOrigin='Selection cleared.';render();};
  for(const name of ['descriptive','cka','supporting'])$('show-'+name).onclick=()=>showMode(name);
  let resizeTimer;window.addEventListener('resize',()=>{clearTimeout(resizeTimer);if(!$('descriptive-report').hidden)resizeTimer=setTimeout(render,220);else if(!$('cka-report').hidden)resizeTimer=setTimeout(()=>window.CKA_REPORT.render(window.DESCRIPTIVE.getState()),220);});
  const initialTop=topProteinData();state.selected=new Set(initialTop.rows.map(row=>row.i));state.focus=initialTop.rows[0].i;state.selectionOrigin=`Top ${initialTop.requested} plus score ties for Qwen ${initialTop.q.size_b}B and ${model(state.a)}, k=${state.k}.`;render();
})();
