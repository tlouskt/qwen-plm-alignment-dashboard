/* Static hosting loader; all analysis code and measurements are unchanged. */
(() => {
  'use strict';
  const expectedSha = 'c75cc4ddf34ccd285de6cbdf9d9becb1c11bc9b5fac31031a453e73d1b8b4f00';
  const expectedCompressedBytes = 22591578;
  const shell = document.getElementById('dashboard-shell');
  const notice = document.getElementById('hosted-loading');
  const message = document.getElementById('hosted-loading-message');
  const progress = document.getElementById('hosted-loading-progress');
  const retry = document.getElementById('hosted-retry');
  const status = window.HOSTED_DASHBOARD = {ready: false, phase: 'loading', expectedDataSha256: expectedSha, dataSha256: null, decoder: null, error: null};
  function say(text, phase) { message.textContent = text; status.phase = phase; }
  function script(path) {
    return new Promise((resolve, reject) => {
      const element = document.createElement('script');
      element.src = path;
      element.onload = resolve;
      element.onerror = () => reject(new Error('A bundled script could not be downloaded. Please check your connection and try again.'));
      document.head.appendChild(element);
    });
  }
  async function downloadData() {
    const response = await fetch('assets/dashboard_data.json.gz');
    if (!response.ok) throw new Error('The bundled data could not be downloaded (HTTP ' + response.status + ').');
    if (!response.body || !response.body.getReader) return new Uint8Array(await response.arrayBuffer());
    const reader = response.body.getReader(), chunks = [];
    let received = 0;
    while (true) {
      const part = await reader.read();
      if (part.done) break;
      chunks.push(part.value); received += part.value.length;
      progress.value = Math.min(received, expectedCompressedBytes);
      message.textContent = 'Downloading the fixed-panel data… ' + (received / 1048576).toFixed(1) + ' MB';
    }
    const bytes = new Uint8Array(received);
    let offset = 0;
    for (const chunk of chunks) { bytes.set(chunk, offset); offset += chunk.length; }
    return bytes;
  }
  async function decompress(compressed) {
    // Some static servers transparently decode files tagged Content-Encoding: gzip.
    if (compressed[0] !== 0x1f || compressed[1] !== 0x8b) { status.decoder = 'server'; return compressed; }
    if (typeof DecompressionStream === 'function') {
      try {
        const stream = new Blob([compressed]).stream().pipeThrough(new DecompressionStream('gzip'));
        const bytes = new Uint8Array(await new Response(stream).arrayBuffer());
        status.decoder = 'native'; return bytes;
      } catch (_) { /* A locally bundled decoder supports browsers without native gzip. */ }
    }
    await script('assets/pako_inflate.min.js');
    status.decoder = 'pako';
    return window.pako.ungzip(compressed);
  }
  async function loadData() {
    const compressed = await downloadData();
    progress.removeAttribute('value');
    say('Unpacking the saved measurements…', 'unpacking');
    const bytes = await decompress(compressed);
    if (window.crypto && window.crypto.subtle) {
      const digest = await window.crypto.subtle.digest('SHA-256', bytes);
      status.dataSha256 = Array.from(new Uint8Array(digest), value => value.toString(16).padStart(2, '0')).join('');
      if (status.dataSha256 !== expectedSha) throw new Error('The downloaded data failed its integrity check. Please reload to download it again.');
    }
    say('Reading the saved measurements…', 'parsing');
    return JSON.parse(new TextDecoder('utf-8', {fatal: true}).decode(bytes));
  }
  async function ready() {
    const start = performance.now();
    while (!window.DESCRIPTIVE || !window.DESCRIPTIVE.ready) {
      if (performance.now() - start > 120000) throw new Error('The initial view did not finish rendering. Please reload, or try a current browser.');
      await new Promise(resolve => setTimeout(resolve, 50));
    }
    if (window.DESCRIPTIVE.error) throw new Error(window.DESCRIPTIVE.error);
    if (window.DASHBOARD && window.DASHBOARD.lastError) throw new Error(window.DASHBOARD.lastError);
  }
  retry.addEventListener('click', () => location.reload());
  async function boot() {
    try {
      if (location.protocol === 'file:') throw new Error('Open this version through its website address. The original single-file report is the version for opening directly from disk.');
      const [data] = await Promise.all([loadData(), script('assets/plotly.min.js')]);
      window.ALIGNMENT_DATA = data;
      say('Drawing the report…', 'rendering');
      shell.hidden = false;
      shell.style.visibility = 'hidden';
      await script('assets/dashboard_ui.js');
      await script('assets/dashboard_cka.js');
      await script('assets/dashboard_descriptive.js');
      await ready();
      shell.style.removeProperty('visibility');
      shell.removeAttribute('aria-busy');
      notice.hidden = true;
      status.phase = 'ready'; status.ready = true;
      document.dispatchEvent(new Event('alignment-hosted-ready'));
    } catch (error) {
      shell.hidden = true;
      progress.hidden = true;
      retry.hidden = false;
      notice.setAttribute('role', 'alert');
      status.error = error instanceof Error ? error.message : String(error);
      say(status.error, 'error');
      console.error('Dashboard load failed:', error);
    }
  }
  boot();
})();
