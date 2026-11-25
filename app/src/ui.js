/* BondX Studio v2 — user interface layer. Depends on window.BondX (engine) and window.BondXFF (force-field data). */
(function () {
  'use strict';
  const B = window.BondX;
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));
  const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[c]);
  const f2 = (v, n = 2) => (v === null || v === undefined || isNaN(v)) ? '—' : Number(v).toFixed(n);
  const elTag = e => `<span class="el ${esc(e)}">${esc(e === 'CL' ? 'Cl' : e === 'BR' ? 'Br' : e)}</span>`;
  const HCOLOR = { I: 'var(--io)', BR: 'var(--br)', CL: 'var(--cl)', F: 'var(--s4)' };
  const HSHAPE = { I: 'circle', BR: 'square', CL: 'triangle', F: 'diamond' };
  const HN = x => x === 'BR' ? 'Br' : x === 'CL' ? 'Cl' : x;
  const DEF = B.defaults();

  /* ------------------------------ state ------------------------------ */
  const ST = {
    runId: 'BXR-' + new Date().toISOString().slice(0, 10).replace(/-/g, '') + '-' + Math.random().toString(36).slice(2, 7).toUpperCase(),
    S0: null, S: null, name: '', srcText: '', srcPQR: false, hits: [], P: B.defaults(),
    traj: null, batch: [], mat: null, cands: [], models: {}, audit: [], mol2: null, mol2Name: ''
  };
  window.BondXState = ST;
  function log(action, detail) { ST.audit.push({ t: new Date().toISOString(), action, detail: detail || '' }); }
  async function sha256(data) {
    try { const buf = typeof data === 'string' ? new TextEncoder().encode(data) : new Uint8Array(data); const h = await crypto.subtle.digest('SHA-256', buf);
      return Array.from(new Uint8Array(h)).map(b => b.toString(16).padStart(2, '0')).join(''); }
    catch (e) { return 'fnv:' + B.hashString(typeof data === 'string' ? data : String(data.byteLength)); }
  }
  function toast(msg, ms = 2800) { const t = $('#toast'); t.textContent = msg; t.classList.add('on'); clearTimeout(t._h); t._h = setTimeout(() => t.classList.remove('on'), ms); }
  function download(name, content, type = 'text/plain') {
    const blob = content instanceof Blob ? content : new Blob([content], { type }); const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }
  const readText = f => new Promise((res, rej) => { const r = new FileReader(); r.onload = () => res(r.result); r.onerror = rej; r.readAsText(f); });
  const readBuf = f => new Promise((res, rej) => { const r = new FileReader(); r.onload = () => res(r.result); r.onerror = rej; r.readAsArrayBuffer(f); });
  function dropZone(el, onFiles, multiple) {
    const inp = el.querySelector('input'); if (multiple) inp.multiple = true;
    el.addEventListener('click', e => { if (e.target !== inp) inp.click(); });
    inp.addEventListener('change', () => { if (inp.files.length) onFiles([...inp.files]); inp.value = ''; });
    el.addEventListener('dragover', e => { e.preventDefault(); el.classList.add('over'); });
    el.addEventListener('dragleave', () => el.classList.remove('over'));
    el.addEventListener('drop', e => { e.preventDefault(); el.classList.remove('over'); if (e.dataTransfer.files.length) onFiles([...e.dataTransfer.files]); });
  }
  const mol2For = S => B.mol2ForStructure(ST.mol2, S);
  const ffAssign = S => B.assignForceField(S, ST.P, { mol2: mol2For(S) });

  /* ------------------------------ ZIP (store) ------------------------------ */
  const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
  function crc32(u8) { let c = 0xFFFFFFFF; for (let i = 0; i < u8.length; i++) c = CRC[(c ^ u8[i]) & 0xFF] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; }
  function zip(files) {
    const enc = new TextEncoder(); const parts = []; const central = []; let off = 0;
    for (const [name, content] of Object.entries(files)) {
      const data = enc.encode(content); const nm = enc.encode(name); const crc = crc32(data);
      const lh = new DataView(new ArrayBuffer(30)); lh.setUint32(0, 0x04034b50, true); lh.setUint16(4, 20, true); lh.setUint16(8, 0, true);
      lh.setUint32(14, crc, true); lh.setUint32(18, data.length, true); lh.setUint32(22, data.length, true); lh.setUint16(26, nm.length, true);
      parts.push(new Uint8Array(lh.buffer), nm, data);
      const ch = new DataView(new ArrayBuffer(46)); ch.setUint32(0, 0x02014b50, true); ch.setUint16(4, 20, true); ch.setUint16(6, 20, true);
      ch.setUint32(16, crc, true); ch.setUint32(20, data.length, true); ch.setUint32(24, data.length, true); ch.setUint16(28, nm.length, true); ch.setUint32(42, off, true);
      central.push(new Uint8Array(ch.buffer), nm); off += 30 + nm.length + data.length;
    }
    const csize = central.reduce((s, p) => s + p.length, 0); const end = new DataView(new ArrayBuffer(22));
    end.setUint32(0, 0x06054b50, true); end.setUint16(8, Object.keys(files).length, true); end.setUint16(10, Object.keys(files).length, true); end.setUint32(12, csize, true); end.setUint32(16, off, true);
    return new Blob([...parts, ...central, new Uint8Array(end.buffer)], { type: 'application/zip' });
  }

  /* ------------------------------ tooltip ------------------------------ */
  const tip = $('#tip');
  function showTip(html, e) { tip.innerHTML = html; tip.style.display = 'block'; const x = e.clientX + 14, y = e.clientY + 14; const w = tip.offsetWidth, hh = tip.offsetHeight;
    tip.style.left = Math.min(x, innerWidth - w - 8) + 'px'; tip.style.top = Math.min(y, innerHeight - hh - 8) + 'px'; }
  function hideTip() { tip.style.display = 'none'; }
  document.addEventListener('mousemove', e => { const t = e.target.closest && e.target.closest('[data-tip]'); if (t) showTip(t.getAttribute('data-tip'), e); else if (!e.target.closest || !e.target.closest('.xhair')) hideTip(); });

  /* ------------------------------ charts (SVG) ------------------------------ */
  const Charts = (() => {
    const nice = (lo, hi, n = 5) => { if (!(hi - lo > 1e-9)) { const c = (hi + lo) / 2; lo = c - Math.max(0.5, Math.abs(c) * 0.1); hi = c + Math.max(0.5, Math.abs(c) * 0.1); } const span = hi - lo; const step0 = span / n; const mag = Math.pow(10, Math.floor(Math.log10(step0))); const st = [1, 2, 2.5, 5, 10].map(s => s * mag).find(s => s >= step0);
      const a = Math.floor(lo / st) * st, b = Math.ceil(hi / st) * st; const ticks = []; for (let v = a; v <= b + st / 2; v += st) ticks.push(+v.toFixed(10)); return { a, b, ticks }; };
    const marker = (shape, x, y, r, fill, tipHtml) => {
      const t = tipHtml ? ` data-tip="${esc(tipHtml)}"` : ''; const ring = `stroke="#fff" stroke-width="2"`;
      if (shape === 'square') return `<rect x="${x - r}" y="${y - r}" width="${2 * r}" height="${2 * r}" rx="1.5" fill="${fill}" ${ring}${t}/>`;
      if (shape === 'triangle') return `<path d="M${x},${y - r * 1.2} L${x + r * 1.1},${y + r * .8} L${x - r * 1.1},${y + r * .8} Z" fill="${fill}" ${ring}${t}/>`;
      if (shape === 'diamond') return `<path d="M${x},${y - r * 1.2} L${x + r * 1.2},${y} L${x},${y + r * 1.2} L${x - r * 1.2},${y} Z" fill="${fill}" ${ring}${t}/>`;
      return `<circle cx="${x}" cy="${y}" r="${r}" fill="${fill}" ${ring}${t}/>`;
    };
    function frame(W, H, m, xs, ys, xLabel, yLabel, fmtX = v => v, fmtY = v => v) {
      let s = '';
      ys.ticks.forEach(v => { const y = m.t + (H - m.t - m.b) * (1 - (v - ys.a) / (ys.b - ys.a)); s += `<line class="gr" x1="${m.l}" x2="${W - m.r}" y1="${y}" y2="${y}"/><text x="${m.l - 6}" y="${y + 3.5}" text-anchor="end">${fmtY(v)}</text>`; });
      xs.ticks.forEach(v => { const x = m.l + (W - m.l - m.r) * (v - xs.a) / (xs.b - xs.a); s += `<text x="${x}" y="${H - m.b + 15}" text-anchor="middle">${fmtX(v)}</text>`; });
      s += `<line class="ax" x1="${m.l}" x2="${W - m.r}" y1="${H - m.b}" y2="${H - m.b}"/>`;
      if (xLabel) s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle" style="font-weight:600">${esc(xLabel)}</text>`;
      if (yLabel) s += `<text transform="translate(12,${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle" style="font-weight:600">${esc(yLabel)}</text>`;
      return s;
    }
    function legendHTML(items) { return `<div class="legend">${items.map(i => `<span><i style="background:${i.color}"></i>${esc(i.name)}</span>`).join('')}</div>`; }
    function scatter(el, o) {
      const W = Math.max(360, el.clientWidth || 620), H = o.height || 320, m = { l: 52, r: 16, t: 14, b: 40 };
      if (!o.points.length) { el.innerHTML = '<p class="muted">No data.</p>'; return; }
      const xv = o.points.map(p => p.x), yv = o.points.map(p => p.y);
      const pad = (a, b) => { const r = (b - a) || Math.max(Math.abs(a) * 0.2, 1); return [a - r * 0.08, b + r * 0.08]; };
      const xd = o.xDomain || pad(Math.min(...xv), Math.max(...xv)), yd = o.yDomain || pad(Math.min(...yv), Math.max(...yv));
      const xs = nice(xd[0], xd[1]), ys = nice(yd[0], yd[1]);
      const X = v => m.l + (W - m.l - m.r) * (v - xs.a) / (xs.b - xs.a), Y = v => m.t + (H - m.t - m.b) * (1 - (v - ys.a) / (ys.b - ys.a));
      let s = frame(W, H, m, xs, ys, o.xLabel, o.yLabel, o.fmtX, o.fmtY);
      (o.refLines || []).filter(r => (r.y === undefined || (r.y >= ys.a && r.y <= ys.b)) && (r.x === undefined || (r.x >= xs.a && r.x <= xs.b))).forEach(r => {
        if (r.y !== undefined) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(r.y)}" y2="${Y(r.y)}" stroke="var(--ink-3)" stroke-dasharray="4 4"/><text x="${W - m.r}" y="${Y(r.y) - 4}" text-anchor="end">${esc(r.label)}</text>`;
        if (r.x !== undefined) s += `<line y1="${m.t}" y2="${H - m.b}" x1="${X(r.x)}" x2="${X(r.x)}" stroke="var(--ink-3)" stroke-dasharray="4 4"/><text x="${X(r.x) + 4}" y="${m.t + 10}">${esc(r.label)}</text>`; });
      o.points.forEach(p => { const g = o.groups.find(g => g.key === p.group) || o.groups[0]; s += marker(g.shape || 'circle', X(p.x), Y(p.y), 5, g.color, p.tip); });
      el.innerHTML = `<div class="chart">${o.groups.length > 1 ? legendHTML(o.groups) : ''}<svg viewBox="0 0 ${W} ${H}">${s}</svg></div>`;
    }
    function hbars(el, o) {
      const rows = o.rows; if (!rows.length) { el.innerHTML = '<p class="muted">No data.</p>'; return; }
      const bh = 18, gap = 8, m = { l: o.labelWidth || 210, r: 70, t: o.threshold !== undefined ? 22 : 8, b: 34 }; const W = Math.max(360, el.clientWidth || 620); const H = m.t + m.b + rows.length * (bh + gap);
      const max = o.max ?? (Math.max(...rows.map(r => r.value), o.threshold || 0) * 1.08 || 1); const xs = nice(0, max, 5);
      const X = v => m.l + (W - m.l - m.r) * v / xs.b;
      let s = '';
      xs.ticks.forEach(v => { s += `<line class="gr" y1="${m.t}" y2="${H - m.b}" x1="${X(v)}" x2="${X(v)}"/><text x="${X(v)}" y="${H - m.b + 15}" text-anchor="middle">${o.fmt ? o.fmt(v) : v}</text>`; });
      rows.forEach((r, i) => { const y = m.t + i * (bh + gap); const w = Math.max(0, X(r.value) - m.l);
        s += `<text x="${m.l - 8}" y="${y + bh / 2 + 4}" text-anchor="end">${esc(r.label.length > 34 ? r.label.slice(0, 33) + '…' : r.label)}</text>`;
        s += `<rect x="${m.l}" y="${y}" width="${Math.max(w, 1)}" height="${bh}" rx="4" fill="${r.color || 'var(--s1)'}" data-tip="${esc(r.tip || (r.label + ': ' + f2(r.value)))}"/>`;
        s += `<text x="${m.l + w + 6}" y="${y + bh / 2 + 4}" style="font-family:var(--mono)">${esc(r.text ?? f2(r.value))}</text>`; });
      if (o.threshold !== undefined) s += `<line x1="${X(o.threshold)}" x2="${X(o.threshold)}" y1="${m.t - 4}" y2="${H - m.b}" stroke="var(--bad)" stroke-dasharray="4 3"/><text x="${X(o.threshold) + 4}" y="${m.t - 8}" style="fill:var(--bad)">${esc(o.thresholdLabel || 'threshold')}</text>`;
      s += `<line class="ax" x1="${m.l}" x2="${m.l}" y1="${m.t - 4}" y2="${H - m.b}"/>`;
      if (o.xLabel) s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle" style="font-weight:600">${esc(o.xLabel)}</text>`;
      el.innerHTML = `<div class="chart">${o.legend ? legendHTML(o.legend) : ''}<svg viewBox="0 0 ${W} ${H}">${s}</svg></div>`;
    }
    function line(el, o) {
      const W = Math.max(360, el.clientWidth || 620), H = o.height || 220, m = { l: 52, r: 16, t: 12, b: 38 };
      const all = o.series.flatMap(s => s.ys); if (!all.length) { el.innerHTML = '<p class="muted">No data.</p>'; return; }
      const xsAll = o.series.flatMap(s => s.xs);
      const xs = nice(o.xDomain ? o.xDomain[0] : Math.min(...xsAll), o.xDomain ? o.xDomain[1] : Math.max(...xsAll)), ys = nice(Math.min(...all, o.refY ?? Infinity), Math.max(...all, o.refY ?? -Infinity));
      const X = v => m.l + (W - m.l - m.r) * (v - xs.a) / (xs.b - xs.a), Y = v => m.t + (H - m.t - m.b) * (1 - (v - ys.a) / (ys.b - ys.a));
      let s = frame(W, H, m, xs, ys, o.xLabel, o.yLabel, o.fmtX, o.fmtY);
      if (o.refY !== undefined) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(o.refY)}" y2="${Y(o.refY)}" stroke="var(--bad)" stroke-dasharray="4 3"/><text x="${W - m.r}" y="${Y(o.refY) - 4}" text-anchor="end" style="fill:var(--bad)">${esc(o.refLabel || '')}</text>`;
      o.series.forEach(se => {
        let d = '', prev = null; se.xs.forEach((x, k) => { const gapped = prev !== null && (x - prev) > (o.gapAbove || 1e9); d += (k === 0 || gapped ? 'M' : 'L') + X(x).toFixed(1) + ',' + Y(se.ys[k]).toFixed(1); prev = x; });
        s += `<path d="${d}" fill="none" stroke="${se.color}" stroke-width="2" stroke-linejoin="round"/>`;
        if (se.xs.length < 80) se.xs.forEach((x, k) => { s += `<circle cx="${X(x)}" cy="${Y(se.ys[k])}" r="2.6" fill="${se.color}"/>`; });
      });
      s += `<g class="xh"><line class="xl" y1="${m.t}" y2="${H - m.b}" stroke="var(--ink-3)" stroke-width="1" style="display:none"/></g><rect class="xhair" x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}" fill="transparent"/>`;
      el.innerHTML = `<div class="chart">${o.series.length > 1 ? legendHTML(o.series) : ''}<svg viewBox="0 0 ${W} ${H}">${s}</svg></div>`;
      const svg = el.querySelector('svg'), rect = el.querySelector('.xhair'), xl = el.querySelector('.xl');
      rect.addEventListener('mousemove', e => {
        const bb = svg.getBoundingClientRect(); const vx = (e.clientX - bb.left) * W / bb.width; const xv = xs.a + (vx - m.l) / (W - m.l - m.r) * (xs.b - xs.a);
        let html = ''; let snapX = null;
        o.series.forEach(se => { let bi = -1, bd = Infinity; se.xs.forEach((x, k) => { const dd = Math.abs(x - xv); if (dd < bd) { bd = dd; bi = k; } });
          if (bi >= 0) { snapX = se.xs[bi]; html += `<div><span style="color:${se.color}">●</span> ${esc(se.name)}: <b>${f2(se.ys[bi], o.digits ?? 3)}</b></div>`; } });
        if (snapX === null) return; xl.style.display = ''; xl.setAttribute('x1', X(snapX)); xl.setAttribute('x2', X(snapX));
        showTip(`<div>${esc(o.xLabel || 'x')} ${f2(snapX, 0)}</div>${html}`, e);
      });
      rect.addEventListener('mouseleave', () => { xl.style.display = 'none'; hideTip(); });
    }
    function diverging(el, o) {
      const rows = o.rows; if (!rows.length) { el.innerHTML = '<p class="muted">No data.</p>'; return; }
      const bh = 22, gap = 12, m = { l: o.labelWidth || 200, r: 24, t: 8, b: 36 }; const W = Math.max(360, el.clientWidth || 660); const H = m.t + m.b + rows.length * (bh + gap);
      let lo = 0, hi = 0; rows.forEach(r => { let n = 0, p = 0; o.comps.forEach(c => { const v = r.v[c.key]; if (v === undefined || v === null) return; if (v < 0) n += v; else p += v; }); lo = Math.min(lo, n); hi = Math.max(hi, p); });
      const xs = nice(lo, hi, 6); const X = v => m.l + (W - m.l - m.r) * (v - xs.a) / (xs.b - xs.a);
      let s = '';
      xs.ticks.forEach(v => { s += `<line class="gr" y1="${m.t}" y2="${H - m.b}" x1="${X(v)}" x2="${X(v)}"/><text x="${X(v)}" y="${H - m.b + 15}" text-anchor="middle">${v}</text>`; });
      rows.forEach((r, i) => {
        const y = m.t + i * (bh + gap); let neg = 0, pos = 0;
        s += `<text x="${m.l - 8}" y="${y + bh / 2 + 4}" text-anchor="end">${esc(r.label.length > 30 ? r.label.slice(0, 29) + '…' : r.label)}</text>`;
        o.comps.forEach(c => { const v = r.v[c.key]; if (v === undefined || v === null || v === 0) return;
          let x0, x1; if (v < 0) { x1 = X(neg); neg += v; x0 = X(neg); } else { x0 = X(pos); pos += v; x1 = X(pos); }
          s += `<rect x="${x0 + 1}" y="${y}" width="${Math.max(0, x1 - x0 - 2)}" height="${bh}" rx="3" fill="${c.color}" data-tip="${esc(`<b>${r.label}</b><br>${c.name}: ${f2(v)} kJ/mol`)}"/>`; });
        if (r.total !== undefined) s += `<path d="M${X(r.total)},${y - 3} l4,-5 h-8 z" fill="var(--ink)" data-tip="${esc('Total interaction: ' + f2(r.total) + ' kJ/mol')}"/>`;
      });
      s += `<line x1="${X(0)}" x2="${X(0)}" y1="${m.t - 4}" y2="${H - m.b}" stroke="var(--ink)" stroke-width="1.2"/>`;
      s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle" style="font-weight:600">${esc(o.xLabel || 'kJ/mol')}</text>`;
      el.innerHTML = `<div class="chart">${legendHTML(o.comps)}<svg viewBox="0 0 ${W} ${H}">${s}</svg></div>`;
    }
    return { scatter, hbars, line, diverging };
  })();

  /* ------------------------------ 3D viewer ------------------------------ */
  const V = {
    viewer: null, ready: false, host: null,
    load() {
      return new Promise(res => {
        if (window.$3Dmol) { this.ready = true; return res(true); }
        const urls = ['https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js', 'https://3Dmol.org/build/3Dmol-min.js'];
        const tryLoad = i => { if (i >= urls.length) return res(false); const s = document.createElement('script'); s.src = urls[i]; s.onload = () => { this.ready = !!window.$3Dmol; res(this.ready); }; s.onerror = () => tryLoad(i + 1); document.head.appendChild(s); };
        tryLoad(0);
      });
    },
    ensure() {
      if (!this.host) { this.host = document.createElement('div'); this.host.id = 'viewerHost'; this.host.innerHTML = '<div class="vmsg">Load a structure to view it in 3D.</div>'; }
      if (this.ready && !this.viewer) {
        this.host.innerHTML = '<div class="vtools"><button class="btn sm" id="vReset">Reset view</button><button class="btn sm" id="vPng">PNG</button></div>';
        this.viewer = $3Dmol.createViewer(this.host, { backgroundColor: 'white', antialias: true });
        this.host.querySelector('#vReset').onclick = () => { this.viewer.zoomTo(); this.viewer.render(); };
        this.host.querySelector('#vPng').onclick = () => { const uri = this.viewer.pngURI(); const a = document.createElement('a'); a.href = uri; a.download = (ST.name || 'BondX') + '_view.png'; a.click(); };
      }
      return this.host;
    },
    mount(slot) { const h = this.ensure(); if (slot && h.parentNode !== slot) slot.appendChild(h); if (this.viewer) setTimeout(() => { this.viewer.resize(); this.viewer.render(); }, 30); },
    offline() { if (!this.ready) { this.ensure().innerHTML = '<div class="vmsg">3D viewer library (3Dmol.js) could not be loaded — you appear to be offline.<br>All calculations still work; only the 3D view is unavailable.</div>'; return true; } return false; },
    structure(S, frameIdx, hits, focus) {
      if (this.offline()) return; this.ensure(); const v = this.viewer; v.clear();
      v.addModel(B.writePDB(S, S.frames[frameIdx || 0]), 'pdb');
      const chains = [...new Set(S.atoms.filter(a => !a.het).map(a => a.chain))]; const cc = ['#8fb8de', '#f2c38b', '#9ed6b6', '#d7b3e6', '#f2a7a7', '#b8c0cc'];
      chains.forEach((c, i) => v.setStyle({ chain: c === '_' ? ' ' : c, hetflag: false }, { cartoon: { color: cc[i % cc.length], opacity: 0.92 } }));
      v.setStyle({ resn: ['HOH', 'WAT', 'TIP3', 'SOL'] }, { sphere: { radius: 0.25, color: '#e2483d', opacity: 0.6 } });
      S.residues.filter(r => r.cls === 'ion').forEach(r => v.setStyle({ resi: r.resSeq, resn: r.resName }, { sphere: { radius: 0.7, color: '#2ca02c' } }));
      S.residues.filter(r => r.cls === 'ligand').forEach(r => v.setStyle({ resi: r.resSeq, resn: r.resName, chain: r.chain === '_' ? ' ' : r.chain }, { stick: { radius: 0.18, colorscheme: 'lightGreyCarbon' }, sphere: { scale: 0.22, colorscheme: 'lightGreyCarbon' } }));
      if (hits && hits.length) {
        const frame = S.frames[frameIdx || 0]; const accRes = new Set(hits.map(h => h.accRes));
        accRes.forEach(r => { const R = S.residues[r]; v.addStyle({ resi: R.resSeq, chain: R.chain === '_' ? ' ' : R.chain, resn: R.resName }, { stick: { radius: 0.14, colorscheme: 'greenCarbon' } }); });
        hits.forEach(h => {
          const px = B.coordOf(frame, h.X); let py;
          if (h.type === 'pi') { const pts = h.ring.map(i => B.coordOf(frame, i)); py = pts.reduce((a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]]).map(x => x / pts.length); v.addSphere({ center: { x: py[0], y: py[1], z: py[2] }, radius: 0.22, color: '#0e7c86' }); }
          else py = B.coordOf(frame, h.Y);
          v.addCylinder({ start: { x: px[0], y: px[1], z: px[2] }, end: { x: py[0], y: py[1], z: py[2] }, radius: 0.07, color: focus === h ? '#e2483d' : '#0e7c86', dashed: true, fromCap: 1, toCap: 1 });
          if (focus === h) v.addLabel(`${f2(h.d)} Å · ${f2(h.theta1, 1)}°`, { position: { x: (px[0] + py[0]) / 2, y: (px[1] + py[1]) / 2, z: (px[2] + py[2]) / 2 }, backgroundColor: '#101722', fontColor: 'white', fontSize: 12, backgroundOpacity: 0.85 });
        });
      }
      if (focus) { const p = B.coordOf(S.frames[frameIdx || 0], focus.X); v.zoomTo({ within: { distance: 7, sel: { resi: S.residues[focus.donorRes].resSeq, resn: S.residues[focus.donorRes].resName } } }); v.center({ x: p[0], y: p[1], z: p[2] }); }
      else if (hits && hits.length) v.zoomTo({ hetflag: true, resn: [...new Set(hits.map(h => S.residues[h.donorRes].resName))] });
      else v.zoomTo();
      v.render();
    },
    qmModel(model, pcs) {
      if (this.offline()) return; this.ensure(); const v = this.viewer; v.clear();
      v.addModel(B.writeXYZ(model), 'xyz'); v.setStyle({}, { stick: { radius: 0.13 }, sphere: { scale: 0.22 } });
      model.qm.forEach(a => { if (a.link) v.addSphere({ center: { x: a.x, y: a.y, z: a.z }, radius: 0.28, color: '#e87ba4', opacity: 0.85 }); });
      model.manifest.xbs.forEach(x => { const A = model.qm[x.X - 1]; let Y;
        if (x.Y) Y = model.qm[x.Y - 1]; else { const pts = x.ring.map(n => model.qm[n - 1]); Y = { x: pts.reduce((s, p) => s + p.x, 0) / pts.length, y: pts.reduce((s, p) => s + p.y, 0) / pts.length, z: pts.reduce((s, p) => s + p.z, 0) / pts.length }; }
        v.addCylinder({ start: { x: A.x, y: A.y, z: A.z }, end: { x: Y.x, y: Y.y, z: Y.z }, radius: 0.07, color: '#0e7c86', dashed: true }); });
      if (pcs && pcs.charges.length <= 4000) pcs.charges.forEach(c => { if (Math.abs(c.q) > 0.3) v.addSphere({ center: { x: c.x, y: c.y, z: c.z }, radius: 0.15, color: c.q > 0 ? '#2a78d6' : '#e2483d', opacity: 0.55 }); });
      v.zoomTo(); v.render();
    },
    material(M, hit) {
      if (this.offline()) return; this.ensure(); const v = this.viewer; v.clear();
      const atoms = M.atoms.map(a => ({ elem: a.el[0] + a.el.slice(1).toLowerCase(), x: a.c[0], y: a.c[1], z: a.c[2] }));
      const m = v.addModel(); m.addAtoms(atoms); v.setStyle({}, { stick: { radius: 0.12 }, sphere: { scale: 0.2 } });
      if (hit) { const X = M.atoms[hit.X].c; v.addSphere({ center: { x: hit.Ypos[0], y: hit.Ypos[1], z: hit.Ypos[2] }, radius: 0.3, color: '#e2483d', opacity: .7 }); v.addCylinder({ start: { x: X[0], y: X[1], z: X[2] }, end: { x: hit.Ypos[0], y: hit.Ypos[1], z: hit.Ypos[2] }, radius: 0.07, color: '#0e7c86', dashed: true }); }
      if (M.M) { try { m.setCrystData(M.cell.a, M.cell.b, M.cell.c, M.cell.alpha, M.cell.beta, M.cell.gamma); v.addUnitCell(m, { box: { color: '#9aa1ad' } }); } catch (e) { /* optional */ } }
      v.zoomTo(); v.render();
    }
  };

  /* ------------------------------ parameters: forms, assumptions ------------------------------ */
  const SPEC = B.PARAM_SPEC;
  const specOf = path => SPEC.find(s => s.path === path);
  const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
  const fmtVal = v => typeof v === 'boolean' ? (v ? 'on' : 'off') : (v && typeof v === 'object') ? '' : String(v);
  function inputFor(s) {
    const v = B.getPath(ST.P, s.path); const d = B.getPath(DEF, s.path);
    if (s.type === 'number') return `<input type="number" data-p="${s.path}" data-t="number" value="${esc(v)}" step="${s.step || 1}">`;
    if (s.type === 'text') return `<input type="text" data-p="${s.path}" data-t="text" value="${esc(v)}" style="width:150px">`;
    if (s.type === 'bool') return `<label class="chk"><input type="checkbox" data-p="${s.path}" data-t="bool" ${v ? 'checked' : ''}> ${v ? 'on' : 'off'}</label>`;
    if (s.type === 'select') { const opts = s.optionsFn ? s.optionsFn(ST.P) : s.options; return `<select data-p="${s.path}" data-t="${typeof d === 'number' ? 'number' : 'text'}">${opts.map(([k, l]) => `<option value="${esc(k)}" ${String(k) === String(v) ? 'selected' : ''}>${esc(l)}</option>`).join('')}</select>`; }
    if (s.type === 'cutoffs') {
      const xs = ['CL', 'BR', 'I'].concat(ST.P.geom.includeF ? ['F'] : []);
      return `<table class="mini"><tr><th>X \\ Y</th><th>N</th><th>O</th><th>S</th><th>π</th></tr>${xs.map(x => `<tr><td>${elTag(x)}</td>${['N', 'O', 'S', 'AR'].map(y => `<td><input type="number" step="0.01" data-p="geom.cutoffs.${x}.${y}" data-t="number" value="${v[x][y]}" title="default ${d[x][y]}"></td>`).join('')}</tr>`).join('')}</table>`;
    }
    if (s.type === 'vdw' || s.type === 'xh') return `<div class="row">${Object.keys(v).map(k => `<label class="f">${HN(k)}<input type="number" step="0.01" data-p="${s.path}.${k}" data-t="number" value="${v[k]}" style="width:66px" title="default ${d[k]}"></label>`).join('')}</div>`;
    if (s.type === 'weights') return `<div class="row">${['G', 'D', 'E', 'C'].map(k => `<label class="f">${k}<input type="number" step="0.05" min="0" data-p="${s.path}.${k}" data-t="number" value="${v[k]}" style="width:64px"></label>`).join('')}</div>`;
    return '';
  }
  function paramRow(s) {
    const changed = !same(B.getPath(ST.P, s.path), B.getPath(DEF, s.path));
    const d = B.getPath(DEF, s.path);
    const defTxt = (d && typeof d === 'object') ? Object.entries(d).map(([k, v]) => typeof v === 'object' ? '' : HN(k) + ' ' + v).filter(Boolean).join(', ') || 'hover a field to see it'
      : (s.type === 'select' ? (((s.optionsFn ? s.optionsFn(DEF) : s.options).find(o => String(o[0]) === String(d)) || [d, d])[1].split(' (')[0].split(' · ')[0]) : fmtVal(d));
    return `<div class="prow ${changed ? 'changed' : ''}" data-row="${s.path}"><div class="ptop"><div class="plab">${esc(s.label)}${s.unit ? ` <span class="muted">(${esc(s.unit)})</span>` : ''}${changed ? '<span class="chg">changed</span>' : ''}</div><div class="pin">${inputFor(s)}</div></div>
      <div class="pmeta">${s.desc ? `<span>${esc(s.desc)}</span>` : ''}<span class="pdef">default: ${esc(defTxt)}</span>${changed ? `<button class="btn xs" data-reset="${s.path}">reset</button>` : ''}</div></div>`;
  }
  function paramTable(paths) { return `<div class="ptable">${paths.map(p => specOf(p)).filter(Boolean).map(paramRow).join('')}</div>`; }
  function nChanged(paths) { return paths.filter(p => !same(B.getPath(ST.P, p), B.getPath(DEF, p))).length; }
  function paramsBlock(id, title, paths, open) {
    const n = nChanged(paths);
    return `<details class="params" data-pblock="${id}" ${open ? 'open' : ''}><summary>${esc(title)} · ${paths.length} setting${paths.length > 1 ? 's' : ''}${n ? `<span class="chg">${n} changed</span>` : ''}</summary>${paramTable(paths)}</details>`;
  }
  const PBLOCKS = {}; // id -> paths (to refresh)
  function pblock(id, title, paths, open) { PBLOCKS[id] = { title, paths }; return `<div data-pwrap="${id}">${paramsBlock(id, title, paths, open)}</div>`; }
  function refreshPBlocks() { $$('[data-pwrap]').forEach(w => { const id = w.dataset.pwrap; const b = PBLOCKS[id]; const open = w.querySelector('details') && w.querySelector('details').open; w.innerHTML = paramsBlock(id, b.title, b.paths, open); }); }
  function assumeBox(stage, title) { return `<details class="assume" open data-stage="${stage}"><summary>${esc(title || 'Assumptions used in this step')}</summary><ul>${B.assumptions(stage, ST.P).map(a => `<li>${esc(a)}</li>`).join('')}</ul></details>`; }
  function refreshAssumptions() { $$('.assume[data-stage]').forEach(d => { d.querySelector('ul').innerHTML = B.assumptions(d.dataset.stage, ST.P).map(a => `<li>${esc(a)}</li>`).join(''); }); }
  let refreshTimer = null;
  function onParamChanged(path) {
    log('param', `${path} = ${JSON.stringify(B.getPath(ST.P, path))}`);
    clearTimeout(refreshTimer); refreshTimer = setTimeout(() => refreshAfterParam(path), 0); // never rebuild DOM inside a change/blur handler
  }
  function refreshAfterParam(path) {
    refreshAssumptions(); if (current === 'params') R.params.render(); else refreshPBlocks();
    if (path.startsWith('geom.includeF')) refreshPBlocks();
    if (current === 'qm' && (path.startsWith('ff.') || path.startsWith('qm.env'))) R.qm.refreshEnv();
    if (current === 'quality') R.quality.render();
  }
  document.addEventListener('change', e => {
    const el = e.target; if (!el.dataset || !el.dataset.p) return;
    const t = el.dataset.t; let v = t === 'bool' ? el.checked : t === 'number' ? parseFloat(el.value) : el.value;
    if (t === 'number' && !Number.isFinite(v)) { toast('Please enter a number'); el.value = B.getPath(ST.P, el.dataset.p); return; }
    B.setPath(ST.P, el.dataset.p, v); onParamChanged(el.dataset.p);
  });
  document.addEventListener('click', e => {
    const b = e.target.closest && e.target.closest('[data-reset]'); if (!b) return;
    const p = b.dataset.reset; B.setPath(ST.P, p, JSON.parse(JSON.stringify(B.getPath(DEF, p)))); onParamChanged(p); toast('Reset ' + p);
  });
  const stagePaths = st => SPEC.filter(s => s.stage === st).map(s => s.path);

  /* ------------------------------ candidate registry ------------------------------ */
  function nextCid() { return 'XB-' + String(ST.cands.length + 1).padStart(3, '0'); }
  function addCandidate(c) {
    const dup = ST.cands.find(x => x.sourceKey === c.sourceKey); if (dup) return dup;
    c.cid = nextCid(); c.state = 'DETECTED'; c.created = new Date().toISOString(); ST.cands.push(c);
    log('registry.add', `${c.cid} ${c.label} (${c.source})`); updateSide(); return c;
  }
  function candFromHit(h, srcKind, structRef, frame, extra) {
    return Object.assign({ sourceKey: `${srcKind}|${ST.name}|${frame}|${h.key}`, source: srcKind === 'traj' ? `trajectory frame ${frame}` : (srcKind === 'cif' ? 'crystal' : 'static structure'),
      structName: ST.name, structRef, frame, hit: h, label: `${h.donorLabel}@${h.Xname} ··· ${h.accLabel}@${h.accName}`, halogen: h.halogen, accClass: h.accClass,
      d: h.d, theta1: h.theta1, theta2: h.theta2, G: h.G }, extra || {});
  }
  function stateBadge(c) { const m = { 'DETECTED': 'mute', 'MODELLED': 'acc', 'NBO-VERIFIED': 'ok', 'NBO-REJECTED': 'bad', 'EDA-DONE': 'ok' }; return `<span class="pill ${m[c.state] || 'mute'}">${c.state}</span>`; }
  function score(c) { return B.computeBXQI({ G: c.G, occupancy: c.occupancy, Dscore: c.Dscore, E2: c.E2, eda: c.eda }, ST.P); }
  function registryTable(opts = {}) {
    if (!ST.cands.length) return '<p class="muted">No candidates yet. Find halogen bonds (step 2, 3 or Materials) and press “Send to registry”.</p>';
    return `<div class="tablewrap"><table class="t"><thead><tr>${opts.select ? '<th><input type="checkbox" class="selAll"></th>' : ''}<th>ID</th><th>Interaction</th><th>X</th><th>Source</th><th>d (Å)</th><th>θ₁</th><th>G</th><th>Occ.</th><th>Model</th><th>E(2)</th><th>E<sub>int</sub></th><th>State</th></tr></thead><tbody>${ST.cands.map(c => `
      <tr data-cid="${c.cid}">${opts.select ? `<td><input type="checkbox" class="selc" value="${c.cid}" ${opts.isChecked && opts.isChecked(c) ? 'checked' : ''}></td>` : ''}<td class="mono">${c.cid}</td><td>${esc(c.label)}</td><td>${elTag(c.halogen)}</td><td class="muted">${esc(c.source)}</td>
      <td class="n">${f2(c.d)}</td><td class="n">${f2(c.theta1, 1)}</td><td class="n">${f2(c.G)}</td><td class="n">${c.occupancy != null ? f2(c.occupancy * 100, 0) + '%' : '—'}</td>
      <td class="mono">${c.modelId || '—'}</td><td class="n">${c.E2 != null ? f2(c.E2) : '—'}</td><td class="n">${c.eda ? f2(c.eda.total, 1) : '—'}</td><td>${stateBadge(c)}</td></tr>`).join('')}</tbody></table></div>`;
  }
  function bindSelectAll(root) { const a = root.querySelector('.selAll'); if (a) a.onchange = () => $$('.selc,.selh,.selt,.selm', root).forEach(c => { c.checked = a.checked; }); }

  /* ------------------------------ navigation ------------------------------ */
  const PAGES = [
    { id: 'structure', n: '1', label: 'Load structure', sub: 'PDB, PQR or RCSB ID', group: 'Workflow' },
    { id: 'detect', n: '2', label: 'Find halogen bonds', sub: 'geometry · G', group: 'Workflow' },
    { id: 'traj', n: '3', label: 'Dynamics', sub: 'optional · D', group: 'Workflow' },
    { id: 'qm', n: '4', label: 'QM model & charges', sub: 'cluster or QM/MM', group: 'Workflow' },
    { id: 'nbo', n: '5', label: 'NBO gate', sub: 'electronic check · E', group: 'Workflow' },
    { id: 'eda', n: '6', label: 'ALMO-EDA', sub: 'POL vs CT · C', group: 'Workflow' },
    { id: 'quality', n: '7', label: 'Score & report', sub: 'BXQI', group: 'Workflow' },
    { id: 'batch', n: 'B', label: 'Batch PDB', sub: 'many entries', group: 'Tools' },
    { id: 'material', n: 'M', label: 'Materials', sub: 'CIF / XYZ crystals', group: 'Tools' },
    { id: 'guide', n: 'T', label: 'Tutorial', sub: 'five pages · scope', group: 'Reference' },
    { id: 'params', n: 'P', label: 'All parameters', sub: 'defaults & changes', group: 'Reference' },
    { id: 'method', n: 'i', label: 'Method', sub: 'equations & data', group: 'Reference' }
  ];
  const DONE = new Set();
  function buildSide() {
    let g = null; let html = '';
    PAGES.forEach(p => { if (p.group !== g) { g = p.group; html += `<h6>${g}</h6>`; } html += `<button data-p="${p.id}"><span class="num">${p.n}</span><span class="lbl">${p.label}<small>${p.sub}</small></span></button>`; });
    $('#side').innerHTML = html; $$('#side button').forEach(b => b.onclick = () => go(b.dataset.p));
  }
  function setDone(id, on = true) { if (on) DONE.add(id); else DONE.delete(id); updateSide(); }
  function updateSide() { $$('#side button').forEach(b => { b.classList.toggle('done', DONE.has(b.dataset.p)); b.classList.toggle('on', b.dataset.p === current); }); }
  let current = null;
  function go(id) {
    current = id; $$('.page').forEach(p => p.classList.toggle('on', p.id === 'pg-' + id)); updateSide();
    const pg = $('#pg-' + id); const slot = pg.querySelector('.vslot'); if (slot) V.mount(slot);
    refreshAssumptions(); refreshPBlocks();
    if (R[id] && R[id].onShow) R[id].onShow();
    window.scrollTo({ top: 0 });
  }
  const nextBtn = (id, label) => `<div class="nextbar"><button class="btn primary" data-go="${id}">${esc(label)} →</button></div>`;
  document.addEventListener('click', e => { const b = e.target.closest && e.target.closest('[data-go]'); if (b) go(b.dataset.go); });
  const head = (kicker, title, text) => `<div class="pagehead"><div class="kicker">${esc(kicker)}</div><h2>${esc(title)}</h2><p>${text}</p></div>`;
  const R = {};


  /* ====================== TUTORIAL (five pages) ====================== */
  const TUT = [
    { t: '1 · What BondX Studio does', b: `
      <p>BondX Studio checks halogen bonds (XBs) in proteins, ligands and crystals. It runs in your browser. No file leaves your computer. The only network call is the download of a PDB entry from RCSB.</p>
      <p>The analysis is split into tiers. Each tier adds one kind of evidence.</p>
      <table class="t tut"><thead><tr><th>Tier</th><th>Step</th><th>Evidence</th></tr></thead><tbody>
        <tr><td class="mono">T0</td><td>1 Load structure</td><td>Atoms, chains, hydrogens, charges, virtual sites</td></tr>
        <tr><td class="mono">T1</td><td>2 Find halogen bonds</td><td>Geometry: distance, two angles, geometric score <i>G</i></td></tr>
        <tr><td class="mono">T2</td><td>3 Dynamics</td><td>Occupancy of each contact over a trajectory</td></tr>
        <tr><td class="mono">T3</td><td>4–5 QM model, NBO gate</td><td>Cluster or QM/MM input files; NBO donor–acceptor energy Σ<i>E</i><sup>(2)</sup></td></tr>
        <tr><td class="mono">T4</td><td>6 ALMO-EDA</td><td>Frozen, polarisation and charge-transfer energies</td></tr>
        <tr><td class="mono">T5</td><td>7 Score &amp; report</td><td>Quality index BXQI and grade A–D</td></tr>
      </tbody></table>
      <div class="note ok" style="margin-top:12px"><b>Quick start.</b> Open <b>1 Load structure</b>. Press <b>Load demo</b>. Then press the button at the bottom right of each page. The demo takes less than one minute.</div>
      <p class="muted">You can stop after any tier. A result from T1 alone is valid. It is then marked as geometry only.</p>` },
    { t: '2 · Load a structure and find contacts', b: `
      <h4>Step 1 · Load structure (T0)</h4>
      <ul>
        <li>Type a four-letter PDB ID and press <b>Fetch</b>. Or drop a local <span class="mono">.pdb</span>, <span class="mono">.ent</span> or <span class="mono">.pqr</span> file.</li>
        <li>A PQR file brings its own partial charges. A PDB file uses the force-field charges of the program.</li>
        <li>Pick the chains to keep. Check the composition table.</li>
        <li>Extra-point atoms (EP, LP, XX) of σ-hole force fields are removed on load. The message says how many were removed.</li>
      </ul>
      <h4>Step 2 · Find halogen bonds (T1)</h4>
      <ul>
        <li>The program finds every C–X···Y contact with X = Cl, Br or I.</li>
        <li>It measures the X···Y distance <i>d</i>, the C–X···Y angle θ<sub>1</sub> and the X···Y–R angle θ<sub>2</sub>.</li>
        <li>It reports the score <i>G</i> between 0 and 1. A larger <i>G</i> means a better geometry.</li>
        <li>Each acceptor is classed as a Lewis base or not. A backbone N–H has no lone pair. Such a contact is flagged as weak.</li>
        <li>Sort or filter the table. Click a row to see the contact in the 3D viewer.</li>
      </ul>
      <div class="note">All cut-offs are on the page <b>All parameters</b>. A changed value is marked in orange. It is also written to the report.</div>` },
    { t: '3 · Dynamics and the QM model', b: `
      <h4>Step 3 · Dynamics (T2, optional)</h4>
      <ul>
        <li>Drop a topology <span class="mono">.pdb</span>. Its atom order must match the trajectory.</li>
        <li>Drop a <span class="mono">.dcd</span>, a multi-model <span class="mono">.pdb</span> or a multi-frame <span class="mono">.xyz</span> file.</li>
        <li>Convert GROMACS or AMBER files once with <span class="mono">gmx trjconv</span> or <span class="mono">cpptraj</span>.</li>
        <li>The program applies the T1 test to every frame. It reports the occupancy of each contact.</li>
      </ul>
      <h4>Step 4 · QM model and charges (T3)</h4>
      <ul>
        <li>Choose a <b>cluster</b> model or a <b>QM/MM</b> model with point charges.</li>
        <li>Cut bonds are capped with hydrogen atoms. Ligand charges can come from a <span class="mono">.mol2</span> file.</li>
        <li>Pick the program: Gaussian, ORCA or Q-Chem. Pick the method, the basis set and the solvent model.</li>
        <li>Download the input files. Each file carries a manifest ID. The ID links the output back to the contact.</li>
        <li>Run the jobs on your own computer or cluster. BondX Studio does not run QM jobs.</li>
      </ul>` },
    { t: '4 · NBO gate, ALMO-EDA and the score', b: `
      <h4>Step 5 · NBO gate (T3)</h4>
      <ul>
        <li>Drop the QM output (<span class="mono">.log</span> or <span class="mono">.out</span>). The program reads the second-order perturbation table.</li>
        <li>It adds the terms LP(Y) → σ*(C–X) and LP(Y) → RY*(X). For a π acceptor it adds the ring bonds.</li>
        <li>Intramolecular terms are excluded. Only the T1 acceptor counts.</li>
        <li>A contact passes the gate if Σ<i>E</i><sup>(2)</sup> ≥ 0.5 kcal/mol.</li>
        <li>For an ionic NBO structure the halogen has a vacant LP*. Donation into LP*(X) is then counted.</li>
      </ul>
      <h4>Step 6 · ALMO-EDA (T4)</h4>
      <ul>
        <li>Drop the Q-Chem EDA output. The program reads FRZ, POL, CT and the total energy.</li>
        <li>It reports χ<sub>CT</sub> = CT / (POL + CT).</li>
      </ul>
      <h4>Step 7 · Score and report (T5)</h4>
      <ul>
        <li>BXQI is a weighted geometric mean of the tier scores. It runs from 0 to 100.</li>
        <li>Grades: A ≥ 75, B ≥ 50, C ≥ 25, else D. Grade A needs a passed NBO gate.</li>
        <li>Press <b>Report</b> for a single HTML file. It lists every input, parameter and result.</li>
      </ul>` },
    { t: '5 · What it can and cannot do', b: `
      <div class="grid g2e">
        <div class="card ok-card"><h2>It can</h2><ul>
          <li>Read PDB, PQR, CIF, XYZ, DCD and multi-model PDB files.</li>
          <li>Find Cl, Br and I halogen bonds to O, N, S and aromatic rings.</li>
          <li>Screen many PDB entries in one batch.</li>
          <li>Build cluster and QM/MM models with ff14SB, ff19SB or CHARMM36 charges.</li>
          <li>Write Gaussian, ORCA and Q-Chem inputs.</li>
          <li>Read Gaussian NBO output and Q-Chem ALMO-EDA output.</li>
          <li>Save and reload a full session as JSON.</li>
          <li>Write a complete HTML report and CSV tables.</li>
        </ul></div>
        <div class="card no-card"><h2>It cannot</h2><ul>
          <li>Run quantum-chemistry jobs. You run them yourself.</li>
          <li>Read GROMACS <span class="mono">.xtc</span> or AMBER <span class="mono">.nc</span> files directly.</li>
          <li>Give a tested result for fluorine. Fluorine is off by default.</li>
          <li>Replace a full free-energy study of binding.</li>
          <li>Decide on a weak contact without the NBO tier.</li>
        </ul></div>
      </div>
      <h4>Tested routes</h4>
      <table class="t tut"><thead><tr><th>Task</th><th>Program</th><th>Status</th></tr></thead><tbody>
        <tr><td>NBO gate</td><td>Gaussian 16 with NBO 3.1</td><td><span class="pill ok">tested on 59 dimers and 9 protein contacts</span></td></tr>
        <tr><td>ALMO-EDA</td><td>Q-Chem 6.2</td><td><span class="pill ok">tested</span></td></tr>
        <tr><td>Batch screening</td><td>local PDB files or PDB IDs</td><td><span class="pill ok">tested on three test-set proteins</span></td></tr>
        <tr><td>Crystal contacts and QM model</td><td>CIF with symmetry operators</td><td><span class="pill ok">tested on a small test CIF</span></td></tr>
        <tr><td>Counterpoise binding</td><td>ORCA 5 or Gaussian 16</td><td><span class="pill mute">used in the paper; not read by the app</span></td></tr>
        <tr><td>NBO gate</td><td>ORCA with NBO</td><td><span class="pill warn">input written; not tested</span></td></tr>
        <tr><td>NBO gate</td><td>Q-Chem with NBO 5</td><td><span class="pill warn">input written; not tested</span></td></tr>
      </tbody></table>
      <p class="muted" style="margin-top:10px">Source code, tests and paper scripts: <a href="https://github.com/compchemsk/BondX-Studio" target="_blank" rel="noopener">github.com/compchemsk/BondX-Studio</a>.</p>` }
  ];
  R.guide = {
    html: () => `${head('Reference', 'Tutorial', 'Five short pages. Read them once before your first run.')}
      <div class="tutwrap"><div class="tutnav" id="tutNav">${TUT.map((p, i) => `<button data-t="${i}">${esc(p.t)}</button>`).join('')}</div>
      <div class="card tutcard"><div id="tutBody"></div>
      <div class="tutfoot"><button class="btn" id="tutPrev">← Previous</button><span class="tutdots" id="tutDots">${TUT.map((_, i) => `<i data-t="${i}"></i>`).join('')}</span><button class="btn primary" id="tutNext">Next →</button></div></div></div>`,
    init() {
      let k = 0;
      const show = i => { k = Math.max(0, Math.min(TUT.length - 1, i)); $('#tutBody').innerHTML = `<h2 class="tuth">${esc(TUT[k].t)}</h2><div class="tutpage">${TUT[k].b}</div><p class="muted tutpg">Page ${k + 1} of ${TUT.length}</p>`;
        $$('#tutNav button').forEach((b, j) => b.classList.toggle('on', j === k)); $$('#tutDots i').forEach((d, j) => d.classList.toggle('on', j === k));
        $('#tutPrev').disabled = k === 0; $('#tutNext').textContent = k === TUT.length - 1 ? 'Start: load a structure →' : 'Next →'; };
      $$('#tutNav button, #tutDots i').forEach(b => b.onclick = () => show(+b.dataset.t));
      $('#tutPrev').onclick = () => show(k - 1); $('#tutNext').onclick = () => { if (k === TUT.length - 1) go('structure'); else show(k + 1); };
      show(0);
    }
  };

  /* ====================== 1 STRUCTURE ====================== */
  R.structure = {
    html: () => `${head('Step 1 · T0', 'Load a structure', 'Fetch a PDB entry from RCSB or open a local PDB/PQR file (a PQR file supplies its own partial charges). Choose chains and check the composition before you continue.')}
    ${assumeBox('T0')}
    <div class="grid g2">
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Load</h2>
          <div class="row"><input type="text" id="pdbId" placeholder="PDB ID, e.g. 1US0" style="width:150px;text-transform:uppercase"><button class="btn primary" id="btnFetch">Fetch from RCSB</button></div>
          <div class="row"><div class="drop" id="dropPdb" style="flex:1">Drop or click to open <b>.pdb / .ent / .pqr</b><input type="file" accept=".pdb,.ent,.pqr,.txt"></div></div>
          <div class="row"><button class="btn sm" id="btnDemo">Load built-in demo complex</button><span class="muted" style="font-size:12px">iodobenzene → backbone C=O, bromobenzene → Phe π</span></div>
          <div id="loadMsg" style="margin-top:8px"></div>
          <div style="margin-top:10px">${pblock('st', 'Reading settings', stagePaths('T0'))}<div class="row" style="margin-top:6px"><button class="btn sm" id="btnReparse" disabled>Re-read file with these settings</button></div></div>
        </div>
        <div class="card" id="chainCard" style="display:none"><h2>Chains</h2><p class="sub">Keep only the chains you need.</p>
          <div id="chainList" class="row"></div>
          <div class="row" style="margin-top:10px"><button class="btn" id="btnApplyChains">Apply selection</button><button class="btn" id="btnDlProcessed">Download processed PDB</button></div>
        </div>
      </div>
      <div class="grid" style="align-content:start">
        <div class="vslot"></div>
        <div class="card" id="reportCard" style="display:none"><h2 id="repTitle">Composition</h2><div id="repKpis" class="kpis"></div><div class="tablewrap" style="max-height:300px"><table class="t" id="repTable"></table></div></div>
      </div>
    </div>${nextBtn('detect', 'Next: find halogen bonds')}`,
    init() {
      $('#btnFetch').onclick = async () => {
        const id = $('#pdbId').value.trim().toUpperCase(); if (!/^[0-9][A-Z0-9]{3}$/.test(id)) return toast('Enter a valid 4-character PDB ID');
        $('#loadMsg').innerHTML = '<span class="muted">Downloading ' + id + ' …</span>';
        try { const r = await fetch(`https://files.rcsb.org/download/${id}.pdb`); if (!r.ok) throw new Error('HTTP ' + r.status + ' (large entries are distributed only as mmCIF)'); await loadStructure(await r.text(), id, 'RCSB ' + id); }
        catch (e) { $('#loadMsg').innerHTML = `<div class="note bad">Could not fetch ${id}: ${esc(e.message)}. Download the PDB-format file manually and open it.</div>`; }
      };
      $('#pdbId').addEventListener('keydown', e => { if (e.key === 'Enter') $('#btnFetch').click(); });
      dropZone($('#dropPdb'), async files => { const f = files[0]; await loadStructure(await readText(f), f.name.replace(/\.(pdb|ent|pqr|txt)$/i, ''), 'file ' + f.name, /\.pqr$/i.test(f.name)); });
      $('#btnDemo').onclick = () => loadStructure($('#demoPdb').textContent.trim() + '\n', 'BondX_DEMO', 'built-in demo');
      $('#btnReparse').onclick = () => { if (ST.srcText) loadStructure(ST.srcText, ST.name, 're-read', ST.srcPQR); };
      $('#btnApplyChains').onclick = () => { const keep = $$('#chainList input:checked').map(i => i.value); ST.S = B.filterChains(ST.S0, keep); ST.hits = []; log('chains.apply', keep.join(',')); renderReport(); V.structure(ST.S); toast('Kept chains ' + keep.join(', ')); };
      $('#btnDlProcessed').onclick = () => download(`${ST.name}_processed.pdb`, B.writePDB(ST.S, null, null, `BondX processed structure (${ST.runId})`));
    }
  };
  async function loadStructure(text, name, origin, pqr) {
    try {
      const S = B.parsePDB(text, { pqr, altloc: ST.P.struct.altloc, bondTol: ST.P.struct.bondTol }); if (!S.atoms.length) throw new Error('no ATOM/HETATM records found');
      ST.S0 = S; ST.S = S; ST.name = name; ST.srcText = text; ST.srcPQR = !!pqr; ST.hits = [];
      const h = await sha256(text); log('structure.load', `${name} from ${origin}; ${S.atoms.length} atoms; ${S.frames.length} model(s); sha256=${h}; altloc=${ST.P.struct.altloc}, bondTol=${ST.P.struct.bondTol}`); ST.srcHash = h;
      $('#loadMsg').innerHTML = `<div class="note ok">Loaded <b>${esc(name)}</b> — ${S.atoms.length.toLocaleString()} atoms, ${S.residues.length.toLocaleString()} residues${S.frames.length > 1 ? `, ${S.frames.length} models (also usable as a trajectory)` : ''}${S.isPQR ? ', <b>PQR charges present</b>' : ''}${S.vsites ? `, <b>${S.vsites.count} virtual site(s) removed</b> (${esc(S.vsites.names.join(', '))}; extra-point or lone-pair charges of the force field, not atoms)` : ''}${S.hasH ? ', hydrogens present' : ', no hydrogens (they will be added where needed)'}.<br><span class="mono" style="font-size:11px">SHA-256 ${h.slice(0, 32)}…</span></div>`;
      $('#chainList').innerHTML = [...new Set(S.atoms.map(a => a.chain))].map(c => `<label class="chk" style="border:1px solid var(--line);border-radius:7px;padding:4px 8px"><input type="checkbox" value="${esc(c)}" checked> Chain ${esc(c === '_' ? '(blank)' : c)}</label>`).join('');
      $('#chainCard').style.display = ''; $('#btnReparse').disabled = false; renderReport(); setDone('structure'); V.structure(S);
    } catch (e) { $('#loadMsg').innerHTML = `<div class="note bad">Could not read the file: ${esc(e.message)}</div>`; }
  }
  function renderReport() {
    const S = ST.S; const rep = B.structureReport(S);
    $('#reportCard').style.display = ''; $('#repTitle').textContent = `Composition — ${ST.name}`;
    const hal = S.residues.filter(r => r.cls === 'ligand' && r.halogens.some(x => x !== 'F'));
    $('#repKpis').innerHTML = [[S.atoms.length.toLocaleString(), 'atoms'], [S.residues.filter(r => r.cls === 'protein').length, 'amino acids'], [S.residues.filter(r => r.cls === 'ligand').length, 'ligands / het groups'], [hal.length, 'halogenated ligands'], [S.residues.filter(r => r.cls === 'water').length, 'waters'], [S.header.resolution ? S.header.resolution + ' Å' : '—', 'resolution']]
      .map(k => `<div class="kpi"><b>${k[0]}</b><span>${k[1]}</span></div>`).join('');
    $('#repTable').innerHTML = `<thead><tr><th>Chain</th><th>Molecule</th><th>Length</th><th>Range</th><th>Breaks</th><th>Waters</th><th>Ligands</th><th>Ions</th><th>Halogenated</th></tr></thead><tbody>${rep.map(r => `<tr><td class="mono">${esc(r.chain)}</td><td>${esc(r.molecule)}</td><td class="n">${r.length || ''}</td><td class="mono">${esc(r.range)}</td><td class="n">${r.breaks || 0}</td><td class="n">${r.waters}</td><td>${esc(r.ligands.join(', ') || '—')}</td><td>${esc(r.ions.join(', ') || '—')}</td><td>${r.halogenLigands.length ? r.halogenLigands.map(x => `<span class="pill warn">${esc(x)}</span>`).join(' ') : '<span class="muted">—</span>'}</td></tr>`).join('')}</tbody>`;
  }

  /* ====================== 2 DETECT ====================== */
  function hitsTable(hits, opts = {}) {
    return `<table class="t"><thead><tr>${opts.select ? '<th><input type="checkbox" class="selAll"></th>' : ''}<th>#</th><th>Donor (A–X)</th><th>Acceptor</th><th>Class</th><th>d (Å)</th><th>Σr<sub>vdW</sub></th><th>R</th><th>θ₁ (°)</th><th>θ₂ / α (°)</th><th>Lewis</th><th>G</th><th>Status</th></tr></thead><tbody>${hits.map((h, k) => `
      <tr data-k="${k}">${opts.select ? `<td><input type="checkbox" class="selh" value="${k}" ${h.status === 'XB' ? 'checked' : ''}></td>` : ''}<td class="mono">${k + 1}</td><td>${elTag(h.halogen)} ${esc(h.donorLabel)} <span class="muted">${esc(h.Aname)}–${esc(h.Xname)}</span></td>
      <td>${esc(h.accLabel)} <b>${esc(h.accName)}</b></td><td class="muted">${esc(h.accClass)}</td><td class="n">${f2(h.d, 3)}</td><td class="n">${f2(h.sumVdw)}</td><td class="n">${f2(h.R, 3)}</td><td class="n">${f2(h.theta1, 1)}</td><td class="n">${f2(h.theta2, 1)}</td>
      <td>${h.lewis === 'weak' ? '<span class="pill warn" data-tip="Delocalised or protonated lone pair: poor XB acceptor. The NBO gate will decide.">weak</span>' : '<span class="pill mute">strong</span>'}</td>
      <td><div class="bar" data-tip="G = ${f2(h.G, 3)} (g_d ${f2(h.gd, 2)}, g_θ ${f2(h.gt, 2)}, g₂ ${f2(h.g2, 2)})"><i style="width:${(h.G * 100).toFixed(0)}%"></i></div></td><td>${h.status === 'XB' ? '<span class="pill ok">XB</span>' : '<span class="pill mute">near-miss</span>'}</td></tr>`).join('')}</tbody></table>`;
  }
  const GEO_MAIN = ['geom.cutoffs', 'geom.theta1Min', 'geom.theta1Max', 'geom.includePi', 'geom.acceptorScope.water', 'geom.reportNearMiss'];
  const GEO_MORE = stagePaths('T1').filter(p => !GEO_MAIN.includes(p) && !p.startsWith('mat.'));
  R.detect = {
    html: () => `${head('Step 2 · T1', 'Find halogen bonds by geometry', 'Every C–X (X = Cl, Br, I) of ligands and modified residues is tested against protein, nucleic-acid (and optionally water or ligand) acceptors. Each contact gets a geometric score G between 0 and 1.')}
    ${assumeBox('T1')}
    <div class="grid g2">
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Criteria</h2>
          ${pblock('geoMain', 'Main criteria', GEO_MAIN, true)}
          <div style="margin-top:8px">${pblock('geoMore', 'More settings (acceptors, radii, scoring)', GEO_MORE)}</div>
          <div class="row" style="margin-top:14px"><button class="btn primary" id="btnDetect">Run detection</button></div></div>
        <div class="card"><h2>Distance–angle map</h2><p class="sub">Shape and colour encode the halogen.</p><div id="dScatter"></div></div>
      </div>
      <div class="grid" style="align-content:start">
        <div class="vslot"></div>
        <div class="card"><div class="row" style="justify-content:space-between"><h2 id="dTitle">Results</h2><div class="row"><button class="btn sm" id="btnToReg" disabled>Send selected to registry →</button><button class="btn sm" id="btnCsv" disabled>CSV</button><button class="btn sm" id="btnHbData" disabled data-tip="Legacy notebook format">hb_data.txt</button></div></div>
          <div id="dKpis" class="kpis"></div><div class="tablewrap" id="dTable"><p class="muted" style="padding:10px">Load a structure and run detection.</p></div></div>
      </div>
    </div>${nextBtn('qm', 'Next: build the QM model')}`,
    init() {
      $('#btnDetect').onclick = () => this.run();
      $('#btnCsv').onclick = () => download(`${ST.name}_XB_T1.csv`, B.toCSV(ST.hits, [{ h: 'donor', k: 'donorLabel' }, { h: 'A', k: 'Aname' }, { h: 'X', k: 'Xname' }, { h: 'halogen', k: 'halogen' }, { h: 'acceptor', k: 'accLabel' }, { h: 'acc_atom', k: 'accName' }, { h: 'class', k: 'accClass' }, { h: 'd_A', k: 'd' }, { h: 'sum_vdW', k: 'sumVdw' }, { h: 'R', k: 'R' }, { h: 'theta1', k: 'theta1' }, { h: 'theta2_or_alpha', k: 'theta2' }, { h: 'lewis', k: 'lewis' }, { h: 'g_d', k: 'gd' }, { h: 'g_theta', k: 'gt' }, { h: 'g2', k: 'g2' }, { h: 'G', k: 'G' }, { h: 'status', k: 'status' }]));
      $('#btnHbData').onclick = () => download('hb_data.txt', 'Ligand_Name Halogen_Atom Residue_Name Residue_ID Heavy_Atom Distance Angle\n' + ST.hits.filter(h => h.status === 'XB').map(h => [ST.S.residues[h.donorRes].resName, h.Xname, ST.S.residues[h.accRes].resName, ST.S.residues[h.accRes].resSeq, h.type === 'pi' ? 'RING_CENTROID' : h.accName, h.d.toFixed(3), h.theta1.toFixed(2)].join(' ')).join('\n') + '\n');
      $('#btnToReg').onclick = () => { const sel = $$('#dTable .selh:checked').map(i => ST.hits[+i.value]); sel.forEach(h => addCandidate(candFromHit(h, 'pdb', 'S', 0))); toast(`${sel.length} candidate(s) in registry`); setDone('detect'); };
    },
    onShow() { if (ST.S && !ST.hits.length) V.structure(ST.S); },
    run() {
      if (!ST.S) return toast('Load a structure first');
      const P = ST.P; const t0 = performance.now(); const { hits } = B.detectXB(ST.S, P); ST.hits = hits;
      log('T1.detect', `${hits.length} contacts; geom=${JSON.stringify(P.geom)}; score=${JSON.stringify(P.score)}`);
      const xb = hits.filter(h => h.status === 'XB');
      $('#dTitle').textContent = `Results — ${ST.name}`;
      const byX = {}; xb.forEach(h => byX[h.halogen] = (byX[h.halogen] || 0) + 1);
      $('#dKpis').innerHTML = [[xb.length, 'halogen bonds'], [xb.filter(h => h.type === 'pi').length, 'X···π'], [xb.filter(h => h.accClass === 'backbone C=O').length, 'to backbone C=O'], [hits.length - xb.length, 'near-misses'], [(performance.now() - t0).toFixed(0) + ' ms', 'runtime']].concat(Object.entries(byX).map(([k, v]) => [v, HN(k) + '-bonds']))
        .map(k => `<div class="kpi"><b>${k[0]}</b><span>${k[1]}</span></div>`).join('');
      if (!hits.length) { $('#dTable').innerHTML = `<div class="note" style="margin:10px">${ST.S.residues.some(r => r.cls === 'ligand' && r.halogens.some(x => x !== 'F')) ? 'A halogenated ligand is present, but no contact meets the criteria.' : 'No halogenated ligand (Cl/Br/I) in the selected chains.'}</div>`; ['#btnToReg', '#btnCsv', '#btnHbData'].forEach(b => $(b).disabled = true); V.structure(ST.S); Charts.scatter($('#dScatter'), { points: [], groups: [] }); return; }
      $('#dTable').innerHTML = hitsTable(hits, { select: true }); bindSelectAll($('#dTable'));
      ['#btnToReg', '#btnCsv', '#btnHbData'].forEach(b => $(b).disabled = false);
      $$('#dTable tbody tr').forEach(tr => tr.addEventListener('click', e => { if (e.target.tagName === 'INPUT') return; $$('#dTable tr').forEach(t => t.classList.remove('sel')); tr.classList.add('sel'); V.structure(ST.S, 0, ST.hits.filter(h => h.status === 'XB'), ST.hits[+tr.dataset.k]); }));
      V.structure(ST.S, 0, xb);
      const groups = ['I', 'BR', 'CL', 'F'].filter(x => hits.some(h => h.halogen === x)).map(x => ({ key: x, name: HN(x), color: HCOLOR[x], shape: HSHAPE[x] }));
      Charts.scatter($('#dScatter'), { points: hits.map(h => ({ x: h.R, y: h.theta1, group: h.halogen, tip: `<b>${esc(h.donorLabel)}@${esc(h.Xname)}</b> → ${esc(h.accLabel)} ${esc(h.accName)}<br>d ${f2(h.d, 3)} Å · R ${f2(h.R, 3)} · θ₁ ${f2(h.theta1, 1)}°<br>G ${f2(h.G, 3)} · ${h.status}` })),
        groups, xDomain: [Math.min(0.75, ...hits.map(h => h.R)), Math.max(1.15, ...hits.map(h => h.R))], yDomain: [Math.min(P.geom.theta1Min - 10, ...hits.map(h => h.theta1)), 180], xLabel: 'R = d / Σr_vdW', yLabel: 'θ₁ A–X···Y (°)', refLines: [{ y: P.geom.theta1Min, label: `θ₁ ≥ ${P.geom.theta1Min}°` }, { x: 1, label: 'Σ r_vdW' }], fmtX: v => v.toFixed(2) });
    }
  };

  /* ====================== BATCH ====================== */
  R.batch = {
    html: () => `${head('Tool · T1', 'Batch screening of many PDB entries', 'Screens a list of PDB IDs (downloaded one after another from RCSB) or several local files with the criteria of step 2. For thousands of entries use the command-line batch tool.')}
    ${assumeBox('T1')}
    <div class="grid g2">
      <div class="card" style="align-self:start"><h2>Entries</h2><p class="sub">Comma, space or newline separated.</p><textarea id="bIds" rows="5" placeholder="1US0, 3OS3, 2YNN"></textarea>
        <div class="row" style="margin-top:8px"><div class="drop" id="dropBatch" style="flex:1">…or drop several local <b>.pdb</b> files<input type="file" accept=".pdb,.ent"></div></div>
        <div class="row" style="margin-top:10px"><button class="btn primary" id="btnBatch">Run batch</button><label class="chk"><input type="checkbox" id="bNear"> include near-misses</label></div>
        <div class="progress"><i id="bProg"></i></div><div id="bLog" class="muted mono" style="font-size:11.5px;max-height:160px;overflow:auto"></div></div>
      <div class="grid" style="align-content:start">
        <div class="card"><div class="row" style="justify-content:space-between"><h2>Aggregated results</h2><button class="btn sm" id="btnBCsv" disabled>CSV</button></div><div id="bKpis" class="kpis"></div><div class="tablewrap" id="bTable"><p class="muted" style="padding:10px">No batch run yet.</p></div></div>
        <div class="card"><h2>Halogen bonds per entry</h2><div id="bChart"></div></div>
      </div>
    </div>`,
    init() {
      let localFiles = [];
      dropZone($('#dropBatch'), files => { localFiles = files; $('#dropBatch').firstChild.textContent = `${files.length} local file(s) queued `; }, true);
      $('#btnBatch').onclick = async () => {
        const ids = $('#bIds').value.toUpperCase().split(/[\s,;]+/).filter(x => /^[0-9][A-Z0-9]{3}$/.test(x));
        const jobs = ids.map(id => ({ id, get: async () => { const r = await fetch(`https://files.rcsb.org/download/${id}.pdb`); if (!r.ok) throw new Error('HTTP ' + r.status); return r.text(); } }))
          .concat(localFiles.map(f => ({ id: f.name.replace(/\.\w+$/, ''), get: () => readText(f) })));
        if (!jobs.length) return toast('Add PDB IDs or files');
        const P = B.mergeParams(ST.P, { geom: { reportNearMiss: $('#bNear').checked } }); ST.batch = []; const logEl = $('#bLog'); logEl.innerHTML = '';
        for (let k = 0; k < jobs.length; k++) {
          const j = jobs[k]; $('#bProg').style.width = (k / jobs.length * 100) + '%';
          try { const S = B.parsePDB(await j.get(), { altloc: P.struct.altloc, bondTol: P.struct.bondTol }); const { hits } = B.detectXB(S, P);
            hits.forEach(h => ST.batch.push({ pdb: j.id, ligand: S.residues[h.donorRes].resName, halogen: h.halogen, X: h.Xname, residue: S.residues[h.accRes].resName, resid: S.residues[h.accRes].resSeq, chain: S.residues[h.accRes].chain, atom: h.accName, cls: h.accClass, d: h.d, R: h.R, theta1: h.theta1, G: h.G, status: h.status }));
            const nHal = S.residues.filter(r => r.cls === 'ligand' && r.halogens.some(x => x !== 'F')).length;
            logEl.innerHTML += `<div>[${k + 1}/${jobs.length}] ${esc(j.id)}: ${S.atoms.length} atoms, ${nHal} halogenated ligand(s), ${hits.filter(h => h.status === 'XB').length} XB</div>`;
          } catch (e) { logEl.innerHTML += `<div style="color:var(--bad)">[${k + 1}/${jobs.length}] ${esc(j.id)}: ${esc(e.message)}</div>`; }
          logEl.scrollTop = 1e9; await new Promise(r => setTimeout(r, 0));
        }
        $('#bProg').style.width = '100%'; log('T1.batch', `${jobs.length} entries; ${ST.batch.length} contacts`); this.render(jobs.map(j => j.id)); setDone('batch');
      };
      $('#btnBCsv').onclick = () => download('BondX_batch_XB.csv', B.toCSV(ST.batch, ['pdb', 'ligand', 'halogen', 'X', 'residue', 'resid', 'chain', 'atom', 'cls', 'd', 'R', 'theta1', 'G', 'status'].map(k => ({ h: k, k }))));
    },
    render(ids) {
      const rows = ST.batch; $('#btnBCsv').disabled = !rows.length; const xb = rows.filter(r => r.status === 'XB');
      $('#bKpis').innerHTML = [[ids.length, 'entries'], [new Set(xb.map(r => r.pdb)).size, 'entries with XB'], [xb.length, 'halogen bonds'], [xb.length ? f2(xb.reduce((s, r) => s + r.d, 0) / xb.length) + ' Å' : '—', 'mean d']].map(k => `<div class="kpi"><b>${k[0]}</b><span>${k[1]}</span></div>`).join('');
      $('#bTable').innerHTML = rows.length ? `<table class="t"><thead><tr><th>PDB</th><th>Ligand</th><th>X</th><th>Acceptor</th><th>Class</th><th>d (Å)</th><th>R</th><th>θ₁</th><th>G</th><th>Status</th></tr></thead><tbody>${rows.map(r => `<tr><td class="mono">${esc(r.pdb)}</td><td>${esc(r.ligand)} <span class="muted">${esc(r.X)}</span></td><td>${elTag(r.halogen)}</td><td>${esc(r.residue)}${r.resid}:${esc(r.chain)} <b>${esc(r.atom)}</b></td><td class="muted">${esc(r.cls)}</td><td class="n">${f2(r.d, 3)}</td><td class="n">${f2(r.R, 3)}</td><td class="n">${f2(r.theta1, 1)}</td><td class="n">${f2(r.G)}</td><td>${r.status === 'XB' ? '<span class="pill ok">XB</span>' : '<span class="pill mute">near-miss</span>'}</td></tr>`).join('')}</tbody></table>` : '<p class="muted" style="padding:10px">No halogen contacts found in these entries.</p>';
      Charts.hbars($('#bChart'), { rows: ids.map(id => { const v = xb.filter(r => r.pdb === id).length; return { label: id, value: v, text: String(v), color: 'var(--s1)', tip: `${id}: ${v} halogen bond(s)` }; }), xLabel: 'halogen bonds (count)', labelWidth: 90, fmt: v => v });
    }
  };

  /* ====================== 3 TRAJECTORY ====================== */
  R.traj = {
    html: () => `${head('Step 3 · T2 (optional)', 'Does the halogen bond last?', 'Re-tests the step-2 criteria in every frame of an MD trajectory and reports occupancy, residence time and the best frame, which is passed on to the QM model.')}
    ${assumeBox('T2')}
    <div class="grid g2">
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Inputs</h2>
          <h3>Topology</h3><div class="row"><label class="chk"><input type="radio" name="topo" value="cur" checked> Current structure</label><label class="chk"><input type="radio" name="topo" value="file"> Upload topology PDB</label></div>
          <div class="drop" id="dropTopo" style="margin-top:6px">Topology <b>.pdb</b> (same atom order as trajectory)<input type="file" accept=".pdb,.pqr"></div>
          <h3>Trajectory</h3><div class="drop" id="dropTraj">Drop <b>.dcd</b> (CHARMM/NAMD/OpenMM), multi-model <b>.pdb</b> or multi-frame <b>.xyz</b><input type="file" accept=".dcd,.pdb,.xyz"></div>
          <p class="muted" style="font-size:11.5px;margin:6px 0 0">GROMACS .xtc / AMBER .nc: convert once with <span class="mono">gmx trjconv -o traj.pdb</span> or <span class="mono">cpptraj … trajout traj.dcd</span>.</p>
          <div style="margin-top:10px">${pblock('dyn', 'Settings', stagePaths('T2'), true)}</div>
          <div class="row" style="margin-top:6px"><button class="btn sm" id="btnTDemo">Generate demo trajectory (synthetic)</button></div>
          <div class="row" style="margin-top:10px"><button class="btn primary" id="btnTraj">Analyse trajectory</button><span id="tInfo" class="muted" style="font-size:12px"></span></div>
          <div class="progress"><i id="tProg"></i></div></div>
        <div class="card"><h2>Occupancy</h2><div id="tOccChart"></div></div>
      </div>
      <div class="grid" style="align-content:start">
        <div class="card"><div class="row" style="justify-content:space-between"><h2>Persistent halogen bonds</h2><div class="row"><button class="btn sm" id="btnTReg" disabled>Send best frames to registry →</button><button class="btn sm" id="btnTCsv" disabled>CSV</button></div></div>
          <div class="tablewrap" id="tTable"><p class="muted" style="padding:10px">No trajectory analysed yet.</p></div></div>
        <div class="card"><h2 id="tSerTitle">Time series</h2><p class="sub">Select a row. Gaps mean the contact failed the criteria in those frames.</p><div id="tSerD"></div><div id="tSerT"></div><div class="row"><button class="btn sm" id="btnTFrame" disabled>Download representative frame (PDB)</button><button class="btn sm" id="btnTView" disabled>View representative frame</button></div></div>
        <div class="vslot"></div>
      </div>
    </div>${nextBtn('qm', 'Next: build the QM model')}`,
    init() {
      let topoFile = null, trajFile = null;
      dropZone($('#dropTopo'), f => { topoFile = f[0]; $('#dropTopo').firstChild.textContent = 'Topology: ' + f[0].name + ' '; $$('input[name=topo]')[1].checked = true; });
      dropZone($('#dropTraj'), f => { trajFile = f[0]; $('#dropTraj').firstChild.textContent = 'Trajectory: ' + f[0].name + ' '; });
      $('#btnTDemo').onclick = () => {
        if (!ST.S || ST.name !== 'BondX_DEMO') return toast('Load the built-in demo complex first (step 1)');
        const S = ST.S; const f0 = S.frames[0]; let seed = 11; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
        const br = S.atoms.filter(a => a.resName === 'BBZ').map(a => a.i); const frames = [];
        for (let f = 0; f < 120; f++) { const fr = Float32Array.from(f0); for (let i = 0; i < fr.length; i++) fr[i] += rnd() * 0.3;
          const off = f > 60 && f < 95 ? Math.min(1.6, 0.12 * (f - 60)) : 0; br.forEach(i => { fr[3 * i + 2] += off; fr[3 * i] += off * 0.6; }); frames.push(fr); }
        ST.demoTraj = frames; trajFile = null; $('#dropTraj').firstChild.textContent = 'Trajectory: synthetic demo (120 frames) '; toast('Synthetic 120-frame trajectory ready — press Analyse');
      };
      $('#btnTraj').onclick = async () => {
        try {
          let T;
          if ($$('input[name=topo]')[1].checked) { if (!topoFile) return toast('Upload a topology PDB'); T = B.parsePDB(await readText(topoFile), { altloc: ST.P.struct.altloc, bondTol: ST.P.struct.bondTol }); } else { if (!ST.S) return toast('Load a structure first (or upload a topology)'); T = ST.S; }
          let frames, boxes, tname;
          if (trajFile) {
            tname = trajFile.name; const h = await sha256(await readBuf(trajFile));
            if (/\.dcd$/i.test(trajFile.name)) { const d = B.parseDCD(await readBuf(trajFile)); frames = B.dropVsites(T, d.frames); boxes = d.boxes; }
            else if (/\.xyz$/i.test(trajFile.name)) { frames = B.dropVsites(T, B.parseXYZFrames(await readText(trajFile), T.vsites ? T.vsites.fileAtoms : T.atoms.length).map(x => x.f)); boxes = null; }
            else { const M = B.parsePDB(await readText(trajFile), { altloc: ST.P.struct.altloc }); if (M.atoms.length !== T.atoms.length) throw new Error(`trajectory has ${M.atoms.length} atoms, topology ${T.atoms.length}`); frames = M.frames; boxes = M.boxes; }
            log('T2.load', `${tname}: ${frames.length} frames; sha256=${h}`);
          } else if (ST.demoTraj && T === ST.S) { frames = ST.demoTraj; boxes = null; tname = 'synthetic demo trajectory'; log('T2.load', 'SYNTHETIC demo trajectory (120 frames)'); }
          else if (T.frames.length > 1) { frames = T.frames; boxes = T.boxes; tname = 'models of ' + ST.name; }
          else return toast('Upload a trajectory (or load a multi-model PDB)');
          const stride = Math.max(1, Math.round(ST.P.dyn.stride) || 1);
          const acc = B.trajAccumulator(T, ST.P, { stride }); const t0 = performance.now();
          for (let f = 0; f < frames.length; f += stride) { acc.add(frames[f], boxes ? boxes[f] : null, f); if ((f / stride) % 50 === 0) { $('#tProg').style.width = (f / frames.length * 100) + '%'; await new Promise(r => setTimeout(r, 0)); } }
          const res = acc.finish(); $('#tProg').style.width = '100%';
          ST.traj = { T, frames, boxes, name: tname, res, dt: ST.P.dyn.dt || 1, stride };
          $('#tInfo').textContent = `${res.framesUsed} frames analysed in ${((performance.now() - t0) / 1000).toFixed(2)} s`;
          log('T2.analyse', `${tname}; ${res.framesUsed} frames (stride ${stride}); ${res.results.length} persistent XB (≥${ST.P.dyn.occCut}%); dyn=${JSON.stringify(ST.P.dyn)}`); this.render(); setDone('traj');
        } catch (e) { toast('Trajectory error: ' + e.message, 5000); }
      };
      $('#btnTCsv').onclick = () => download(`${ST.name}_XB_T2_dynamics.csv`, B.toCSV(ST.traj.res.results, [{ h: 'donor', k: 'donorLabel' }, { h: 'X', k: 'Xname' }, { h: 'acceptor', k: 'accLabel' }, { h: 'acc_atom', k: 'accName' }, { h: 'occupancy_pct', k: r => r.occupancy * 100 }, { h: 'mean_d', k: 'meanD' }, { h: 'sd_d', k: 'sdD' }, { h: 'mean_theta1', k: 'meanT1' }, { h: 'sd_theta1', k: 'sdT1' }, { h: 'best_frame', k: 'bestFrame' }, { h: 'best_d', k: 'bestD' }, { h: 'best_theta1', k: 'bestT1' }, { h: 'mean_residence_frames', k: 'meanResidence' }, { h: 'max_residence_frames', k: 'maxResidence' }, { h: 'events', k: 'nEvents' }, { h: 'D_score', k: 'D' }]));
      $('#btnTReg').onclick = () => { const sel = $$('#tTable .selt:checked').map(i => ST.traj.res.results[+i.value]); const T = ST.traj.T; const prep = B.prepareDetection(T, ST.P);
        sel.forEach(r => { const frame = r.bestFrame; const hit = B.detectFrame(T, prep, ST.traj.frames[frame], ST.traj.boxes ? ST.traj.boxes[frame] : null).find(h => h.key === r.key);
          if (!hit) return; B.decorate(T, hit); addCandidate(Object.assign(candFromHit(Object.assign({}, r, { d: hit.d, theta1: hit.theta1, theta2: hit.theta2, G: hit.G }), 'traj', 'T', frame), { occupancy: r.occupancy, Dscore: r.D, meanD: r.meanD })); });
        toast(`${sel.length} representative frame(s) in registry`); };
      $('#btnTFrame').onclick = () => { const r = this.sel; if (!r) return; download(`${ST.name}_frame${r.bestFrame}.pdb`, B.writePDB(ST.traj.T, ST.traj.frames[r.bestFrame], null, `BondX representative frame ${r.bestFrame} for ${r.donorLabel}@${r.Xname}...${r.accLabel}@${r.accName}`)); };
      $('#btnTView').onclick = () => { const r = this.sel; if (!r) return; const T = ST.traj.T; const tmp = Object.assign({}, T, { frames: [ST.traj.frames[r.bestFrame]] }); const hit = B.detectFrame(T, B.prepareDetection(T, ST.params || ST.P), ST.traj.frames[r.bestFrame], null).find(h => h.key === r.key); if (hit) B.decorate(T, hit); V.structure(tmp, 0, hit ? [hit] : [], hit || null); };
    },
    render() {
      const res = ST.traj.res.results; const dt = ST.traj.dt;
      $('#btnTCsv').disabled = !res.length; $('#btnTReg').disabled = !res.length;
      $('#tTable').innerHTML = res.length ? `<table class="t"><thead><tr><th><input type="checkbox" class="selAll" checked></th><th>Interaction</th><th>Occ.</th><th>⟨d⟩ ± σ (Å)</th><th>⟨θ₁⟩ ± σ</th><th>Best frame</th><th>τ (ps)</th><th>Events</th><th>D</th></tr></thead><tbody>${res.map((r, k) => `
        <tr data-k="${k}"><td><input type="checkbox" class="selt" value="${k}" checked></td><td>${elTag(r.halogen)} ${esc(r.donorLabel)}@${esc(r.Xname)} → ${esc(r.accLabel)} <b>${esc(r.accName)}</b></td><td class="n">${f2(r.occupancy * 100, 1)}%</td><td class="n">${f2(r.meanD)} ± ${f2(r.sdD)}</td><td class="n">${f2(r.meanT1, 1)} ± ${f2(r.sdT1, 1)}</td><td class="n">${r.bestFrame} <span class="muted">(${f2(r.bestD)} Å, ${f2(r.bestT1, 0)}°)</span></td><td class="n">${f2(r.meanResidence * dt, 1)}</td><td class="n">${r.nEvents}</td><td class="n">${f2(r.D)}</td></tr>`).join('')}</tbody></table>` : '<p class="muted" style="padding:10px">No halogen bond reaches the occupancy cut-off.</p>';
      bindSelectAll($('#tTable'));
      $$('#tTable tbody tr').forEach(tr => tr.addEventListener('click', e => { if (e.target.tagName === 'INPUT') return; $$('#tTable tr').forEach(t => t.classList.remove('sel')); tr.classList.add('sel'); this.series(res[+tr.dataset.k]); }));
      Charts.hbars($('#tOccChart'), { rows: res.map(r => ({ label: `${r.donorLabel.split(':')[0]}@${r.Xname}→${r.accLabel.split(':')[0]} ${r.accName}`, value: r.occupancy * 100, text: f2(r.occupancy * 100, 1) + '%', color: HCOLOR[r.halogen], tip: `<b>${esc(r.donorLabel)}@${esc(r.Xname)} → ${esc(r.accLabel)} ${esc(r.accName)}</b><br>occupancy ${f2(r.occupancy * 100, 1)}% · τ ${f2(r.meanResidence * dt, 1)} ps` })), max: 100, xLabel: 'occupancy (% of analysed frames)', fmt: v => v,
        legend: ['I', 'BR', 'CL'].filter(x => res.some(r => r.halogen === x)).map(x => ({ name: HN(x), color: HCOLOR[x] })) });
      if (res.length) this.series(res[0]);
    },
    series(r) {
      this.sel = r; $('#btnTFrame').disabled = false; $('#btnTView').disabled = false; const dt = ST.traj.dt;
      $('#tSerTitle').textContent = `Time series — ${r.donorLabel}@${r.Xname} → ${r.accLabel} ${r.accName}`;
      const xs = r.series.frames.map(f => f * dt); const cut = (ST.P.geom.cutoffs[r.halogen] || {})[r.type === 'pi' ? 'AR' : r.accEl];
      const span = [0, (ST.traj.frames.length - 1) * dt];
      Charts.line($('#tSerD'), { series: [{ name: 'd(X···Y)', color: HCOLOR[r.halogen], xs, ys: r.series.d }], xLabel: 'time (ps)', yLabel: 'd (Å)', refY: cut, refLabel: `cut-off ${f2(cut)} Å`, height: 190, gapAbove: ST.traj.stride * dt * 1.01, xDomain: span });
      Charts.line($('#tSerT'), { series: [{ name: 'θ₁', color: 'var(--s7)', xs, ys: r.series.t1 }], xLabel: 'time (ps)', yLabel: 'θ₁ (°)', refY: ST.P.geom.theta1Min, refLabel: `θ₁ min ${ST.P.geom.theta1Min}°`, height: 190, gapAbove: ST.traj.stride * dt * 1.01, digits: 1, xDomain: span });
    }
  };

  /* ====================== MATERIALS ====================== */
  R.material = {
    html: () => `${head('Tool · T1', 'Materials and crystal structures', 'Reads a CIF (expanded by its symmetry operators, with periodic images) or an XYZ cluster, finds C–X···Y halogen bonds and builds a two-molecule QM model for any contact. That model then goes through the same NBO → ALMO-EDA → BXQI steps.')}
    ${assumeBox('T1')}
    <div class="grid g2">
      <div class="card" style="align-self:start"><h2>Input</h2><div class="drop" id="dropCif">Drop <b>.cif</b> or <b>.xyz</b><input type="file" accept=".cif,.xyz"></div><div id="mInfo" style="margin-top:8px"></div>
        <div style="margin-top:10px">${pblock('mat', 'Crystal settings', ['mat.images', 'mat.minDist'], true)}</div>
        <div class="row" style="margin-top:10px"><button class="btn primary" id="btnMat" disabled>Find halogen bonds</button></div></div>
      <div class="grid" style="align-content:start">
        <div class="vslot"></div>
        <div class="card"><div class="row" style="justify-content:space-between"><h2>Contacts</h2><button class="btn sm" id="btnMReg" disabled>Send selected to registry →</button></div><div class="tablewrap" id="mTable"><p class="muted" style="padding:10px">No material loaded.</p></div></div>
      </div>
    </div>`,
    init() {
      dropZone($('#dropCif'), async f => {
        try { const txt = await readText(f[0]); const M = /\.cif$/i.test(f[0].name) ? B.parseCIF(txt) : B.parseXYZ(txt); ST.mat = { M, name: f[0].name.replace(/\.\w+$/, ''), hits: [] };
          $('#mInfo').innerHTML = `<div class="note">${esc(f[0].name)}: ${M.atoms.length} atoms${M.cell ? ` in the unit cell (${M.asymCount} asymmetric, ${M.ops.length} symmetry operators; a=${M.cell.a}, b=${M.cell.b}, c=${M.cell.c} Å)` : ''}.</div>`;
          $('#btnMat').disabled = false; log('material.load', `${f[0].name}; ${M.atoms.length} atoms; sha256=${await sha256(txt)}`); V.material(M);
        } catch (e) { $('#mInfo').innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
      });
      $('#btnMat').onclick = () => {
        const hits = B.detectMaterialXB(ST.mat.M, ST.P); ST.mat.hits = hits; log('T1.material', `${hits.length} contacts`);
        $('#btnMReg').disabled = !hits.length;
        $('#mTable').innerHTML = hits.length ? `<table class="t"><thead><tr><th><input type="checkbox" class="selAll"></th><th>A–X</th><th>Y (image)</th><th>d (Å)</th><th>R</th><th>θ₁</th><th>G</th><th title="Symmetry-equivalent copies in the unit cell">Copies</th><th>Status</th></tr></thead><tbody>${hits.map((h, k) => `<tr data-k="${k}"><td><input type="checkbox" class="selm" value="${k}" ${h.status === 'XB' ? 'checked' : ''}></td><td>${elTag(h.halogen)} ${esc(h.Alabel)}–${esc(h.Xlabel)}</td><td>${esc(h.Ylabel)} <span class="muted mono">[${esc(h.image)}]</span></td><td class="n">${f2(h.d, 3)}</td><td class="n">${f2(h.R, 3)}</td><td class="n">${f2(h.theta1, 1)}</td><td class="n">${f2(h.G)}</td><td class="n">${h.mult || 1}</td><td>${h.status === 'XB' ? '<span class="pill ok">XB</span>' : '<span class="pill mute">near-miss</span>'}</td></tr>`).join('')}</tbody></table>` : '<p class="muted" style="padding:10px">No halogen bonds found.</p>';
        bindSelectAll($('#mTable')); $$('#mTable tbody tr').forEach(tr => tr.addEventListener('click', e => { if (e.target.tagName !== 'INPUT') V.material(ST.mat.M, hits[+tr.dataset.k]); }));
        setDone('material');
      };
      $('#btnMReg').onclick = () => { const sel = $$('#mTable .selm:checked').map(i => ST.mat.hits[+i.value]);
        sel.forEach(h => addCandidate({ sourceKey: `cif|${ST.mat.name}|${h.X}>${h.Y}|${h.image}`, source: 'crystal ' + ST.mat.name, structName: ST.mat.name, structRef: 'M', frame: 0, mhit: h, mat: ST.mat,
          label: `${h.Alabel}–${h.Xlabel} ··· ${h.Ylabel} [${h.image}]`, halogen: h.halogen, accClass: 'crystal ' + h.accEl, d: h.d, theta1: h.theta1, G: h.G }));
        toast(`${sel.length} candidate(s) in registry`); };
    }
  };

  /* ====================== 4 QM MODEL ====================== */
  const QM_REGION = ['qm.mode', 'qm.shell', 'qm.shellWater', 'qm.shellIons', 'qm.forceH', 'qm.linkScheme', 'qm.xh'];
  const QM_FF = ['ff.name', 'ff.water', 'ff.hisDefault', 'ff.termini', 'ff.chainBreaks', 'ff.ssCut', 'ff.chargeSource', 'ff.foldH', 'ff.missing', 'qm.ligandCharge'];
  const QM_CLUSTER = ['qm.pcm', 'qm.eps'];
  const QM_MM = ['mm.radius', 'mm.selection', 'mm.boundary', 'mm.neutralizeCut'];
  function progNote() {
    const q = ST.P.qm; const t = [];
    if (q.program === 'orca') t.push(B.ORCA_NBO_CAUTION);
    if (q.program === 'qchem') t.push(B.QCHEM_NBO_CAUTION);
    if (q.program === 'gaussian' && q.nbo7) t.push(B.G16_NBO7_CAUTION);
    return t.map(x => `<div class="note warn" style="margin-top:8px"><b>Caution.</b> ${esc(x)}</div>`).join('');
  }
  const QM_PROG = ['qm.program', 'qm.method', 'qm.basis', 'qm.nproc', 'qm.mem', 'qm.nbo7'];
  const METHODS = { gaussian: { method: 'wB97XD', basis: 'def2TZVP' }, orca: { method: 'wB97X-D3', basis: 'def2-TZVP' }, qchem: { method: 'wB97X-D', basis: 'def2-TZVP' } };
  R.qm = {
    html: () => `${head('Step 4 · T3a', 'Build the QM model with force-field protonation and charges', 'Cuts the halogenated ligand and its partners out of the structure, protonates residues from force-field templates and ligands from one Kekulé structure (checked against the PDB formula), caps cut bonds and writes Gaussian, ORCA or Q-Chem inputs. In QM/MM mode the rest of the protein enters as force-field point charges. Every QM atom is traced back to the PDB in a manifest.')}
    ${assumeBox('T3')}
    <div class="card" style="margin-bottom:16px"><div class="row" style="justify-content:space-between"><h2>1 · Choose candidates</h2><span class="muted" style="font-size:12px">Candidates sharing one donor ligand and frame are modelled together.</span></div><div id="qmReg" style="margin-top:8px"></div></div>
    <div class="grid g2">
      <div class="grid" style="align-content:start">
        <div class="card"><h2>2 · Environment</h2>
          <div class="seg" id="qEnvSeg"><button data-env="cluster">QM cluster</button><button data-env="qmmm">QM/MM (point charges)</button></div>
          <div id="qEnvParams" style="margin-top:10px"></div>
        </div>
        <div class="card"><h2>3 · Protonation & charges</h2><p class="sub">Charges and protonation templates come from the selected force field; the same choice is used for the QM fragments and for the MM point charges.</p>
          ${pblock('qmff', 'Force-field settings', QM_FF, true)}
          <h3>Ligand charges (optional)</h3><div class="drop" id="dropMol2">Drop a <b>.mol2</b> with partial charges (e.g. AM1-BCC from antechamber) for ligands / cofactors<input type="file" accept=".mol2"></div><div id="mol2Msg" style="margin-top:6px"></div>
          <h3>Protonation states near the candidates</h3><div id="protTable"><p class="muted">Select candidates above to list nearby titratable residues.</p></div>
        </div>
        <div class="card"><h2>4 · QM region & program</h2>
          ${pblock('qmreg', 'QM region', QM_REGION)}
          <div style="margin-top:8px">${pblock('qmprog', 'Program & level of theory (NBO single point)', QM_PROG, true)}</div>
          <div id="qProgNote">${progNote()}</div>
          <h3>Manual overrides (optional)</h3>
          <div class="row"><label class="f">Donor fragment charge<input type="text" id="qQ1" placeholder="auto" style="width:70px"></label><label class="f">Acceptor fragment charge<input type="text" id="qQ2" placeholder="auto" style="width:70px"></label><label class="f">Multiplicity<input type="text" id="qMult" placeholder="auto" style="width:70px"></label></div>
          <div class="row" style="margin-top:14px"><button class="btn primary" id="btnBuild">Build QM model for selected</button></div>
          <div id="qMsg" style="margin-top:8px"></div>
        </div>
      </div>
      <div class="grid" style="align-content:start">
        <div class="vslot"></div>
        <div class="card" id="qOut" style="display:none"><div class="row" style="justify-content:space-between"><h2 id="qOutTitle">Model</h2><div class="row"><button class="btn sm primary" id="btnZip">Download job bundle (.zip)</button></div></div>
          <div id="qKpis" class="kpis"></div><div id="qWarn"></div>
          <h3>Fragment charges</h3><div class="tablewrap" id="qCharges" style="max-height:220px"></div>
          <div id="qEmb"></div>
          <h3>Files</h3><div class="tabs2" id="qTabs"></div><pre class="code" id="qCode"></pre></div>
      </div>
    </div>${nextBtn('nbo', 'Next: NBO gate')}`,
    init() {
      $$('#qEnvSeg button').forEach(b => b.onclick = () => { ST.P.qm.env = b.dataset.env; onParamChanged('qm.env'); });
      $('#btnBuild').onclick = () => this.build();
      $('#btnZip').onclick = () => { const m = ST.models[this.cur]; if (m) { download(`${m.manifest.id}_bundle.zip`, zip(m.files)); log('T3.bundle', m.manifest.id); } };
      dropZone($('#dropMol2'), async f => {
        try { const txt = await readText(f[0]); ST.mol2 = B.parseMOL2(txt, true); ST.mol2Name = f[0].name;
          $('#mol2Msg').innerHTML = `<div class="note ok">${esc(f[0].name)}: ${Object.entries(ST.mol2).map(([k, v]) => `${esc(k)} (${v._nAtoms} atoms, total ${f2(v._total, 3)} e)`).join(', ')}. Charges are matched to residues of the same name by atom name.</div>`;
          log('mol2.load', `${f[0].name}; residues ${Object.keys(ST.mol2).join(',')}; sha256=${await sha256(txt)}`); this.renderProt();
        } catch (e) { $('#mol2Msg').innerHTML = `<div class="note bad">${esc(e.message)}</div>`; }
      });
      document.addEventListener('change', e => { if (e.target.classList && e.target.classList.contains('selc') && current === 'qm') this.renderProt(); });
    },
    refreshEnv() {
      $$('#qEnvSeg button').forEach(b => b.classList.toggle('on', b.dataset.env === ST.P.qm.env));
      $('#qEnvParams').innerHTML = ST.P.qm.env === 'qmmm' ? pblock('qmmm', 'QM/MM embedding', QM_MM, true) + '<p class="muted" style="font-size:12px;margin:8px 0 0">Charges: PQR file if loaded, else the force field below. QM–MM van der Waals terms are not included (single-point NBO).</p>' : pblock('qmcl', 'Cluster environment', QM_CLUSTER, true);
      this.renderProt();
    },
    onShow() { $('#qmReg').innerHTML = registryTable({ select: true, isChecked: c => c.state === 'DETECTED' }); bindSelectAll($('#qmReg')); this.refreshEnv(); },
    selected() { return $$('#qmReg .selc:checked').map(i => ST.cands.find(c => c.cid === i.value)).filter(Boolean); },
    renderProt() {
      const el = $('#protTable'); if (!el) return; const S = ST.S; const sel = this.selected().filter(c => c.structRef !== 'M');
      if (!S || !sel.length) { el.innerHTML = '<p class="muted">Select candidates above to list nearby titratable residues.</p>'; return; }
      const T = sel[0].structRef === 'T' && ST.traj ? ST.traj.T : S; const frame = T.frames[0];
      const centres = sel.map(c => B.coordOf(frame, c.hit.X)); const rad = ST.P.qm.env === 'qmmm' ? ST.P.mm.radius : 8;
      const ffa = ffAssign(T); const list = B.titratableResidues(T, ST.P, ffa).filter(t => T.residues[t.idx].atoms.some(i => centres.some(p => B.dist(p, B.coordOf(frame, i)) <= rad)));
      if (!list.length) { el.innerHTML = `<p class="muted">No titratable residue within ${rad} Å.</p>`; return; }
      el.innerHTML = `<p class="sub">Residues with an atom within ${rad} Å of the selected halogens. Change a state to override the default rules (recorded in the manifest).</p><div class="tablewrap" style="max-height:240px"><table class="t"><thead><tr><th>Residue</th><th>State</th><th>Template</th><th>Charge</th><th>Why</th></tr></thead><tbody>${list.map(t => `<tr><td class="mono">${esc(t.resName)} ${esc(t.id)}</td><td>${t.choices.length > 1 ? `<select data-prot="${esc(t.id)}"><option value="">auto (${esc(t.variant)})</option>${t.choices.map(c => `<option value="${c}" ${ST.P.ff.overrides[t.id] === c ? 'selected' : ''}>${c}</option>`).join('')}</select>` : esc(t.variant)}</td><td class="mono">${esc(t.tpl || '—')}</td><td class="n">${f2(t.net, 3)}</td><td class="muted" style="font-size:11.5px">${esc(t.why)}</td></tr>`).join('')}</tbody></table></div>`;
      $$('[data-prot]', el).forEach(s => s.onchange = () => { if (s.value) ST.P.ff.overrides[s.dataset.prot] = s.value; else delete ST.P.ff.overrides[s.dataset.prot]; log('protonation', `${s.dataset.prot} → ${s.value || 'auto'}`); setTimeout(() => { refreshAssumptions(); this.renderProt(); }, 0); });
    },
    build() {
      const sel = this.selected(); if (!sel.length) return toast('Select candidates in the registry');
      const groups = {}; sel.forEach(c => { const k = c.structRef + '|' + c.structName + '|' + c.frame + '|' + (c.structRef === 'M' ? c.cid : c.hit.donorRes); (groups[k] = groups[k] || []).push(c); });
      $('#qMsg').innerHTML = Object.keys(groups).length > 1 ? `<div class="note">Selected candidates come from ${Object.keys(groups).length} different ligands/structures/frames — building ${Object.keys(groups).length} models (one halogen-bond donor fragment per ALMO-EDA).</div>` : '';
      const P = ST.P; const prog = P.qm.program; const memGB = parseFloat(P.qm.mem) || 16;
      const o = { nproc: P.qm.nproc, mem: P.qm.mem, maxcore: Math.floor(memGB * 1000 * 0.75 / Math.max(1, P.qm.nproc)) || 2000, memMB: Math.round(memGB * 1000), method: P.qm.method.trim(), basis: P.qm.basis.trim(), solvent: P.qm.pcm, eps: P.qm.eps, nbo7: P.qm.nbo7 };
      const qmmm = P.qm.env === 'qmmm'; let lastId = null;
      const chg = {}; if ($('#qQ1').value.trim() !== '') chg[1] = parseInt($('#qQ1').value); if ($('#qQ2').value.trim() !== '') chg[2] = parseInt($('#qQ2').value);
      const mult = $('#qMult').value.trim() ? parseInt($('#qMult').value) : undefined;
      for (const cs of Object.values(groups)) {
        const c0 = cs[0]; let model, pcs = null;
        try {
          if (c0.structRef === 'M') { model = B.buildMaterialQM(c0.mat.M, Object.assign({ cid: c0.cid }, c0.mhit), { chargeOverride: chg, multiplicity: mult, sourceName: c0.structName, P }); if (qmmm) toast('QM/MM is not available for crystal dimers; cluster model written.'); }
          else {
            const S = c0.structRef === 'T' ? ST.traj && ST.traj.T : ST.S; if (!S) { toast('Structure for ' + c0.cid + ' is no longer loaded'); continue; }
            const hits = cs.map(c => Object.assign({ cid: c.cid }, c.hit));
            const Sx = c0.structRef === 'T' ? Object.assign({}, S, { frames: [ST.traj.frames[c0.frame]] }) : S;
            const mol2 = mol2For(S); const ffa = B.assignForceField(Sx, P, { mol2 });
            model = B.buildQMModel(Sx, 0, hits, { P, ffa, mol2, chargeOverride: Object.keys(chg).length ? chg : null, multiplicity: mult, sourceName: c0.structName });
            model.manifest.frame = c0.frame;
            if (qmmm) {
              pcs = B.embeddingCharges(Sx, 0, model.include, { P, ffa, qmResidueCharges: model.qmResCharge });
              model.manifest.embedding = { source: pcs.source, water: pcs.water, n: pcs.charges.length, net: +pcs.net.toFixed(6), radius: pcs.radius, selection: pcs.selection, boundary: pcs.boundary, M1: pcs.M1, cutResidues: pcs.cut, foldedH: pcs.foldedH, notes: pcs.notes };
              pcs.notes.forEach(n => model.manifest.warnings.push(n));
            }
          }
        } catch (e) { $('#qMsg').innerHTML += `<div class="note bad">${esc(c0.cid)}: ${esc(e.message)}</div>`; continue; }
        const mf = model.manifest; const id = mf.id; const job = id.replace('BX-', 'bx_');
        mf.env = qmmm && pcs ? `QM/MM electrostatic embedding (${pcs.source}, ${pcs.charges.length} charges, ${pcs.boundary})` : (o.solvent ? `QM cluster, implicit ε=${o.eps}` : 'QM cluster, gas phase');
        mf.program = prog; mf.level = `${o.method}/${o.basis}`; mf.assumptions = B.assumptions('T3', P); mf.parameters = { qm: P.qm, mm: P.mm, ff: P.ff, struct: P.struct };
        const files = {};
        if (prog === 'gaussian') files[`${job}_nbo.com`] = B.writeGaussian(model, Object.assign({}, o, { job: job + '_nbo' }), pcs);
        if (prog === 'orca') { const w = B.writeORCA(model, Object.assign({}, o, { job: job + '_nbo' }), pcs); files[`${job}_nbo.inp`] = w.inp; if (w.pc) files[`${job}_nbo.pc`] = w.pc; }
        if (prog === 'qchem') files[`${job}_nbo.in`] = B.writeQChemNBO(model, o, pcs);
        if (prog === 'orca') mf.warnings.push(B.ORCA_NBO_CAUTION);
        if (prog === 'gaussian' && o.nbo7) mf.warnings.push(B.G16_NBO7_CAUTION);
        files[`${job}_cluster.xyz`] = B.writeXYZ(model);
        if (pcs) files[`${job}_mm_charges.csv`] = B.chargesCSV(pcs);
        files[`${job}_manifest.json`] = JSON.stringify(mf, null, 1);
        files['README.txt'] = `BondX job bundle ${id}\nRun ${ST.runId}\n\n1. Submit the *_nbo input (${prog}). NBO second-order perturbation output is required.\n2. Return to BondX → NBO gate, select model ${id} (or load ${job}_manifest.json) and open the output file.\n3. Only interactions passing the gate go on to ALMO-EDA (Q-Chem), whose input is written on the ALMO-EDA page.\n\nAtom order in all inputs equals manifest order. Do not re-order atoms or allow symmetry re-orientation.\n\nAssumptions\n${mf.assumptions.map(a => ' - ' + a).join('\n')}\n\nWarnings\n${mf.warnings.length ? mf.warnings.map(w => ' - ' + w).join('\n') : ' - none'}\n`;
        ST.models[id] = { manifest: mf, qm: model.qm, files, candIds: cs.map(c => c.cid), created: new Date().toISOString(), pcs };
        cs.forEach(c => { c.modelId = id; if (c.state === 'DETECTED') c.state = 'MODELLED'; });
        log('T3.model', `${id}: ${mf.nAtoms} atoms (${mf.addedH} H added, ${mf.linkH} link), charges ${JSON.stringify(mf.fragCharge)}, mult ${mf.multiplicity}, ${mf.env}, ${prog} ${o.method}/${o.basis}; ff ${P.ff.name}; candidates ${cs.map(c => c.cid).join(',')}`);
        lastId = id;
      }
      if (lastId) { this.show(lastId); $('#qmReg').innerHTML = registryTable({ select: true, isChecked: c => c.state === 'DETECTED' }); bindSelectAll($('#qmReg')); setDone('qm'); }
    },
    show(id) {
      this.cur = id; const m = ST.models[id]; const mf = m.manifest; $('#qOut').style.display = '';
      $('#qOutTitle').innerHTML = `Model <span class="mono">${id}</span>`;
      $('#qKpis').innerHTML = [[mf.nAtoms, 'QM atoms'], [mf.addedH, 'H added'], [mf.linkH, 'link atoms'], [`${mf.fragCharge[1] ?? 0} / ${mf.fragCharge[2] ?? 0}`, 'charge donor / acceptor'], [mf.multiplicity, 'multiplicity'], [mf.embedding ? mf.embedding.n : '—', 'MM point charges']].map(k => `<div class="kpi"><b>${k[0]}</b><span>${k[1]}</span></div>`).join('');
      $('#qWarn').innerHTML = mf.warnings.map(w => `<div class="note warn" style="margin-bottom:6px">${esc(w)}</div>`).join('') + `<div class="note" style="margin-bottom:8px">${esc(mf.env)} · ${mf.xbs.map(x => `${esc(x.id)}: X=${x.X}, A=${x.A}, ${x.Y ? 'Y=' + x.Y : 'ring=' + x.ring.join(',')}`).join(' · ')}</div>`;
      $('#qCharges').innerHTML = mf.residueCharges ? `<table class="t"><thead><tr><th>Residue</th><th>Fragment</th><th>Template</th><th>Charge</th><th>Source</th></tr></thead><tbody>${mf.residueCharges.map(r => `<tr><td class="mono">${esc(r.resName)} ${esc(r.res)}</td><td class="n">${r.frag}</td><td class="mono">${esc(r.tpl || '—')}</td><td class="n">${r.charge}</td><td class="muted" style="font-size:11.5px">${esc(r.how)}</td></tr>`).join('')}</tbody></table>` : '<p class="muted" style="padding:8px">Crystal dimer: both fragments neutral unless overridden.</p>';
      const e = mf.embedding; $('#qEmb').innerHTML = e ? `<h3>Embedding</h3><div class="kpis">${[[e.n, 'point charges'], [f2(e.net, 3), 'net MM charge (e)'], [e.radius + ' Å', 'radius'], [e.boundary, 'boundary'], [e.M1, 'M1 atoms'], [e.cutResidues, 'cut residues']].map(k => `<div class="kpi"><b>${k[0]}</b><span>${k[1]}</span></div>`).join('')}</div><p class="muted" style="font-size:12px">Source: ${esc(e.source)}${e.source !== 'PQR' ? ', water ' + esc(e.water) : ''}${e.foldedH ? '; missing H charges folded onto heavy atoms' : ''}. Full list in <span class="mono">*_mm_charges.csv</span>.</p>` : '';
      const names = Object.keys(m.files); $('#qTabs').innerHTML = names.map((n, k) => `<button data-n="${esc(n)}" class="${k === 0 ? 'on' : ''}">${esc(n)}</button>`).join('');
      const showF = n => { const t = m.files[n]; $('#qCode').textContent = t.length > 200000 ? t.slice(0, 200000) + '\n… (truncated in preview; full file in the bundle)' : t; $$('#qTabs button').forEach(b => b.classList.toggle('on', b.dataset.n === n)); };
      $$('#qTabs button').forEach(b => b.onclick = () => showF(b.dataset.n)); showF(names[0]);
      V.qmModel({ qm: m.qm, manifest: mf }, m.pcs);
    }
  };
  // program change: update method/basis defaults for that program
  document.addEventListener('change', e => { if (e.target.dataset && e.target.dataset.p === 'qm.program') { const m = METHODS[ST.P.qm.program]; ST.P.qm.method = m.method; ST.P.qm.basis = m.basis; log('param', `qm.method/basis → ${m.method}/${m.basis} (defaults for ${ST.P.qm.program})`); }
    if (e.target.dataset && /^qm\.(program|nbo7)$/.test(e.target.dataset.p || '')) setTimeout(() => { const n = $('#qProgNote'); if (n) n.innerHTML = progNote(); }, 0); });

  /* ====================== 5 NBO GATE ====================== */
  R.nbo = {
    html: () => `${head('Step 5 · T3b', 'NBO gate: is there real charge transfer?', 'Reads the NBO second-order perturbation table and adds up donation into the halogen σ-hole. Contacts that reach the threshold are <b>electronically verified</b> and go on to ALMO-EDA; the rest are rejected as geometric false positives.')}
    ${assumeBox('T3b')}
    <div class="grid g2">
      <div class="card" style="align-self:start"><h2>Inputs</h2>
        <label class="f">QM model<select id="nModel"></select></label>
        <div class="row" style="margin-top:6px"><div class="drop" id="dropMan" style="flex:1">…or load a <b>manifest.json</b><input type="file" accept=".json"></div></div>
        <div class="row" style="margin-top:8px"><div class="drop" id="dropNbo" style="flex:1">Drop the QM output (<b>.log / .out</b>) with NBO analysis<input type="file" accept=".log,.out,.txt,.nbo"></div></div>
        <div style="margin-top:10px">${pblock('nbo', 'Gate settings', stagePaths('T3b'), true)}</div>
        <div id="nMsg"></div>
        <details style="margin-top:10px"><summary>Testing aid</summary><p class="muted" style="font-size:12px">Writes a <b>synthetic</b> Gaussian-format NBO output for the selected model so the gate can be tried without a QM job. Never use it for real conclusions.</p><button class="btn sm" id="btnFakeNbo">Download synthetic NBO output</button></details></div>
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Gate decisions</h2><div id="nRes"><p class="muted">No NBO output analysed yet.</p></div></div>
        <div class="card"><h2>ΣE(2) per interaction</h2><div id="nChart"></div></div>
        <div class="card"><h2>Contributing NBO terms</h2><div class="tablewrap" id="nRows"><p class="muted" style="padding:10px">—</p></div></div>
      </div>
    </div>${nextBtn('eda', 'Next: ALMO-EDA')}`,
    init() {
      dropZone($('#dropMan'), async f => { try { const mf = JSON.parse(await readText(f[0])); if (!mf.atoms || !mf.xbs) throw new Error('not a BondX manifest'); if (!ST.models[mf.id]) ST.models[mf.id] = { manifest: mf, qm: mf.atoms.map(a => ({ el: a.el, x: a.x, y: a.y, z: a.z, frag: a.frag })), files: {}, candIds: [] }; this.onShow(); $('#nModel').value = mf.id; toast('Manifest ' + mf.id + ' loaded'); } catch (e) { toast('Manifest error: ' + e.message); } });
      dropZone($('#dropNbo'), async f => this.run(f[0]));
      $('#btnFakeNbo').onclick = () => { const m = ST.models[$('#nModel').value]; if (!m) return toast('Select a model'); download(m.manifest.id + '_SYNTHETIC_nbo.log', syntheticNBO(m.manifest)); };
    },
    onShow() { const ids = Object.keys(ST.models); const cur = $('#nModel').value; $('#nModel').innerHTML = ids.length ? ids.map(id => `<option value="${id}">${id} — ${esc(ST.models[id].manifest.xbs.map(x => x.label).join(' ; ').slice(0, 70))}</option>`).join('') : '<option value="">(build a QM model first or load a manifest)</option>'; if (cur && ST.models[cur]) $('#nModel').value = cur; },
    async run(file) {
      const id = $('#nModel').value; const m = ST.models[id]; if (!m) return toast('Select a QM model / manifest first');
      const txt = await readText(file); const h = await sha256(txt); const nbo = B.parseNBO(txt); const NP = ST.P.nbo;
      if (nbo.error) { $('#nMsg').innerHTML = `<div class="note bad">${esc(nbo.error)}</div>`; return; }
      const ver = B.verifyGeometry(m.manifest, nbo.geometry, NP.geomTol); const gate = B.nboGate(m.manifest, nbo, { threshold: NP.threshold, includeRY: NP.includeRY }); const fq = B.fragmentNPA(m.manifest, nbo);
      m.nbo = { file: file.name, sha256: h, program: nbo.program, terminated: nbo.terminated, verification: ver, gate, fragNPA: fq, nEntries: nbo.entries.length, settings: NP };
      const blocked = ver.checked && ver.ok === false && ver.elMis;
      $('#nMsg').innerHTML = `<div class="note ${blocked ? 'bad' : ver.ok ? 'ok' : 'warn'}" style="margin-top:8px">${esc(nbo.program)} output · ${nbo.entries.length} E(2) entries · ${nbo.terminated ? 'normal termination' : '<b>no normal-termination marker</b>'}${nbo.threshold != null ? ` · NBO print threshold ${nbo.threshold} kcal/mol` : ''}<br>${esc(ver.note)}</div>${nbo.program === 'ORCA' ? `<div class="note warn" style="margin-top:6px">${esc(B.ORCA_NBO_CAUTION)} Check the parsed E(2) entries against the output by hand.</div>` : ''}`;
      if (blocked) { $('#nRes').innerHTML = '<div class="note bad">The output does not match this manifest — gate not applied.</div>'; log('T3.nbo', `${id}: ${file.name} REFUSED (${ver.note})`); return; }
      gate.forEach(g => { const c = ST.cands.find(c => c.cid === g.id) || ST.cands.find(c => c.modelId === id && c.label === g.label); if (!c) return; c.E2 = g.E2_core; c.nbo = { contrib: g.contrib, npaX: g.npaX }; c.verified = g.verified; c.state = g.verified ? 'NBO-VERIFIED' : 'NBO-REJECTED'; });
      log('T3.nbo', `${id}: ${file.name} (${nbo.program}, sha256=${h}); ${ver.note}; threshold ${NP.threshold}, RY* ${NP.includeRY}; ${gate.map(g => `${g.id}=${g.E2_core.toFixed(2)} ${g.verified ? 'PASS' : 'FAIL'}`).join(', ')}`);
      $('#nRes').innerHTML = `<table class="t"><thead><tr><th>ID</th><th>Interaction</th><th>n→σ*</th><th>n→RY*</th><th>π→σ*</th><th>ΣE(2)</th><th>q<sub>NPA</sub>(X)</th><th>Decision</th></tr></thead><tbody>${gate.map(g => `<tr><td class="mono">${esc(g.id)}</td><td>${esc(g.label)}</td><td class="n">${f2(g.contrib.lp_sigma)}</td><td class="n">${f2(g.contrib.lp_ry)}${g.includeRY ? '' : ' <span class="muted">(not counted)</span>'}</td><td class="n">${f2(g.contrib.pi_sigma)}</td><td class="n"><b>${f2(g.E2_core)}</b></td><td class="n">${f2(g.npaX, 3)}</td><td>${g.verified ? '<span class="pill ok">VERIFIED → EDA</span>' : '<span class="pill bad">REJECTED</span>'}</td></tr>`).join('')}</tbody></table>
        ${fq ? `<div class="note" style="margin-top:8px">NPA inter-fragment charge transfer: the donor fragment gains <b>${f2(fq.ct, 4)} e</b> relative to its formal charge (fragment charges ${f2(fq.frag1, 3)} / ${f2(fq.frag2, 3)}).</div>` : ''}`;
      Charts.hbars($('#nChart'), { rows: gate.map(g => ({ label: g.id + ' ' + (g.label.split('···')[1] || ''), value: g.E2_core, color: g.verified ? 'var(--s1)' : 'var(--neutral)', tip: `<b>${esc(g.label)}</b><br>ΣE(2) ${f2(g.E2_core)} kcal/mol · ${g.verified ? 'verified' : 'rejected'}` })), threshold: NP.threshold, thresholdLabel: `gate ${NP.threshold} kcal/mol`, xLabel: 'ΣE(2) (kcal/mol)', fmt: v => v });
      const rows = gate.flatMap(g => g.rows.map(r => Object.assign({ id: g.id }, r)));
      const side = s => `${s.type}(${s.idx}) ${s.atoms.map(a => a.el[0] + a.el.slice(1).toLowerCase() + a.n).join('–')}`;
      $('#nRows').innerHTML = rows.length ? `<table class="t"><thead><tr><th>XB</th><th>Donor NBO</th><th>Acceptor NBO</th><th>Term</th><th>E(2)</th><th>ΔE (a.u.)</th><th>F(i,j)</th><th>Spin</th></tr></thead><tbody>${rows.map(r => `<tr><td class="mono">${esc(r.id)}</td><td class="mono">${esc(side(r.donor))}</td><td class="mono">${esc(side(r.acceptor))}</td><td>${esc(r.tag)}</td><td class="n">${f2(r.e2)}</td><td class="n">${f2(r.dE)}</td><td class="n">${f2(r.F, 3)}</td><td>${esc(r.spin)}</td></tr>`).join('')}</tbody></table>` : '<p class="muted" style="padding:10px">No NBO term couples the acceptor into the halogen σ-hole above the print threshold.</p>';
      setDone('nbo');
    }
  };

  /* ====================== 6 ALMO-EDA ====================== */
  const EDA_COMPS = [{ key: 'elec', name: 'Electrostatics', color: 'var(--s1)' }, { key: 'pauli', name: 'Pauli repulsion', color: 'var(--neutral)' }, { key: 'disp', name: 'Dispersion', color: 'var(--s3)' }, { key: 'pol', name: 'Polarization', color: 'var(--s2)' }, { key: 'ct', name: 'Charge transfer', color: 'var(--s7)' }];
  R.eda = {
    html: () => `${head('Step 6 · T4', 'ALMO-EDA: polarization or charge transfer?', 'Writes a Q-Chem EDA2 input only for models that passed the NBO gate (fragment 1 = halogenated ligand, fragment 2 = acceptor environment) and reads the energy decomposition.')}
    ${assumeBox('T4')}
    <div class="grid g2">
      <div class="card" style="align-self:start"><h2>Q-Chem EDA2 input</h2>
        <label class="f">Model<select id="eModel"></select></label>
        <div style="margin-top:10px">${pblock('eda', 'EDA settings', stagePaths('T4'), true)}</div>
        <div class="row" style="margin-top:8px"><label class="chk"><input type="checkbox" id="eOverride"> allow EDA for NBO-rejected/unverified models (logged as gate override)</label></div>
        <div class="row" style="margin-top:10px"><button class="btn primary" id="btnEdaIn">Generate EDA input</button></div>
        <div id="eGateMsg" style="margin-top:8px"></div>
        <pre class="code" id="eCode" style="display:none;margin-top:10px"></pre>
        <h3>Results</h3>
        <div class="drop" id="dropEda">Drop the Q-Chem EDA output (.out) for the selected model<input type="file" accept=".out,.log,.txt"></div>
        <details style="margin-top:8px"><summary>…or type the terms (kJ/mol)</summary>
          <div class="row" style="margin-top:8px">${['frz', 'elec', 'pauli', 'disp', 'pol', 'ct', 'total'].map(k => `<label class="f">${k.toUpperCase()}<input type="number" step="any" id="em_${k}"></label>`).join('')}</div>
          <div class="row" style="margin-top:6px"><button class="btn sm" id="btnEdaManual">Use these values</button></div></details>
        <div id="eMsg" style="margin-top:8px"></div>
        <details style="margin-top:10px"><summary>Testing aid</summary><p class="muted" style="font-size:12px">Writes a <b>synthetic</b> Q-Chem EDA2 summary for the selected model to try the parser and plots.</p><button class="btn sm" id="btnFakeEda">Download synthetic EDA output</button></details></div>
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Energy decomposition</h2><p class="sub">Attractive terms left of zero, repulsive terms right; ▼ marks the total interaction energy.</p><div id="eChart"><p class="muted">No EDA results yet.</p></div></div>
        <div class="card"><h2>Charge-transfer character</h2><p class="sub">χ<sub>CT</sub> = E<sub>CT</sub> / (E<sub>POL</sub> + E<sub>CT</sub>).</p><div id="eChi"></div></div>
        <div class="card"><h2>Summary</h2><div class="tablewrap" id="eTable"><p class="muted" style="padding:10px">—</p></div></div>
      </div>
    </div>${nextBtn('quality', 'Next: score & report')}`,
    init() {
      $('#btnEdaIn').onclick = () => {
        const id = $('#eModel').value; const m = ST.models[id]; if (!m) return toast('Select a model');
        const verified = ST.cands.filter(c => c.modelId === id).some(c => c.verified);
        if (!verified && !$('#eOverride').checked) { $('#eGateMsg').innerHTML = `<div class="note warn">Model ${id} has no NBO-verified interaction. Run the NBO gate first, or tick the override (it will be recorded).</div>`; return; }
        if (!verified) log('T4.gate-override', `${id}: EDA input generated without NBO verification`);
        const E = ST.P.eda;
        const txt = B.writeQChemEDA({ qm: m.qm, manifest: m.manifest }, { method: E.method, basis: E.basis, eda2: +E.preset, bsse: E.bsse, covp: E.covp, solvent: E.pcm, eps: E.eps, memMB: E.memMB });
        const name = id.replace('BX-', 'bx_') + '_eda.in'; m.files[name] = txt; $('#eCode').style.display = ''; $('#eCode').textContent = txt; download(name, txt);
        $('#eGateMsg').innerHTML = `<div class="note">${verified ? 'Gate passed' : '<b>Gate overridden</b>'} — ${esc(name)} written (also added to the model bundle).</div>`; log('T4.input', `${name}: ${E.method}/${E.basis}, EDA2=${E.preset}, BSSE=${E.bsse}`);
      };
      dropZone($('#dropEda'), async f => { const txt = await readText(f[0]); const r = B.parseEDA(txt); if (r.error) { $('#eMsg').innerHTML = `<div class="note bad">${esc(r.error)}</div>`; return; } this.apply(r, f[0].name, await sha256(txt)); });
      $('#btnFakeEda').onclick = () => { const m = ST.models[$('#eModel').value]; if (!m) return toast('Select a model'); download(m.manifest.id + '_SYNTHETIC_eda.out', syntheticEDA(m.manifest)); };
      $('#btnEdaManual').onclick = () => { const r = { units: 'kJ/mol' }; ['frz', 'elec', 'pauli', 'disp', 'pol', 'ct', 'total'].forEach(k => { const v = $('#em_' + k).value; if (v !== '') r[k] = parseFloat(v); }); if (r.total === undefined) r.total = (r.frz || 0) + (r.pol || 0) + (r.ct || 0); this.apply(r, 'manual entry', '-'); };
    },
    onShow() { const ids = Object.keys(ST.models); $('#eModel').innerHTML = ids.length ? ids.map(id => { const v = ST.cands.some(c => c.modelId === id && c.verified); return `<option value="${id}">${v ? '✔ ' : '✖ '}${id} — ${esc(ST.models[id].manifest.xbs.map(x => x.label).join(' ; ').slice(0, 60))}</option>`; }).join('') : '<option value="">(no models)</option>'; this.render(); },
    apply(r, fname, hash) {
      const id = $('#eModel').value; const m = ST.models[id]; if (!m) return toast('Select the model this output belongs to');
      m.eda = Object.assign({ file: fname, sha256: hash }, r); const d = B.edaDescriptors(r, ST.P);
      ST.cands.filter(c => c.modelId === id).forEach(c => { c.eda = r; if (c.state !== 'NBO-REJECTED' || $('#eOverride').checked) c.state = 'EDA-DONE'; });
      log('T4.eda', `${id}: ${fname} sha256=${hash}; FRZ=${f2(r.frz)} POL=${f2(r.pol)} CT=${f2(r.ct)} TOTAL=${f2(r.total)} kJ/mol; ${d.mechanism}`);
      $('#eMsg').innerHTML = `<div class="note">${esc(fname)} → ${esc(d.mechanism)}${d.character ? '; dominant attraction: ' + esc(d.character) : ''}.${r.terminated === false ? ' <b>No normal-termination marker found.</b>' : ''}</div>`;
      this.render(); setDone('eda');
    },
    render() {
      const done = Object.entries(ST.models).filter(([, m]) => m.eda); if (!done.length) return;
      const rows = done.map(([id, m]) => ({ label: id + ' ' + (m.manifest.xbs[0] ? m.manifest.xbs[0].label.split('···')[0].trim() : ''), v: m.eda, total: m.eda.total }));
      Charts.diverging($('#eChart'), { rows, comps: EDA_COMPS.filter(c => rows.some(r => r.v[c.key] !== undefined)).concat(rows.some(r => r.v.elec === undefined && r.v.frz !== undefined) ? [{ key: 'frz', name: 'Frozen (undecomposed)', color: 'var(--s5)' }] : []), xLabel: 'energy (kJ/mol)' });
      const chis = done.map(([id, m]) => ({ id, d: B.edaDescriptors(m.eda, ST.P) })); const E = ST.P.eda;
      Charts.hbars($('#eChi'), { rows: chis.map(c => ({ label: c.id, value: c.d.chiCT ?? 0, text: f2(c.d.chiCT, 2), color: (c.d.chiCT ?? 0) >= E.chiHigh ? 'var(--s7)' : (c.d.chiCT ?? 0) <= E.chiLow ? 'var(--s2)' : 'var(--s1)', tip: `<b>${c.id}</b><br>χ_CT ${f2(c.d.chiCT, 3)} · ${esc(c.d.mechanism)}` })), max: 1, threshold: E.chiHigh, thresholdLabel: `CT-driven ≥ ${E.chiHigh}`, xLabel: 'χ_CT', fmt: v => v.toFixed(1), labelWidth: 130 });
      $('#eTable').innerHTML = `<table class="t"><thead><tr><th>Model</th><th>FRZ</th><th>ELEC</th><th>PAULI</th><th>DISP</th><th>POL</th><th>CT</th><th>TOTAL</th><th>CT/POL</th><th>Mechanism</th></tr></thead><tbody>${done.map(([id, m]) => { const e = m.eda, d = B.edaDescriptors(e, ST.P); return `<tr><td class="mono">${id}</td>${['frz', 'elec', 'pauli', 'disp', 'pol', 'ct', 'total'].map(k => `<td class="n">${f2(e[k])}</td>`).join('')}<td class="n">${f2(d.ctPolRatio)}</td><td>${esc(d.mechanism)}</td></tr>`; }).join('')}</tbody></table>`;
    }
  };

  /* ====================== 7 QUALITY ====================== */
  R.quality = {
    html: () => `${head('Step 7 · T5', 'Score and report', 'Combines every tier that was run for a candidate into one 0–100 index (a weighted geometric mean, so one weak tier cannot be hidden) and states which evidence it rests on (G·D·E·C).')}
    ${assumeBox('T5')}
    <div class="grid g2">
      <div class="card" style="align-self:start"><h2>Index settings</h2>${pblock('bxqi', 'BXQI settings', stagePaths('T5'), true)}
        <h3>Exports</h3>
        <div class="row"><button class="btn primary" id="btnRep">Report (HTML)</button><button class="btn" id="btnRunJson">Run record (JSON)</button><button class="btn" id="btnQCsv">Scores (CSV)</button></div>
      </div>
      <div class="grid" style="align-content:start">
        <div class="card"><h2>Ranking</h2><div id="qRank"><p class="muted">No candidates yet.</p></div></div>
        <div class="card"><h2>Orbital vs variational charge transfer</h2><p class="sub">NBO ΣE(2) against −E<sub>CT</sub> (ALMO) for candidates with both tiers.</p><div id="qScatter"></div></div>
        <div class="card"><h2>Scorecard</h2><div class="tablewrap" id="qTable"></div></div>
        <div class="card"><h2>Audit trail</h2><div class="tablewrap" id="qAudit" style="max-height:260px"></div></div>
      </div>
    </div>`,
    init() {
      $('#btnRep').onclick = () => { download(`BondX_report_${ST.runId}.html`, buildReport(), 'text/html'); log('report', 'HTML report generated'); };
      $('#btnRunJson').onclick = () => download(`BondX_run_${ST.runId}.json`, JSON.stringify(runRecord(), null, 1), 'application/json');
      $('#btnQCsv').onclick = () => download(`BondX_scores_${ST.runId}.csv`, B.toCSV(ST.cands.map(c => Object.assign({}, c, { s: score(c) })), [{ h: 'id', k: 'cid' }, { h: 'interaction', k: 'label' }, { h: 'halogen', k: 'halogen' }, { h: 'source', k: 'source' }, { h: 'd', k: 'd' }, { h: 'theta1', k: 'theta1' }, { h: 'G', k: 'G' }, { h: 'occupancy', k: 'occupancy' }, { h: 'D', k: 'Dscore' }, { h: 'E2_kcal', k: 'E2' }, { h: 'Eint_kJ', k: c => c.eda && c.eda.total }, { h: 'Epol_kJ', k: c => c.eda && c.eda.pol }, { h: 'Ect_kJ', k: c => c.eda && c.eda.ct }, { h: 'chi_CT', k: c => c.eda && B.edaDescriptors(c.eda, ST.P).chiCT }, { h: 'BXQI', k: c => c.s && c.s.score }, { h: 'grade', k: c => c.s && c.s.grade }, { h: 'provenance', k: c => c.s && c.s.provenance }, { h: 'state', k: 'state' }]));
    },
    onShow() { this.render(); },
    render() {
      const sc = ST.cands.map(c => ({ c, s: score(c) })).filter(x => x.s).sort((a, b) => b.s.score - a.s.score);
      if (sc.length) Charts.hbars($('#qRank'), { rows: sc.map(({ c, s }) => ({ label: `${c.cid} ${c.label.split('···')[1] || ''}`, value: s.score, text: `${s.score.toFixed(0)} · ${s.provenance}`, color: c.state === 'NBO-REJECTED' ? 'var(--neutral)' : HCOLOR[c.halogen], tip: `<b>${c.cid}</b> ${esc(c.label)}<br>BXQI ${s.score.toFixed(1)} — ${esc(s.grade)}<br>tiers: ${Object.entries(s.tiers).map(([k, v]) => k + ' ' + v.toFixed(2)).join(', ')}` })), max: 100, xLabel: 'BXQI (0–100)', fmt: v => v, labelWidth: 230,
        legend: ['I', 'BR', 'CL'].filter(x => sc.some(e => e.c.halogen === x)).map(x => ({ name: HN(x), color: HCOLOR[x] })).concat(sc.some(e => e.c.state === 'NBO-REJECTED') ? [{ name: 'NBO-rejected', color: 'var(--neutral)' }] : []) });
      const both = ST.cands.filter(c => c.E2 != null && c.eda && c.eda.ct != null);
      Charts.scatter($('#qScatter'), { points: both.map(c => ({ x: c.E2, y: -c.eda.ct, group: c.halogen, tip: `<b>${c.cid}</b> ${esc(c.label)}<br>ΣE(2) ${f2(c.E2)} kcal/mol · −E_CT ${f2(-c.eda.ct)} kJ/mol` })), groups: ['I', 'BR', 'CL'].filter(x => both.some(c => c.halogen === x)).map(x => ({ key: x, name: HN(x), color: HCOLOR[x], shape: HSHAPE[x] })), xLabel: 'NBO ΣE(2) (kcal/mol)', yLabel: '−E_CT ALMO (kJ/mol)', height: 260 });
      $('#qTable').innerHTML = ST.cands.length ? `<table class="t"><thead><tr><th>ID</th><th>Interaction</th><th>G</th><th>D</th><th>E</th><th>C</th><th>BXQI</th><th>Grade</th><th>Evidence</th><th>Mechanism</th><th>Consistency</th></tr></thead><tbody>${ST.cands.map(c => { const s = score(c); const d = c.eda ? B.edaDescriptors(c.eda, ST.P) : null; const t = s ? s.tiers : {};
        const tb = k => t[k] !== undefined ? `<div class="bar" data-tip="${k} = ${t[k].toFixed(3)}"><i style="width:${(t[k] * 100).toFixed(0)}%"></i></div>` : '<span class="muted">—</span>';
        const B_ = ST.P.bxqi;
        return `<tr><td class="mono">${c.cid}</td><td>${elTag(c.halogen)} ${esc(c.label)}</td><td>${tb('G')}</td><td>${tb('D')}</td><td>${tb('E')}</td><td>${tb('C')}</td><td class="n"><b>${s ? s.score.toFixed(1) : '—'}</b></td><td>${s ? `<span class="pill ${s.score >= B_.gradeA ? 'ok' : s.score >= B_.gradeB ? 'acc' : s.score >= B_.gradeC ? 'warn' : 'bad'}">${esc(s.grade.split(' · ')[0])}</span>` : ''}</td><td class="mono">${s ? s.provenance : ''}</td><td style="font-size:12px">${d ? esc(d.mechanism) : (c.state === 'NBO-REJECTED' ? '<span class="pill bad">rejected at NBO gate</span>' : '<span class="muted">—</span>')}</td><td style="font-size:11.5px" class="muted">${s && s.consistency ? esc(s.consistency) : ''}</td></tr>`; }).join('')}</tbody></table>` : '<p class="muted" style="padding:10px">No candidates.</p>';
      $('#qAudit').innerHTML = `<table class="t"><thead><tr><th>Time (UTC)</th><th>Action</th><th>Detail</th></tr></thead><tbody>${ST.audit.slice().reverse().map(a => `<tr><td class="mono" style="white-space:nowrap">${a.t.slice(11, 19)}</td><td class="mono">${esc(a.action)}</td><td style="font-size:11.5px">${esc(a.detail)}</td></tr>`).join('')}</tbody></table>`;
      if (sc.length) setDone('quality');
    }
  };

  /* ====================== PARAMETERS ====================== */
  const STAGE_NAMES = { T0: 'Step 1 · Structure input', T1: 'Step 2 · Geometric detection (and Batch / Materials)', T2: 'Step 3 · Dynamics', T3: 'Step 4 · QM model, protonation & charges', T3b: 'Step 5 · NBO gate', T4: 'Step 6 · ALMO-EDA', T5: 'Step 7 · Quality index' };
  R.params = {
    html: () => `${head('Reference', 'All parameters', 'Every number BondX uses, with its default. Changed values are highlighted and recorded in the audit trail, manifests and the report. Save them to reuse the same settings later.')}
    <div class="card"><div class="row" style="justify-content:space-between"><div class="row"><input type="text" class="search" id="pSearch" placeholder="Search parameters…"><span id="pCount" class="muted" style="font-size:12px"></span></div>
      <div class="row"><button class="btn sm" id="pExport">Export parameters (JSON)</button><label class="btn sm" style="margin:0">Import…<input type="file" id="pImport" accept=".json" hidden></label><button class="btn sm" id="pResetAll">Reset all to defaults</button></div></div>
      <div id="pList" style="margin-top:10px"></div></div>`,
    init() {
      $('#pSearch').oninput = () => this.render();
      $('#pExport').onclick = () => download(`BondX_parameters_${ST.runId}.json`, JSON.stringify({ kind: 'BondX-parameters', version: B.VERSION, parameters: ST.P }, null, 1), 'application/json');
      $('#pImport').onchange = async e => { const f = e.target.files[0]; if (!f) return; try { const j = JSON.parse(await readText(f)); if (!j.parameters) throw new Error('not a BondX parameter file'); ST.P = B.mergeParams(B.defaults(), j.parameters); log('params.import', f.name); onParamChanged('all'); toast('Parameters imported'); } catch (er) { toast('Import error: ' + er.message); } e.target.value = ''; };
      $('#pResetAll').onclick = () => { ST.P = B.defaults(); log('params.reset', 'all'); onParamChanged('all'); toast('All parameters reset to defaults'); };
    },
    onShow() { this.render(); },
    render() {
      const q = ($('#pSearch').value || '').toLowerCase(); let total = 0, changed = 0;
      const html = Object.keys(STAGE_NAMES).map(st => { const rows = SPEC.filter(s => s.stage === st && (!q || (s.label + ' ' + s.desc + ' ' + s.path).toLowerCase().includes(q)));
        if (!rows.length) return ''; total += rows.length; changed += nChanged(rows.map(r => r.path));
        return `<h3>${esc(STAGE_NAMES[st])}</h3><div class="ptable pgrid">${rows.map(paramRow).join('')}</div>`; }).join('');
      $('#pList').innerHTML = html || '<p class="muted">No parameter matches.</p>'; $('#pCount').textContent = `${total} parameters shown · ${changed} changed`;
    }
  };

  /* ====================== METHOD ====================== */
  R.method = {
    html: () => { const F = B.FF; return `${head('Reference', 'Method', 'How BondX decides that a halogen contact is a halogen bond, and what each number means.')}
    <div class="card method">
      <div class="flow">
        <div class="box"><b>T1 · Geometry</b><span>σ-hole criteria → G</span></div>
        <div class="box"><b>T2 · Dynamics</b><span>occupancy, residence, best frame → D</span></div>
        <div class="box"><b>T3a · QM model</b><span>force-field protonation & charges, link atoms, manifest</span></div>
        <div class="box gate"><b>T3b · NBO gate</b><span>ΣE(2) n/π→σ*(A–X) ≥ threshold → E</span></div>
        <div class="box"><b>T4 · ALMO-EDA</b><span>FRZ/POL/CT, χ<sub>CT</sub> → C</span></div>
        <div class="box"><b>T5 · BXQI</b><span>weighted geometric fusion</span></div>
      </div>
      <h4>T1 — geometric score</h4>
      <div class="eq">R = d(X···Y) / (r<sub>vdW</sub>(X) + r<sub>vdW</sub>(Y)),  g<sub>d</sub> = clip[(R<sub>cut</sub> − R) / (R<sub>cut</sub> − R<sub>best</sub>)],  g<sub>θ</sub> = clip[(θ₁ − θ<sub>min</sub>) / (180° − θ<sub>min</sub>)]<br>g₂ = exp[−((θ₂ − θ₂<sub>opt</sub>)/w)²] for lone-pair acceptors;  g₂ = 1 − α/α<sub>max</sub> for π acceptors<br>G = [g<sub>d</sub> · g<sub>θ</sub> · (½ + ½ g₂)]<sup>1/3</sup>; G = 0 when a criterion fails.</div>
      <h4>T2 — persistence</h4>
      <div class="eq">f = N<sub>frames with XB</sub> / N<sub>frames</sub>;  τ = mean run length;  D = f · (1 − e<sup>−τ/τ₀</sup>)<sup>p</sup>;  best frame = argmin d, then max θ₁</div>
      <h4>T3a — force-field protonation and charges</h4>
      <p>Residue templates and partial charges are taken from the OpenMM distributions of ${esc(F.ff14SB.label)}, ${esc(F.ff19SB.label)} and ${esc(F.charmm36.label)} (CHARMM patches NTER, GLYP, PROP, CTER, ASPP, GLUP, LSN and DISU applied when the data were built). Waters: ${esc(Object.keys(F.ff14SB.waters).join(', '))} (AMBER), TIP3P (CHARMM). Ions: Joung–Cheatham / Li–Merz (AMBER) and CHARMM ion parameters; other ions take their formal charge. The tables were checked against the independent PDB2PQR AMBER and CHARMM data files (all atoms identical except the CHARMM36 tryptophan update).</p>
      <div class="eq">Hydrogen count per atom = number of H in the residue template (variant chosen from protonation state and termini).<br>Missing H: q(heavy) ← q(heavy) + Σ q(H bonded to it) (residue charge unchanged).<br>QM/MM boundary (M1 bonded to QM): RC — q(M1) → q(M1)/n at each M1–M2 midpoint; Z1/Z2/Z3 — delete q of M1 / M1–M2 / M1–M3.<br>Cut residue: MM atoms shifted equally so that Σq(MM part) = q(residue) − q(QM part).</div>
      <p><b>References:</b> ff14SB — ${esc(F.ff14SB.ref)} · ff19SB — ${esc(F.ff19SB.ref)} · CHARMM36 — ${esc(F.charmm36.ref)} · RC boundary — Lin & Truhlar, J. Phys. Chem. A 2005, 109, 3991.</p>
      <h4>T3b — electronic gate</h4>
      <div class="eq">ΣE(2) = Σ E(2)[n(Y) → σ*(A–X)] (+ Σ E(2)[n(Y) → RY*(X)] if enabled);  aromatic: Σ E(2)[π(ring) → σ*(A–X)];  verified ⇔ ΣE(2) ≥ E<sub>gate</sub></div>
      <h4>T4 — ALMO-EDA descriptors</h4>
      <div class="eq">ΔE<sub>int</sub> = ΔE<sub>FRZ</sub> + ΔE<sub>POL</sub> + ΔE<sub>CT</sub>;  χ<sub>CT</sub> = ΔE<sub>CT</sub> / (ΔE<sub>POL</sub> + ΔE<sub>CT</sub>)</div>
      <h4>T5 — fusion</h4>
      <div class="eq">BXQI = 100 · exp[ Σ<sub>k</sub> (w<sub>k</sub>/Σw) ln s<sub>k</sub> ],  s<sub>G</sub> = G, s<sub>D</sub> = D, s<sub>E</sub> = 1 − e<sup>−ΣE(2)/E₀</sup>, s<sub>C</sub> = 1 − e<sup>−|ΔE<sub>int</sub>|/C₀</sup></div>
      <h4>Key references</h4>
      <p style="font-size:12.5px">Desiraju et al., <i>Pure Appl. Chem.</i> 85, 1711 (2013) · Bondi, <i>J. Phys. Chem.</i> 68, 441 (1964) · Reed, Curtiss & Weinhold, <i>Chem. Rev.</i> 88, 899 (1988) · Horn, Mao & Head-Gordon, <i>Phys. Chem. Chem. Phys.</i> 18, 23067 (2016) · Karmakar & Mishra, <i>J. Chem. Inf. Model.</i> (2026).</p>
    </div>`; },
    init() { }
  };

  /* ------------------------------ synthetic outputs (testing aid only) ------------------------------ */
  function syntheticNBO(mf) {
    const p = (v, n) => String(v).padStart(n); const el = n => { const e = mf.atoms[n - 1].el; return e[0] + e.slice(1).toLowerCase(); };
    let out = ' SYNTHETIC BondX TEST FILE — NOT A REAL CALCULATION\n Entering Gaussian System\n Input orientation:\n ' + '-'.repeat(69) + '\n Center     Atomic      Atomic             Coordinates (Angstroms)\n Number     Number       Type             X           Y           Z\n ' + '-'.repeat(69) + '\n';
    mf.atoms.forEach(a => { out += p(a.n, 7) + p(B.Z[a.el], 11) + p(0, 12) + p(a.x.toFixed(6), 16) + p(a.y.toFixed(6), 12) + p(a.z.toFixed(6), 12) + '\n'; });
    out += ' ' + '-'.repeat(69) + '\n Summary of Natural Population Analysis:\n\n    Atom  No    Charge         Core      Valence    Rydberg      Total\n ' + '-'.repeat(71) + '\n';
    mf.atoms.forEach(a => { const q = ['I', 'BR', 'CL'].includes(a.el) ? 0.12 : a.el === 'O' ? -0.6 : a.el === 'N' ? -0.55 : a.el === 'H' ? 0.21 : -0.08; out += `    ${el(a.n).padEnd(2)}${p(a.n, 4)}${p(q.toFixed(5), 12)}     10.00000     5.00000    0.01000    15.00000\n`; });
    out += ' ' + '='.repeat(71) + '\n\n Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis\n\n     Threshold for printing:   0.50 kcal/mol\n ' + '='.repeat(79) + '\n from unit  1 to unit  2\n';
    let k = 10; mf.xbs.forEach((x, i) => { const e = (x.halogen === 'I' ? 3.1 : x.halogen === 'BR' ? 1.6 : 0.7) * (1 - 0.25 * i);
      if (x.Y) { out += `${p(k++, 5)}. LP (   1) ${el(x.Y)}${p(x.Y, 4)}                / ${p(200 + i, 3)}. BD*(   1) ${el(x.A)}${p(x.A, 4)} - ${el(x.X)}${p(x.X, 4)}         ${p((e * 0.35).toFixed(2), 7)}    0.98    0.031\n`;
        out += `${p(k++, 5)}. LP (   2) ${el(x.Y)}${p(x.Y, 4)}                / ${p(200 + i, 3)}. BD*(   1) ${el(x.A)}${p(x.A, 4)} - ${el(x.X)}${p(x.X, 4)}         ${p((e * 0.65).toFixed(2), 7)}    0.62    0.035\n`; }
      else { const r = x.ring; out += `${p(k++, 5)}. BD (   2) ${el(r[0])}${p(r[0], 4)} - ${el(r[1])}${p(r[1], 4)}       / ${p(200 + i, 3)}. BD*(   1) ${el(x.A)}${p(x.A, 4)} - ${el(x.X)}${p(x.X, 4)}         ${p((e * 0.6).toFixed(2), 7)}    0.55    0.022\n`;
        out += `${p(k++, 5)}. BD (   2) ${el(r[2])}${p(r[2], 4)} - ${el(r[3])}${p(r[3], 4)}       / ${p(200 + i, 3)}. BD*(   1) ${el(x.A)}${p(x.A, 4)} - ${el(x.X)}${p(x.X, 4)}         ${p((e * 0.3).toFixed(2), 7)}    0.55    0.016\n`; } });
    return out + '\n NATURAL BOND ORBITALS (Summary):\n Normal termination of Gaussian 16 (SYNTHETIC)\n';
  }
  function syntheticEDA(mf) {
    const h = mf.xbs[0] ? mf.xbs[0].halogen : 'I'; const s = h === 'I' ? 1 : h === 'BR' ? 0.7 : 0.45;
    const pauli = 38.4 * s, elec = -35.1 * s, disp = -11.7 * s, pol = -4.2 * s, ct = -8.3 * s * s, frz = pauli + elec + disp, tot = frz + pol + ct;
    return `SYNTHETIC BondX TEST FILE — NOT A REAL CALCULATION\n================================\n        Results of EDA2\n================================\nSimplified EDA Summary (kJ/mol)\n--------------------\n PREPARATION      0.0000\n FROZEN         ${frz.toFixed(4)}\n [PAULI         ${pauli.toFixed(4)}]\n [ELEC         ${elec.toFixed(4)}]\n [DISP         ${disp.toFixed(4)}]\n POLARIZATION   ${pol.toFixed(4)}\n CHARGE TRANSFER ${ct.toFixed(4)}\n TOTAL          ${tot.toFixed(4)}     (${tot.toFixed(4)})\n--------------------\n Thank you very much for using Q-Chem. (SYNTHETIC)\n`;
  }

  /* ------------------------------ session, run record, report ------------------------------ */
  function changedParams() { return SPEC.map(s => s.path).filter(p => !same(B.getPath(ST.P, p), B.getPath(DEF, p))).map(p => ({ path: p, value: B.getPath(ST.P, p), default: B.getPath(DEF, p) }))
    .concat(Object.keys(ST.P.ff.overrides || {}).length ? [{ path: 'ff.overrides', value: ST.P.ff.overrides, default: {} }] : []); }
  function runRecord() {
    return { tool: 'BondX Studio', version: B.VERSION, runId: ST.runId, exported: new Date().toISOString(), structure: ST.name ? { name: ST.name, sha256: ST.srcHash, atoms: ST.S && ST.S.atoms.length } : null,
      parameters: ST.P, changedParameters: changedParams(), assumptions: Object.fromEntries(Object.keys(STAGE_NAMES).map(st => [st, B.assumptions(st, ST.P)])), mol2: ST.mol2Name || null,
      candidates: ST.cands.map(c => { const s = score(c); return { id: c.cid, label: c.label, source: c.source, frame: c.frame, halogen: c.halogen, acceptorClass: c.accClass, d: c.d, theta1: c.theta1, theta2: c.theta2, G: c.G, occupancy: c.occupancy, D: c.Dscore, model: c.modelId, E2: c.E2, nbo: c.nbo, eda: c.eda, edaDescriptors: c.eda ? B.edaDescriptors(c.eda, ST.P) : null, bxqi: s, state: c.state }; }),
      models: Object.fromEntries(Object.entries(ST.models).map(([id, m]) => [id, { manifest: m.manifest, nbo: m.nbo || null, eda: m.eda || null }])), audit: ST.audit };
  }
  function saveSession() {
    const sess = { kind: 'BondX-session', version: B.VERSION, runId: ST.runId, name: ST.name, srcText: ST.srcText, srcPQR: ST.srcPQR, P: ST.P, mol2: ST.mol2, mol2Name: ST.mol2Name,
      cands: ST.cands.map(c => { const o = Object.assign({}, c); delete o.mat; if (o.structRef === 'T') o.hit = Object.assign({}, o.hit, { series: undefined }); return o; }),
      models: Object.fromEntries(Object.entries(ST.models).map(([id, m]) => [id, { manifest: m.manifest, qm: m.qm, files: m.files, candIds: m.candIds, nbo: m.nbo, eda: m.eda }])), audit: ST.audit };
    log('session.save', `${ST.cands.length} candidates, ${Object.keys(ST.models).length} models`);
    download(`BondX_session_${ST.runId}.json`, JSON.stringify(sess), 'application/json');
  }
  async function loadSession(file) {
    try { const s = JSON.parse(await readText(file)); if (s.kind !== 'BondX-session') throw new Error('not a BondX session file');
      ST.P = B.mergeParams(B.defaults(), s.P || {}); ST.mol2 = s.mol2 || null; ST.mol2Name = s.mol2Name || '';
      if (s.srcText) await loadStructure(s.srcText, s.name, 'session', s.srcPQR);
      Object.assign(ST, { runId: s.runId, cands: s.cands || [], models: s.models || {}, audit: s.audit || [] });
      ST.cands.forEach(c => { if (c.structRef === 'T' && !ST.traj) c.source += ' (trajectory not reloaded)'; });
      $('#runId').textContent = ST.runId; log('session.load', file.name); toast('Session restored: ' + ST.cands.length + ' candidates'); onParamChanged('session'); if (current && R[current].onShow) R[current].onShow();
    } catch (e) { toast('Session error: ' + e.message, 5000); }
  }
  function buildReport() {
    const rr = runRecord(); const styles = document.querySelector('style').textContent;
    const charts = ['#qRank', '#qScatter', '#eChart', '#eChi', '#nChart', '#dScatter', '#tOccChart'].map(s => $(s)).filter(e => e && e.querySelector('svg')).map(e => `<div class="card" style="margin-bottom:12px">${e.closest('.card').querySelector('h2').outerHTML}${e.innerHTML}</div>`).join('');
    const candRows = rr.candidates.map(c => `<tr><td class="mono">${c.id}</td><td>${esc(c.label)}</td><td>${esc(c.source)}</td><td class="n">${f2(c.d, 3)}</td><td class="n">${f2(c.theta1, 1)}</td><td class="n">${f2(c.G)}</td><td class="n">${c.occupancy != null ? f2(c.occupancy * 100, 1) + '%' : '—'}</td><td class="n">${f2(c.E2)}</td><td class="n">${c.eda ? f2(c.eda.pol) : '—'}</td><td class="n">${c.eda ? f2(c.eda.ct) : '—'}</td><td class="n">${c.eda ? f2(c.eda.total) : '—'}</td><td>${c.edaDescriptors ? esc(c.edaDescriptors.mechanism) : ''}</td><td class="n"><b>${c.bxqi ? c.bxqi.score.toFixed(1) : '—'}</b></td><td>${c.bxqi ? esc(c.bxqi.grade) + ' [' + c.bxqi.provenance + ']' : ''}</td><td>${esc(c.state)}</td></tr>`).join('');
    const models = Object.entries(rr.models).map(([id, m]) => `<tr><td class="mono">${id}</td><td>${m.manifest.nAtoms}</td><td>${esc(m.manifest.env || '')}</td><td>${esc(m.manifest.level || '')}</td><td>${esc(m.manifest.forceField || '')}</td><td>${JSON.stringify(m.manifest.fragCharge)} / ${m.manifest.multiplicity}</td><td>${m.nbo ? esc(m.nbo.file) + ' · ' + esc(m.nbo.verification.note) : '—'}</td><td>${m.eda ? esc(m.eda.file) : '—'}</td></tr>`).join('');
    const ch = rr.changedParameters; const asm = Object.entries(rr.assumptions).map(([st, a]) => `<h3>${esc(STAGE_NAMES[st])}</h3><ul>${a.map(x => `<li>${esc(x)}</li>`).join('')}</ul>`).join('');
    return `<!DOCTYPE html><html><head><meta charset="utf-8"><title>BondX report ${rr.runId}</title><style>${styles} body{background:#fff} main{max-width:1150px;margin:0 auto} @media print{.card{break-inside:avoid}}</style></head><body><main>
      <h1 style="margin:0">BondX halogen-bond characterisation report</h1><p class="muted">Run ${rr.runId} · generated ${rr.exported} · BondX Studio${rr.structure ? ` · structure ${esc(rr.structure.name)} (SHA-256 ${esc((rr.structure.sha256 || '').slice(0, 24))}…)` : ''}</p>
      <div class="card" style="margin:12px 0"><h2>Candidate scorecard</h2><table class="t"><thead><tr><th>ID</th><th>Interaction</th><th>Source</th><th>d</th><th>θ₁</th><th>G</th><th>Occ.</th><th>ΣE(2)</th><th>E<sub>POL</sub></th><th>E<sub>CT</sub></th><th>E<sub>int</sub></th><th>Mechanism</th><th>BXQI</th><th>Grade</th><th>State</th></tr></thead><tbody>${candRows}</tbody></table><p class="muted" style="font-size:11.5px">Distances Å, angles °, E(2) kcal/mol, EDA kJ/mol.</p></div>
      ${charts}
      <div class="card" style="margin-bottom:12px"><h2>QM models</h2><table class="t"><thead><tr><th>Model</th><th>Atoms</th><th>Environment</th><th>Level</th><th>Force field</th><th>Charges / mult.</th><th>NBO output</th><th>EDA output</th></tr></thead><tbody>${models}</tbody></table></div>
      <div class="card" style="margin-bottom:12px"><h2>Changed parameters</h2>${ch.length ? `<table class="t"><thead><tr><th>Parameter</th><th>Value</th><th>Default</th></tr></thead><tbody>${ch.map(c => `<tr><td class="mono">${esc(c.path)}</td><td class="mono">${esc(JSON.stringify(c.value))}</td><td class="mono">${esc(JSON.stringify(c.default))}</td></tr>`).join('')}</tbody></table>` : '<p>All parameters at their defaults.</p>'}</div>
      <div class="card" style="margin-bottom:12px"><h2>Assumptions</h2>${asm}</div>
      <div class="card" style="margin-bottom:12px"><h2>All parameters</h2><pre class="code">${esc(JSON.stringify(rr.parameters, null, 1))}</pre></div>
      <div class="card"><h2>Audit trail</h2><table class="t"><thead><tr><th>Time (UTC)</th><th>Action</th><th>Detail</th></tr></thead><tbody>${rr.audit.map(a => `<tr><td class="mono">${a.t}</td><td class="mono">${esc(a.action)}</td><td>${esc(a.detail)}</td></tr>`).join('')}</tbody></table></div>
    </main></body></html>`;
  }

  /* ------------------------------ visitor counter ------------------------------ */
  // Anonymous page-view counter (GoatCounter). Active only when the page is served from the web (https) AND a
  // counter address is set in <meta name="bondx-counter">. It sends the page address, browser type and screen size;
  // never any structure, file or result. Opening the HTML file from disk (file://) never counts.
  function startVisitorCounter() {
    const m = document.querySelector('meta[name="bondx-counter"]'); const url = m && m.content.trim();
    if (!url || location.protocol !== 'https:') return;
    const s = document.createElement('script'); s.async = true; s.src = 'https://gc.zgo.at/count.js'; s.setAttribute('data-goatcounter', url); document.head.appendChild(s);
  }

  /* ------------------------------ boot ------------------------------ */
  function boot() {
    $('#runId').textContent = ST.runId; buildSide();
    $('#main').innerHTML = PAGES.map(p => `<section class="page" id="pg-${p.id}">${R[p.id].html()}</section>`).join('');
    PAGES.forEach(p => R[p.id].init());
    $('#btnSaveSession').onclick = saveSession; $('#inSession').onchange = e => { if (e.target.files[0]) loadSession(e.target.files[0]); e.target.value = ''; };
    $('#btnReport').onclick = () => { download(`BondX_report_${ST.runId}.html`, buildReport(), 'text/html'); log('report', 'HTML report generated'); };
    $('#btnParams').onclick = () => go('params');
    $('#btnTutorial').onclick = () => go('guide');
    log('session.start', `BondX Studio; ${navigator.userAgent}`);
    startVisitorCounter();
    go('structure');
    V.load().then(ok => { if (ok) { V.viewer = null; const slot = V.host && V.host.parentNode; if (V.host) V.host.remove(); V.host = null; V.ensure(); if (slot) V.mount(slot); if (ST.S) V.structure(ST.S); } else V.offline(); });
  }
  boot();
})();
