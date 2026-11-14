/* AeroMargin demo core + UI. Team Cheese Pizza, NASA Space Apps 2026.
   Terrain is synthetic (procedural), shaped after Melas Chasma landforms.
   The oxygen model, thresholds and PNR rule follow the team's technical specification. */
var AeroCore = (function () {
  'use strict';
  var W = 256, H = 160, CELL = 40, N = W * H;
  var MG = 195 * 3.721; // suited mass 195 kg in Mars gravity = 725.6 N
  var P = {
    Wbase: 115, c1: 0.28, cup: 1.95, cdown: 0.45, Wsink: 38, Ibed: 452, exp: 1.65,
    k: 0.00297, leak: 0.035, psiPerL: 4.38, reserve: 45, kappa: 0.55,
    doseRef: 26.7,            // microSv/h, = 0.64 mSv/day (Curiosity RAD, Gale crater)
    vNom: 0.8, vShed: 0.65,   // walking speed m/s (0.8 is a demo assumption; 0.65 from the spec)
    brake: 1.2,               // braking penalty below -12 deg (demo assumption)
    stride: 0.6,              // suited step length in metres (demo assumption)
    sigmaNoise: 0.25, lambdaRad: 0.02
  };

  function hash(ix, iy, s) {
    var h = (Math.imul(ix, 374761393) + Math.imul(iy, 668265263) + Math.imul(s, 1442695041)) | 0;
    h = Math.imul(h ^ (h >>> 13), 1274126177);
    h ^= h >>> 16;
    return (h >>> 0) / 4294967296;
  }
  function vnoise(x, y, s) {
    var ix = Math.floor(x), iy = Math.floor(y), fx = x - ix, fy = y - iy;
    var sx = fx * fx * (3 - 2 * fx), sy = fy * fy * (3 - 2 * fy);
    var a = hash(ix, iy, s), b = hash(ix + 1, iy, s), c = hash(ix, iy + 1, s), d = hash(ix + 1, iy + 1, s);
    return a + (b - a) * sx + (c - a) * sy + (a - b - c + d) * sx * sy;
  }
  function fbm(x, y, s, oct) {
    var v = 0, a = 0.5, f = 1, t = 0;
    for (var i = 0; i < oct; i++) { v += a * vnoise(x * f, y * f, s + i * 17); t += a; a *= 0.5; f *= 2; }
    return v / t;
  }
  function sstep(e0, e1, x) { var t = Math.min(1, Math.max(0, (x - e0) / (e1 - e0))); return t * t * (3 - 2 * t); }
  function mix(a, b, t) { return a + (b - a) * t; }

  var HAB = { x: 40, y: 112 };
  var MOUND = { x: 150, y: 92, s: 16, h: 260 };

  function rawHeight(x, y, out) {
    var xm = x * CELL, ym = y * CELL, v = y / H;
    var h = (fbm(xm / 2200, ym / 2200, 11, 3) - 0.5) * 150 + (fbm(xm / 520, ym / 520, 29, 2) - 0.5) * 26;
    // north canyon wall with a wavy foot and a talus apron
    var edge = 0.27 + 0.12 * (fbm(xm / 1800, 3.3, 5, 2) - 0.5);
    var wall = sstep(edge, edge - 0.24, v);
    wall *= 1 + 0.15 * (vnoise(xm / 300, ym / 900, 7) - 0.5);
    var talus = sstep(edge + 0.14, edge, v);
    h += 1900 * wall + 180 * talus;
    // layered sulfate mound
    var dx = x - MOUND.x, dy = y - MOUND.y;
    var g = Math.exp(-(dx * dx + dy * dy) / (2 * MOUND.s * MOUND.s));
    var hm = MOUND.h * g, u = hm / 40, q = (Math.floor(u) + sstep(0.35, 0.65, u - Math.floor(u))) * 40;
    h += mix(hm, q, 0.3);
    // block field in the south-east
    var bx = x - 214, by = y - 126, bf = Math.exp(-(bx * bx + by * by) / (2 * 20 * 20));
    var blk = sstep(0.56, 0.62, vnoise(xm / 260, ym / 260, 41)) * sstep(0.15, 0.4, bf);
    h += 45 * blk;
    // two graben (fault troughs) in the south
    var gr = 0;
    for (var i = 0; i < 2; i++) {
      var ox = 96 + i * 4, oy = 148 - i * 15, ux = 0.944, uy = -0.33;
      var rx = x - ox, ry = y - oy;
      var along = rx * ux + ry * uy, across = Math.abs(rx * -uy + ry * ux) * CELL;
      gr = Math.max(gr, sstep(300, 100, across) * sstep(-14, -4, along) * sstep(52, 40, along));
    }
    h -= 50 * gr;
    if (out) { out.wall = wall; out.talus = talus; out.g = g; out.blk = blk; out.gr = gr; }
    return h;
  }

  function build() {
    var z = new Float32Array(N), I = new Float32Array(N), slope = new Float32Array(N), psi = new Float32Array(N);
    var block = new Uint8Array(N), hazard = new Uint8Array(N), rsl = new Uint8Array(N), dune = new Float32Array(N), land = new Uint8Array(N);
    var o = {}, x, y, i;
    var zHab = rawHeight(HAB.x, HAB.y, null);
    for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
      i = y * W + x;
      var h = rawHeight(x, y, o), xm = x * CELL, ym = y * CELL;
      var dh = Math.hypot(x - HAB.x, y - HAB.y), pad = sstep(10, 4, dh);
      h = mix(h, zHab, pad);
      z[i] = h - zHab;
      var ti = 295 + (fbm(xm / 900, ym / 900, 53, 2) - 0.5) * 120;
      var dn = sstep(0.53, 0.61, fbm(xm / 1500, ym / 1100, 77, 2)) * (1 - sstep(0.02, 0.2, o.wall)) * (1 - sstep(0.1, 0.35, o.g)) * (1 - o.blk);
      ti = mix(ti, 160, dn);
      ti = mix(ti, 452, sstep(0.6, 0.68, fbm(xm / 1300 + 9, ym / 1300, 91, 2)) * (1 - dn));
      ti = mix(ti, 240, o.talus * (1 - sstep(0.1, 0.4, o.wall)) * 0.8);
      ti = mix(ti, 430, sstep(0.25, 0.6, o.wall));
      ti = mix(ti, 410, sstep(0.15, 0.5, o.g));
      ti = mix(ti, 452, o.blk);
      ti = mix(ti, 400, pad);
      I[i] = ti; dune[i] = dn * (1 - pad);
      // landform class for the terrain photo preview: 0 plain, 1 wall, 2 talus, 3 mound, 4 blocks, 5 graben, 6 dunes
      land[i] = o.wall > 0.25 ? 1 : o.blk > 0.5 ? 4 : o.gr > 0.4 ? 5 : o.g > 0.2 ? 3 : dn > 0.5 ? 6 : o.talus > 0.5 ? 2 : 0;
    }
    for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
      i = y * W + x;
      var xl = Math.max(0, x - 1), xr = Math.min(W - 1, x + 1), yu = Math.max(0, y - 1), yd = Math.min(H - 1, y + 1);
      var gx = (z[y * W + xr] - z[y * W + xl]) / ((xr - xl) * CELL), gy = (z[yd * W + x] - z[yu * W + x]) / ((yd - yu) * CELL);
      slope[i] = Math.atan(Math.hypot(gx, gy)) * 180 / Math.PI;
    }
    // recurring slope lineae sites: steep wall slopes in a few bands
    var tmp = new Uint8Array(N);
    for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
      i = y * W + x;
      if (slope[i] > 25 && y > H * 0.06 && y < H * 0.36 && fbm(x * CELL / 620, y * CELL / 1100, 63, 2) > 0.64) tmp[i] = 1;
    }
    for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
      i = y * W + x; var r = tmp[i];
      if (!r) for (var k = 0; k < 4 && !r; k++) {
        var nx = x + [1, -1, 0, 0][k], ny = y + [0, 0, 1, -1][k];
        if (nx >= 0 && ny >= 0 && nx < W && ny < H && tmp[ny * W + nx]) r = 1;
      }
      rsl[i] = r;
    }
    // sky-view factor from terrain horizons (8 directions, 800 m reach), clipped to [0.45, 1]
    var D8 = [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1]];
    for (y = 0; y < H; y++) for (x = 0; x < W; x++) {
      i = y * W + x; var acc = 0;
      for (var d = 0; d < 8; d++) {
        var mt = 0, ddx = D8[d][0], ddy = D8[d][1], sl = Math.hypot(ddx, ddy) * CELL;
        for (var st = 2; st <= 20; st += 2) {
          var px = x + ddx * st, py = y + ddy * st;
          if (px < 0 || py < 0 || px >= W || py >= H) break;
          var t = (z[py * W + px] - z[i]) / (st * sl);
          if (t > mt) mt = t;
        }
        acc += mt / Math.sqrt(1 + mt * mt);
      }
      psi[i] = Math.min(1, Math.max(0.45, 1 - 1.6 * acc / 8));
      block[i] = (slope[i] > 30 || rsl[i]) ? 1 : 0;
      hazard[i] = (I[i] < 200 || slope[i] > 18) ? 1 : 0;
    }
    return { z: z, I: I, slope: slope, psi: psi, block: block, hazard: hazard, rsl: rsl, dune: dune, land: land, block0: new Uint8Array(block), hab: HAB.y * W + HAB.x };
  }

  function metW(T, a, b, d, v) {
    var dz = T.z[b] - T.z[a], s = dz / Math.hypot(d, dz);
    var Wm = P.Wbase + P.c1 * MG * v + P.cup * MG * v * Math.max(0, s) + P.Wsink * Math.pow(P.Ibed / T.I[b], P.exp);
    if (s < 0) {
      var sd = -s;
      Wm -= P.cdown * MG * v * Math.min(sd, 0.1736);
      if (sd > 0.2079) Wm += P.brake * MG * v * (sd - 0.2079);
    }
    return Math.max(P.Wbase, Wm);
  }
  function flow(Wm) { return P.k * Wm + P.leak; } // L/min STP
  function edge(T, a, b, d, v) {
    var Wm = metW(T, a, b, d, v), len = Math.hypot(d, T.z[b] - T.z[a]), t = len / v;
    return { o2: flow(Wm) * t / 60, t: t, W: Wm, len: len };
  }
  var DX = [1, 1, 0, -1, -1, -1, 0, 1], DY = [0, 1, 1, 1, 0, -1, -1, -1];

  // pen: optional per-cell cost multiplier (safest-route search). T.uz: optional mask of cells closed by crew
  // hazard reports; they can be left but never entered, so an astronaut standing inside one can still walk out.
  function dijkstra(T, src, v, reverse, pen) {
    var uz = T.uz || null;
    var dist = new Float64Array(N), prev = new Int32Array(N), i;
    for (i = 0; i < N; i++) { dist[i] = Infinity; prev[i] = -1; }
    var hk = [], hv = [];
    function push(k, val) {
      var n = hk.length; hk.push(k); hv.push(val);
      while (n > 0) { var p = (n - 1) >> 1; if (hk[p] <= k) break; hk[n] = hk[p]; hv[n] = hv[p]; n = p; }
      hk[n] = k; hv[n] = val;
    }
    function pop() {
      var k0 = hk[0], v0 = hv[0], k = hk.pop(), val = hv.pop(), n = hk.length;
      if (n > 0) {
        var p = 0;
        for (;;) {
          var c = 2 * p + 1; if (c >= n) break;
          if (c + 1 < n && hk[c + 1] < hk[c]) c++;
          if (hk[c] >= k) break;
          hk[p] = hk[c]; hv[p] = hv[c]; p = c;
        }
        hk[p] = k; hv[p] = val;
      }
      return [k0, v0];
    }
    dist[src] = 0; push(0, src);
    while (hk.length) {
      var top = pop(), d = top[0], u = top[1];
      if (d > dist[u]) continue;
      var ux = u % W, uy = (u / W) | 0;
      for (var k = 0; k < 8; k++) {
        var nx = ux + DX[k], ny = uy + DY[k];
        if (nx < 0 || ny < 0 || nx >= W || ny >= H) continue;
        var n = ny * W + nx;
        if (T.block[n] && !(uz && uz[n])) continue;
        if (uz && (reverse ? (uz[u] && !uz[n]) : (uz[n] && !uz[u]))) continue;
        var dd = (k & 1) ? CELL * 1.41421356 : CELL;
        var c = reverse ? edge(T, n, u, dd, v).o2 : edge(T, u, n, dd, v).o2;
        if (pen) c *= pen[n];
        var nd = d + c;
        if (nd < dist[n]) { dist[n] = nd; prev[n] = u; push(nd, n); }
      }
    }
    return { dist: dist, prev: prev };
  }
  // chain of cells following prev pointers from `from` until `to` (reverse field: prev = next hop to hab)
  function chain(prev, from, to) {
    var out = [from], c = from, guard = 0;
    while (c !== to && prev[c] >= 0 && guard++ < N) { c = prev[c]; out.push(c); }
    return out;
  }
  function beta(F) { return 1 / (1 + P.kappa * F); }

  return { W: W, H: H, CELL: CELL, N: N, P: P, MG: MG, HAB: HAB, MOUND: MOUND, build: build, metW: metW, flow: flow,
    edge: edge, dijkstra: dijkstra, chain: chain, beta: beta, sstep: sstep };
})();
if (typeof module !== 'undefined') module.exports = AeroCore;

if (typeof document !== 'undefined') (function () {
  'use strict';
  var C = AeroCore, W = C.W, H = C.H, CELL = C.CELL, N = C.N, P = C.P, S = 5;
  var $ = function (id) { return document.getElementById(id); };
  var cv = $('map'); if (!cv) return;
  cv.width = W * S; cv.height = H * S;
  var ctx = cv.getContext('2d');
  var T = C.build();
  T.cls = new Uint8Array(N);      // safety class per cell: 0 safe (green), 1 caution (yellow), 2 no entry (red)
  T.pen = new Float32Array(N);    // cost multiplier used by the safest-route search
  T.zoneAt = new Int16Array(N);   // index of the crew report covering the cell, or -1
  T.uz = new Uint8Array(N);       // cells closed by a crew report (not by terrain)
  var baseCv = document.createElement('canvas'); baseCv.width = W; baseCv.height = H;
  var overCv = document.createElement('canvas'); overCv.width = W; overCv.height = H;
  var bctx = baseCv.getContext('2d'), octx = overCv.getContext('2d');
  var reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var PEN_TERRAIN = 4;            // safest route: caution terrain costs 4x

  function findSpot(x, y) { // nearest walkable, gentle cell
    for (var r = 0; r < 30; r++) for (var dy = -r; dy <= r; dy++) for (var dx = -r; dx <= r; dx++) {
      if (Math.max(Math.abs(dx), Math.abs(dy)) !== r) continue;
      var px = x + dx, py = y + dy; if (px < 1 || py < 1 || px >= W - 1 || py >= H - 1) continue;
      var i = py * W + px; if (!T.block[i] && T.slope[i] < 12 && T.I[i] > 220) return i;
    }
    return y * W + x;
  }
  // Science targets are numbered S1 to S5 so that the letters A and B stay free for the route start and destination.
  var SCI0 = [
    { id: 'S1', name: 'Layered sulfates', ru: 'слоистые сульфаты', note: 'Hydrated sulfate layers, a sign of past water', cls: 'sulfate', s2: 4.0, x: 131, y: 96 },
    { id: 'S2', name: 'Layer contact', ru: 'контакт слоёв', note: 'Where layered deposits meet older rock', cls: 'sulfate', s2: 3.0, x: 171, y: 84 },
    { id: 'S3', name: 'Block field', ru: 'поле блоков', note: 'Bright 100–500 m blocks in dark material', cls: 'bedrock', s2: 2.0, x: 196, y: 122 },
    { id: 'S4', name: 'RSL standoff', ru: 'точка наблюдения RSL', note: 'Viewpoint below seasonal slope streaks', cls: 'rsl', s2: 5.0, x: 108, y: 62 },
    { id: 'S5', name: 'Graben scarp', ru: 'уступ грабена', note: 'Fault trough cutting the canyon floor', cls: 'tectonic', s2: 1.5, x: 104, y: 136 }
  ];
  var SCI = [];
  var st = {};
  var ret = null;
  var out = null, outKey = '';
  var safe = null, safeKey = '';
  var userLayer = 'safety';
  var mode = 'safe';              // 'safe' = safest route, 'fast' = lowest oxygen cost
  var clickMode = 'B';            // what a map tap sets: 'A', 'B' or 'haz'
  var vlang = /^(ru|kk|uk|be)/i.test(navigator.language || '') ? 'ru' : 'en';

  /* ---------- crew hazard reports: storage ---------- */
  var LS_KEY = 'aeromargin.hazards.v1';
  var DB_URL = String(window.AEROMARGIN_DB_URL || '').replace(/\/+$/, '');
  var HZ_TYPES = {
    crevasse: ['Crevasse', 'трещина'], pit: ['Sinkhole or pit', 'провал'], crater: ['Crater', 'кратер'],
    sand: ['Loose sand', 'рыхлый песок'], rocks: ['Boulders or rockfall', 'камни, обвал'],
    slope: ['Steep slope', 'крутой склон'], other: ['Other hazard', 'другая опасность']
  };
  var LEVELS = ['', 'Level 1 · caution', 'Level 2 · avoid', 'Level 3 · no entry'];
  var SEEDS = [
    { id: 'seed-1', x: 74, y: 124, r: 4, type: 'crevasse', level: 3, by: 'EV-2', t: 0, photo: '', seed: true, desc: 'Sample report: open crack about 2 m wide, across the direct line to the graben.' },
    { id: 'seed-2', x: 90, y: 104, r: 5, type: 'rocks', level: 2, by: 'EV-3', t: 0, photo: '', seed: true, desc: 'Sample report: loose boulders, unstable footing.' }
  ];
  var zones = [], zver = 0, dbState = DB_URL ? 'connecting' : 'local', pendingZones = null;
  function clamp(v, a, b) { return Math.min(b, Math.max(a, v)); }
  function uid() { return 'hz-' + Date.now().toString(36) + '-' + Math.floor(Math.random() * 1e6).toString(36); }
  function cleanZone(z) { // everything read from storage or the network goes through here
    if (!z || typeof z !== 'object') return null;
    var x = Math.round(+z.x), y = Math.round(+z.y);
    if (!(x >= 0 && x < W && y >= 0 && y < H)) return null;
    var ph = (typeof z.photo === 'string' && z.photo.length < 200000 && /^data:image\/(jpeg|png|webp);base64,[A-Za-z0-9+\/=]+$/.test(z.photo)) ? z.photo : '';
    return {
      id: String(z.id || '').replace(/[^\w-]/g, '').slice(0, 40) || uid(), x: x, y: y,
      r: clamp(Math.round(+z.r) || 3, 1, 12), type: HZ_TYPES.hasOwnProperty(z.type) ? z.type : 'other',
      level: clamp(Math.round(+z.level) || 2, 1, 3), desc: String(z.desc || '').slice(0, 300),
      by: String(z.by || 'EV-?').slice(0, 20), t: +z.t || 0, photo: ph, seed: !!z.seed, local: !!z.local
    };
  }
  function loadLocal() {
    try { var raw = localStorage.getItem(LS_KEY); if (raw) { var a = JSON.parse(raw); if (Array.isArray(a)) return a.map(cleanZone).filter(Boolean); } } catch (e) {}
    return null;
  }
  function saveLocal() { try { localStorage.setItem(LS_KEY, JSON.stringify(zones)); } catch (e) {} }
  function remote(method, path, body) {
    return fetch(DB_URL + path + '.json', { method: method, headers: body ? { 'Content-Type': 'application/json' } : undefined, body: body ? JSON.stringify(body) : undefined })
      .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); });
  }
  function wire(z) { var o = {}; for (var k in z) if (k !== 'local' && k !== 'seed') o[k] = z[k]; return o; }
  function pushZone(z) {
    if (!DB_URL) return;
    remote('PUT', '/hazards/' + z.id, wire(z)).then(function () { z.local = false; dbState = 'online'; saveLocal(); zoneList(); })
      .catch(function () { dbState = 'offline'; zoneList(); });
  }
  function pull() {
    if (!DB_URL) return;
    remote('GET', '/hazards').then(function (obj) {
      var list = [], seen = {};
      if (obj && typeof obj === 'object') Object.keys(obj).forEach(function (k) { var z = cleanZone(obj[k]); if (z) { z.local = false; z.seed = false; list.push(z); seen[z.id] = 1; } });
      zones.forEach(function (z) { if (z.local && !seen[z.id]) { list.push(z); pushZone(z); } }); // reports made while offline
      list.sort(function (a, b) { return a.t - b.t; });
      dbState = 'online';
      if (JSON.stringify(list) !== JSON.stringify(zones)) setZones(list); else zoneList();
    }).catch(function () { dbState = 'offline'; zoneList(); });
  }
  function setZones(list, noSave) {
    if (st.walking) { pendingZones = list; return; }
    zones = list; if (!noSave) saveLocal();
    applyZones(); update();
  }
  // Rebuild the safety classes, the route penalties and the closed cells from terrain plus crew reports.
  function applyZones() {
    var i, hx = C.HAB.x, hy = C.HAB.y;
    for (i = 0; i < N; i++) {
      var b = T.block0[i]; T.block[i] = b; T.uz[i] = 0; T.zoneAt[i] = -1;
      T.cls[i] = b ? 2 : T.hazard[i] ? 1 : 0; T.pen[i] = T.hazard[i] ? PEN_TERRAIN : 1;
    }
    zones.forEach(function (z, zi) {
      var R2 = z.r + (z.level === 3 ? 2 : 0); // a red zone gets a two-cell yellow buffer
      for (var dy = -R2; dy <= R2; dy++) for (var dx = -R2; dx <= R2; dx++) {
        var px = z.x + dx, py = z.y + dy; if (px < 0 || py < 0 || px >= W || py >= H) continue;
        var d = Math.hypot(dx, dy); if (d > R2 + 0.01) continue;
        if (Math.hypot(px - hx, py - hy) < 5) continue; // the airlock pad always stays open
        var q = py * W + px, cur = T.zoneAt[q];
        if (d <= z.r + 0.01) {
          if (z.level === 3) { if (!T.block0[q]) T.uz[q] = 1; T.block[q] = 1; T.cls[q] = 2; }
          else { if (T.cls[q] < 1) T.cls[q] = 1; T.pen[q] = Math.max(T.pen[q], z.level === 2 ? 8 : 3); }
          if (cur < 0 || zones[cur].level <= z.level) T.zoneAt[q] = zi;
        } else {
          if (T.cls[q] < 1) T.cls[q] = 1; T.pen[q] = Math.max(T.pen[q], 3);
          if (cur < 0) T.zoneAt[q] = zi;
        }
      }
    });
    zver++;
    ret = { nom: C.dijkstra(T, T.hab, P.vNom, true), shed: C.dijkstra(T, T.hab, P.vShed, true) };
    outKey = ''; safeKey = '';
  }

  function reset() {
    SCI = SCI0.map(function (s) { return { id: s.id, name: s.name, ru: s.ru, note: s.note, cls: s.cls, s2: s.s2, cell: findSpot(s.x, s.y), done: false }; });
    st = { V: 1500 / P.psiPerL, F: 0.2, pos: T.hab, target: -1, sci: null, elapsed: 0, dose: 0, walked: 0, walking: false, msg: '', walkCells: null, walkIdx: 0, ping: -1, pingT: 0 };
    st.ema = C.flow(C.metW(T, T.hab, T.hab + 1, CELL, P.vNom));
    st.W = 0;
    $('tank').value = 1500; $('fat').value = 20;
    outKey = ''; safeKey = '';
  }
  function shedMode() { return st.F > 0.55; }
  function speed() { return shedMode() ? P.vShed : P.vNom; }
  function R() { return shedMode() ? ret.shed : ret.nom; }
  function ensureOut() {
    var key = st.pos + ':' + (shedMode() ? 1 : 0) + ':' + zver;
    if (key !== outKey) { out = C.dijkstra(T, st.pos, speed(), false); outKey = key; }
  }
  function ensureSafe() {
    var key = st.pos + ':' + (shedMode() ? 1 : 0) + ':' + zver;
    if (key !== safeKey) { safe = C.dijkstra(T, st.pos, speed(), false, T.pen); safeKey = key; }
  }

  /* ---------- layers ---------- */
  function ramp(stops, t) {
    t = Math.min(1, Math.max(0, t)); var n = stops.length - 1, f = t * n, i = Math.min(n - 1, Math.floor(f)), u = f - i;
    return [stops[i][0] + (stops[i + 1][0] - stops[i][0]) * u, stops[i][1] + (stops[i + 1][1] - stops[i][1]) * u, stops[i][2] + (stops[i + 1][2] - stops[i][2]) * u];
  }
  var RAMPS = {
    relief: [[74, 44, 33], [150, 94, 62], [204, 152, 108], [233, 201, 164]],
    thermal: [[52, 104, 168], [120, 140, 150], [205, 150, 96], [244, 150, 60]],
    rad: [[48, 30, 84], [140, 62, 130], [232, 120, 110], [250, 214, 170]],
    o2: [[34, 110, 104], [112, 176, 120], [236, 204, 98], [226, 110, 62], [150, 40, 44]]
  };
  var SAFE_RGB = [[70, 150, 96], [222, 196, 84], [200, 52, 44]];
  var LEG = {
    safety: ['Safety zones: green safe, yellow caution, red no entry', 'safe', 'no entry', [SAFE_RGB[0], SAFE_RGB[0], SAFE_RGB[1], SAFE_RGB[2], SAFE_RGB[2]]],
    relief: ['Synthetic shaded relief', 'low', 'high', RAMPS.relief],
    slope: ['Slope, degrees', '0°', '25°+', [[70, 150, 96], [222, 196, 84], [226, 128, 52], [200, 52, 44]]],
    thermal: ['Thermal inertia, TIU', '150 loose sand', '452 bedrock', RAMPS.thermal],
    rad: ['Dose rate, µSv/h', '12 sheltered', '26.7 open sky', RAMPS.rad],
    o2: ['O₂ to walk home, litres', '0', '200+', RAMPS.o2]
  };
  function reliefColor(i) {
    var zz = T.z[i], c = ramp(RAMPS.relief, zz < 320 ? (zz + 90) / 410 * 0.66 : 0.66 + Math.min(1, (zz - 320) / 1700) * 0.34);
    var dn = T.dune[i]; return [c[0] + (62 - c[0]) * dn * 0.7, c[1] + (52 - c[1]) * dn * 0.7, c[2] + (56 - c[2]) * dn * 0.7];
  }
  function paintBase(layer) {
    var img = bctx.createImageData(W, H), d = img.data, Rr = R().dist;
    for (var y = 0; y < H; y++) for (var x = 0; x < W; x++) {
      var i = y * W + x, xl = Math.max(0, x - 1), xr = Math.min(W - 1, x + 1), yu = Math.max(0, y - 1), yd = Math.min(H - 1, y + 1);
      var gx = (T.z[y * W + xr] - T.z[y * W + xl]) / ((xr - xl) * CELL), gy = (T.z[yd * W + x] - T.z[yu * W + x]) / ((yd - yu) * CELL);
      var sh = (gx * 0.62 - gy * 0.62 + 0.48) / Math.sqrt(gx * gx + gy * gy + 1) / 0.48; // light from south-west, so the south-facing wall is lit
      sh = Math.min(1.3, Math.max(0.5, sh));
      var c, k = 1;
      if (layer === 'relief') { c = reliefColor(i); k = sh; }
      else if (layer === 'safety') {
        var rc = reliefColor(i), sc = SAFE_RGB[T.cls[i]], m = T.cls[i] ? 0.62 : 0.42;
        c = [rc[0] + (sc[0] - rc[0]) * m, rc[1] + (sc[1] - rc[1]) * m, rc[2] + (sc[2] - rc[2]) * m]; k = 0.35 + 0.65 * sh;
      }
      else if (layer === 'slope') { c = ramp(LEG.slope[3], T.slope[i] / 25); k = 0.55 + 0.45 * sh; }
      else if (layer === 'thermal') { c = ramp(RAMPS.thermal, (T.I[i] - 150) / 302); k = 0.55 + 0.45 * sh; }
      else if (layer === 'rad') { c = ramp(RAMPS.rad, (P.doseRef * T.psi[i] - 12) / 14.7); k = 0.6 + 0.4 * sh; }
      else {
        var e = Rr[i];
        if (!isFinite(e)) c = [44, 34, 32]; else { c = ramp(RAMPS.o2, e / 200); if ((e % 20) < 1.6) c = [c[0] * 0.6, c[1] * 0.6, c[2] * 0.6]; }
        k = 0.6 + 0.4 * sh;
      }
      var o = i * 4; d[o] = c[0] * k; d[o + 1] = c[1] * k; d[o + 2] = c[2] * k; d[o + 3] = 255;
    }
    bctx.putImageData(img, 0, 0);
  }
  var paintedKey = '';

  /* ---------- derived state ---------- */
  function fmtT(sec) { if (!isFinite(sec) || sec < 0) sec = 0; var m = Math.round(sec / 60); return Math.floor(m / 60) + ' h ' + ('0' + (m % 60)).slice(-2) + ' min'; }
  function fmtM(m, ru) { return m < 1000 ? Math.round(m / 10) * 10 + (ru ? ' м' : ' m') : (m / 1000).toFixed(2) + (ru ? ' км' : ' km'); }
  function stepLen(a, b) { return (a % W !== b % W && ((a / W) | 0) !== ((b / W) | 0)) ? CELL * 1.41421356 : CELL; }
  function straightHaz(a, b) {
    var x0 = a % W, y0 = (a / W) | 0, x1 = b % W, y1 = (b / W) | 0, n = Math.max(Math.abs(x1 - x0), Math.abs(y1 - y0)), c = 0;
    for (var i = 1; i < n; i++) { var x = Math.round(x0 + (x1 - x0) * i / n), y = Math.round(y0 + (y1 - y0) * i / n), q = y * W + x; if (T.hazard[q] || T.block[q] || T.cls[q]) c++; }
    return { haz: c, dist: Math.hypot(x1 - x0, y1 - y0) * CELL };
  }
  function pathLen(cells) { var L = 0; for (var i = 1; i < cells.length; i++) { var a = cells[i - 1], b = cells[i]; L += Math.hypot(stepLen(a, b), T.z[b] - T.z[a]); } return L; }
  // Resources for a route: oxygen, time, length, step count and how much of it lies on caution ground.
  function routeStats(cells, v) {
    var o2 = 0, t = 0, len = 0, caution = 0;
    for (var i = 1; i < cells.length; i++) {
      var e = C.edge(T, cells[i - 1], cells[i], stepLen(cells[i - 1], cells[i]), v);
      o2 += e.o2; t += e.t; len += e.len; if (T.cls[cells[i]] === 1) caution += e.len;
    }
    return { o2: o2, t: t, len: len, caution: caution, steps: Math.round(len / P.stride) };
  }
  // Warnings along a route: stretches of sand or steep ground, reported caution zones, and close passes by red zones.
  function routeWarn(cells) {
    var res = [], cum = [0], run = null, i;
    function cause(q) {
      if (T.cls[q] !== 1) return '';
      var zi = T.zoneAt[q];
      if (zi >= 0) return zones[zi].level === 3 ? '' : 'z' + zi;
      return T.I[q] < 200 ? 'sand' : 'slope';
    }
    for (i = 1; i < cells.length; i++) {
      var a = cells[i - 1], b = cells[i], l = Math.hypot(stepLen(a, b), T.z[b] - T.z[a]), k = cause(b);
      if (run && run.k !== k) { res.push(run); run = null; }
      if (k) { if (!run) run = { k: k, at: cum[i - 1], len: 0 }; run.len += l; }
      cum.push(cum[i - 1] + l);
    }
    if (run) res.push(run);
    var merged = [];
    res.forEach(function (w) {
      var prev = null; for (var j = merged.length - 1; j >= 0; j--) if (merged[j].k === w.k) { prev = merged[j]; break; }
      if (prev && w.at - (prev.at + prev.len) < 200) prev.len = w.at + w.len - prev.at; else merged.push(w);
    });
    merged = merged.filter(function (w) { return w.k.charAt(0) === 'z' || w.len >= 60; });
    zones.forEach(function (z, zi) {
      if (z.level !== 3) return;
      var best = 1e9, at = 0;
      for (var j = 0; j < cells.length; j++) { var g = Math.hypot((cells[j] % W) - z.x, ((cells[j] / W) | 0) - z.y) - z.r; if (g < best) { best = g; at = cum[j]; } }
      if (best <= 4) merged.push({ k: 'n' + zi, at: at, gap: Math.max(0, best) * CELL });
    });
    merged.sort(function (a, b) { return a.at - b.at; });
    return merged;
  }
  function warnText(w, ru) {
    var c = w.k.charAt(0), z = (c === 'z' || c === 'n') ? zones[+w.k.slice(1)] : null, at = fmtM(w.at, ru);
    if (c === 'n' && z) return ru ? 'Через ' + at + ': красная зона рядом, ' + fmtM(w.gap, 1) + ' до границы. ' + HZ_TYPES[z.type][1] + ', сообщил ' + z.by + '.'
      : 'After ' + at + ': passes ' + fmtM(w.gap) + ' from a red zone. ' + HZ_TYPES[z.type][0] + ', reported by ' + z.by + '.';
    if (c === 'z' && z) return ru ? 'Через ' + at + ': зона осторожности, ' + fmtM(w.len, 1) + '. ' + HZ_TYPES[z.type][1] + ', сообщил ' + z.by + '.'
      : 'After ' + at + ': caution zone for ' + fmtM(w.len) + '. ' + HZ_TYPES[z.type][0] + ', reported by ' + z.by + '.';
    if (w.k === 'sand') return ru ? 'Через ' + at + ': рыхлый песок, ' + fmtM(w.len, 1) + '.' : 'After ' + at + ': loose sand for ' + fmtM(w.len) + '.';
    return ru ? 'Через ' + at + ': склон круче 18°, ' + fmtM(w.len, 1) + '.' : 'After ' + at + ': slope over 18° for ' + fmtM(w.len) + '.';
  }
  function sameRoute(a, b) { if (a.length !== b.length) return false; for (var i = 0; i < a.length; i++) if (a[i] !== b[i]) return false; return true; }

  function derive() {
    ensureOut();
    var Rf = R(), b = C.beta(st.F), p = st.V * P.psiPerL, budget = b * (st.V - P.reserve);
    var retHere = Rf.dist[st.pos], margin = budget - retHere;
    var breach = (st.V - retHere <= P.reserve) || p < 650;
    var D = { beta: b, p: p, budget: budget, retHere: retHere, margin: margin, breach: breach, shed: shedMode(), tgt: null };
    if (st.target >= 0) {
      var route = null, alt = null;
      if (st.walking && st.walkCells) route = st.walkCells.slice(Math.max(0, st.walkIdx - 1));
      else if (isFinite(out.dist[st.target])) {
        var fast = C.chain(out.prev, st.target, st.pos).reverse(), sf = fast;
        ensureSafe();
        if (isFinite(safe.dist[st.target])) sf = C.chain(safe.prev, st.target, st.pos).reverse();
        route = mode === 'safe' ? sf : fast; alt = mode === 'safe' ? fast : sf;
        if (sameRoute(route, alt)) alt = null;
      }
      if (route) {
        var s = routeStats(route, speed()), er = Rf.dist[st.target], sh = straightHaz(st.pos, st.target);
        D.tgt = { eo: s.o2, er: er, after: b * (st.V - s.o2 - P.reserve) - er, route: route, back: C.chain(Rf.prev, st.target, T.hab), len: s.len, time: s.t, steps: s.steps,
          caution: s.caution, haz: sh.haz, straight: sh.dist, warns: routeWarn(route), alt: alt, altStats: alt ? routeStats(alt, speed()) : null };
        D.tgt.ok = D.tgt.after >= 0;
      }
    }
    if (breach) D.code = 'PNR_BREACH_ABORT';
    else if (D.shed) D.code = 'COGNITIVE_SHED';
    else if (D.tgt && D.tgt.haz > 0) D.code = 'SAND_TRAP_WARN';
    else if (p > 1200 && margin > 80) D.code = 'NOMINAL_EVA';
    else D.code = 'CAUTION';
    D.egress = C.chain(Rf.prev, st.pos, T.hab);
    return D;
  }
  var STATE_TEXT = {
    NOMINAL_EVA: ['ok', 'All layers available. Survival contour drawn in green. Science targets ranked by return on oxygen.'],
    CAUTION: ['warn', 'Tank is at or below 1,200 psia, or the margin is under 80 L. Science still allowed inside the contour.'],
    SAND_TRAP_WARN: ['warn', 'The straight course crosses loose sand, slopes over 18° or a reported hazard zone. Route redrawn over safer ground.'],
    COGNITIVE_SHED: ['shed', 'Fatigue above 0.55. Contour shrunk by β, science layers hidden, pace cut to 0.65 m/s. Follow the golden line.'],
    PNR_BREACH_ABORT: ['bad', 'Past the point of no return or tank below 650 psia. Science targets locked. Return to the airlock now.']
  };

  /* ---------- drawing ---------- */
  function cx(i) { return ((i % W) + 0.5) * S; } function cy(i) { return (((i / W) | 0) + 0.5) * S; }
  function strokePath(cells, color, width, dash) {
    if (cells.length < 2) return;
    ctx.save(); ctx.lineJoin = 'round'; ctx.lineCap = 'round'; ctx.setLineDash(dash || []);
    ctx.beginPath(); ctx.moveTo(cx(cells[0]), cy(cells[0]));
    for (var i = 1; i < cells.length; i++) ctx.lineTo(cx(cells[i]), cy(cells[i]));
    ctx.strokeStyle = 'rgba(20,10,6,.75)'; ctx.lineWidth = width + 3; ctx.stroke();
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.stroke(); ctx.restore();
  }
  // the active route: white on safe ground, yellow where it crosses caution ground
  function strokeRoute(cells, ok) {
    if (cells.length < 2) return;
    if (!ok) { strokePath(cells, '#ef5a45', 3.5); return; }
    strokePath(cells, '#ffffff', 3.5);
    ctx.save(); ctx.lineJoin = 'round'; ctx.lineCap = 'round'; ctx.strokeStyle = '#f2c14e'; ctx.lineWidth = 3.5; ctx.beginPath();
    for (var i = 1; i < cells.length; i++) if (T.cls[cells[i]] === 1) { ctx.moveTo(cx(cells[i - 1]), cy(cells[i - 1])); ctx.lineTo(cx(cells[i]), cy(cells[i])); }
    ctx.stroke(); ctx.restore();
  }
  function label(txt, x, y, color, align) {
    ctx.save(); ctx.font = '600 22px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = align || 'left'; ctx.textBaseline = 'middle';
    ctx.lineWidth = 5; ctx.strokeStyle = 'rgba(20,10,6,.85)'; ctx.strokeText(txt, x, y); ctx.fillStyle = color; ctx.fillText(txt, x, y); ctx.restore();
  }
  function pin(letter, X, Y, fill) {
    ctx.save(); ctx.fillStyle = fill; ctx.strokeStyle = '#1c1410'; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.moveTo(X, Y - 16); ctx.lineTo(X - 9, Y - 30); ctx.lineTo(X + 9, Y - 30); ctx.closePath(); ctx.fill(); ctx.stroke();
    ctx.beginPath(); ctx.arc(X, Y - 42, 17, 0, 6.3); ctx.fill(); ctx.stroke();
    ctx.font = '700 22px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillStyle = '#1c1410'; ctx.fillText(letter, X, Y - 41);
    ctx.restore();
  }
  function drawZone(z, draft) {
    var X = (z.x + 0.5) * S, Y = (z.y + 0.5) * S, rgb = z.level === 3 ? '239,90,69' : '242,193,78';
    ctx.save();
    ctx.beginPath(); ctx.arc(X, Y, (z.r + 0.5) * S, 0, 6.3);
    ctx.fillStyle = 'rgba(' + rgb + ',' + (draft ? 0.18 : 0.3) + ')'; ctx.fill();
    ctx.setLineDash(draft ? [10, 8] : z.level === 1 ? [4, 6] : []); ctx.lineWidth = 3; ctx.strokeStyle = 'rgba(20,10,6,.8)'; ctx.stroke();
    ctx.lineWidth = 2; ctx.strokeStyle = draft ? '#fff' : 'rgb(' + rgb + ')'; ctx.stroke(); ctx.setLineDash([]);
    ctx.beginPath(); ctx.moveTo(X, Y - 15); ctx.lineTo(X + 14, Y + 10); ctx.lineTo(X - 14, Y + 10); ctx.closePath();
    ctx.fillStyle = 'rgb(' + rgb + ')'; ctx.fill(); ctx.lineWidth = 2.5; ctx.strokeStyle = '#1c1410'; ctx.stroke();
    ctx.font = '700 17px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillStyle = '#1c1410'; ctx.fillText('!', X, Y + 2);
    ctx.restore();
  }
  var draft = null; // hazard report being placed: {x, y}
  function draw(D) {
    var forced = D.shed || D.breach, layer = forced ? 'relief' : userLayer;
    var key = layer + (layer === 'o2' ? (D.shed ? 1 : 0) : '') + ':' + zver;
    if (key !== paintedKey) { paintBase(layer); paintedKey = key; }
    ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(baseCv, 0, 0, W * S, H * S);
    // overlay: outside-contour shading + RSL zones
    var Rf = R().dist, img = octx.createImageData(W, H), d = img.data, bud = D.budget, b = D.beta, g = new Float32Array(W * H);
    for (var i = 0; i < W * H; i++) {
      var f = Rf[i] + b * out.dist[i] - bud; if (!isFinite(f)) f = 1e3; g[i] = f;
      var o = i * 4;
      if (T.rsl[i]) { var hx = ((i % W) + ((i / W) | 0)) % 3 === 0; d[o] = 105; d[o + 1] = 200; d[o + 2] = 232; d[o + 3] = hx ? 200 : 70; }
      else if (f > 0) { d[o] = 30; d[o + 1] = 6; d[o + 2] = 4; d[o + 3] = layer === 'safety' ? 95 : 150; }
    }
    octx.putImageData(img, 0, 0);
    ctx.drawImage(overCv, 0, 0, W * S, H * S);
    // crew hazard reports are shown on every layer
    zones.forEach(function (z) { drawZone(z, false); });
    if (draft) drawZone({ x: draft.x, y: draft.y, r: +$('hz-r').value, level: +$('hz-level').value }, true);
    // survival contour (marching squares)
    ctx.save(); ctx.beginPath();
    for (var y = 0; y < H - 1; y++) for (var x = 0; x < W - 1; x++) {
      var a = g[y * W + x], bb = g[y * W + x + 1], c = g[(y + 1) * W + x + 1], dd = g[(y + 1) * W + x];
      var ia = a <= 0, ib = bb <= 0, ic = c <= 0, id = dd <= 0;
      if (ia === ib && ib === ic && ic === id) continue;
      var pts = [];
      if (ia !== ib) pts.push([x + 0.5 + a / (a - bb), y + 0.5]);
      if (ib !== ic) pts.push([x + 1.5, y + 0.5 + bb / (bb - c)]);
      if (ic !== id) pts.push([x + 0.5 + dd / (dd - c), y + 1.5]);
      if (id !== ia) pts.push([x + 0.5, y + 0.5 + a / (a - dd)]);
      for (var q = 0; q + 1 < pts.length; q += 2) { ctx.moveTo(pts[q][0] * S, pts[q][1] * S); ctx.lineTo(pts[q + 1][0] * S, pts[q + 1][1] * S); }
    }
    ctx.lineCap = 'round'; ctx.strokeStyle = 'rgba(12,30,18,.8)'; ctx.lineWidth = 6; ctx.stroke();
    ctx.strokeStyle = D.breach ? '#ef5a45' : '#7fe0a2'; ctx.lineWidth = 3; ctx.stroke(); ctx.restore();
    // routes
    if (D.tgt && !forced) {
      strokePath(D.tgt.back, '#f0b24a', 2.5, [8, 8]);
      if (D.tgt.alt) strokePath(D.tgt.alt, '#9fb9d0', 2, [3, 9]);
      strokeRoute(D.tgt.route, D.tgt.ok);
    }
    if (st.pos !== T.hab) strokePath(D.egress, '#f0b24a', forced ? 6 : 3.5);
    // science targets
    if (!forced) SCI.forEach(function (s) {
      var X = cx(s.cell), Y = cy(s.cell);
      ctx.save(); ctx.translate(X, Y); ctx.rotate(Math.PI / 4);
      ctx.fillStyle = s.done ? 'rgba(60,48,42,.9)' : '#1c1410'; ctx.strokeStyle = s.done ? '#8a7667' : '#f3e6d6'; ctx.lineWidth = 2.5;
      ctx.fillRect(-16, -16, 32, 32); ctx.strokeRect(-16, -16, 32, 32); ctx.restore();
      ctx.save(); ctx.font = '700 18px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillStyle = s.done ? '#8a7667' : '#f3e6d6'; ctx.fillText(s.id, X, Y + 1); ctx.restore();
    });
    // custom waypoint
    if (st.target >= 0 && !st.sci && !forced) {
      ctx.save(); ctx.strokeStyle = '#fff'; ctx.lineWidth = 2.5; ctx.beginPath(); ctx.arc(cx(st.target), cy(st.target), 13, 0, 6.3); ctx.stroke(); ctx.restore();
    }
    // airlock
    var hx2 = cx(T.hab), hy2 = cy(T.hab);
    ctx.save(); ctx.fillStyle = '#f0b24a'; ctx.strokeStyle = '#1c1410'; ctx.lineWidth = 3;
    ctx.beginPath(); ctx.rect(hx2 - 14, hy2 - 14, 28, 28); ctx.fill(); ctx.stroke(); ctx.restore();
    label('AIRLOCK', hx2, hy2 + 32, '#f0b24a', 'center');
    // astronaut
    var ax = cx(st.pos), ay = cy(st.pos);
    ctx.save(); ctx.beginPath(); ctx.arc(ax, ay, 20, 0, 6.3); ctx.strokeStyle = 'rgba(255,255,255,.55)'; ctx.lineWidth = 2; ctx.stroke();
    ctx.beginPath(); ctx.arc(ax, ay, 10, 0, 6.3); ctx.fillStyle = '#fff'; ctx.fill(); ctx.strokeStyle = '#1c1410'; ctx.lineWidth = 3; ctx.stroke(); ctx.restore();
    if (st.pos !== T.hab) label('EV-1', ax + 28, ay, '#fff');
    // route pins: A is where the walk starts (the astronaut), B is the destination
    if (st.target >= 0 && !forced) { pin('B', cx(st.target), cy(st.target), '#f0b24a'); if (!st.walking) pin('A', ax, ay, '#69c8e8'); }
    // highlight ring for "nearest hazard" and "locate"
    if (st.ping >= 0 && Date.now() - st.pingT < 6000) {
      ctx.save(); ctx.strokeStyle = '#fff'; ctx.lineWidth = 4; ctx.beginPath(); ctx.arc(cx(st.ping), cy(st.ping), 30, 0, 6.3); ctx.stroke();
      ctx.strokeStyle = '#69c8e8'; ctx.lineWidth = 3; ctx.setLineDash([6, 6]); ctx.beginPath(); ctx.arc(cx(st.ping), cy(st.ping), 44, 0, 6.3); ctx.stroke(); ctx.restore();
    }
    // map furniture: scale bar and north arrow
    var bx = 28, by = H * S - 30, km = 1000 / CELL * S;
    ctx.save(); ctx.strokeStyle = '#f3e6d6'; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(bx, by - 7); ctx.lineTo(bx, by); ctx.lineTo(bx + km, by); ctx.lineTo(bx + km, by - 7); ctx.stroke(); ctx.restore();
    label('1 km', bx + km + 10, by - 3, '#f3e6d6');
    label('N ↑', W * S - 28, 32, '#f3e6d6', 'right');
    label('CANYON WALL', W * S / 2, 32, 'rgba(243,230,214,.8)', 'center');
  }

  /* ---------- AR visor preview: first-person view along the active route (concept for a helmet display) ---------- */
  var ar = $('ar'), actx = ar ? ar.getContext('2d') : null;
  function destName(ru) {
    if (st.target < 0) return '';
    if (st.target === T.hab) return ru ? 'база (шлюз)' : 'Base airlock';
    if (st.sci) return ru ? st.sci.ru : st.sci.id + ' · ' + st.sci.name;
    return ru ? 'точка на карте' : 'Waypoint';
  }
  function drawAR(D) {
    if (!actx) return;
    var w = ar.width, h = ar.height, hz = h * 0.4, forced = D.shed || D.breach;
    var toB = D.tgt && !forced, path = toB ? D.tgt.route : D.egress;
    var px = st.pos % W, py = (st.pos / W) | 0, aim = path[Math.min(path.length - 1, 5)];
    var fx = (aim % W) - px, fy = ((aim / W) | 0) - py, fl = Math.hypot(fx, fy);
    if (fl < 0.5) { fx = 0; fy = -1; fl = 1; } fx /= fl; fy /= fl;
    function sx(l, f) { return w / 2 + l / f * w * 0.6; } function sy(f) { return hz + (h - hz) * 0.9 / f; }
    var sky = actx.createLinearGradient(0, 0, 0, hz); sky.addColorStop(0, '#20130d'); sky.addColorStop(1, '#a76a48');
    actx.fillStyle = sky; actx.fillRect(0, 0, w, hz);
    var gr = actx.createLinearGradient(0, hz, 0, h); gr.addColorStop(0, '#6a4230'); gr.addColorStop(1, '#2a1912');
    actx.fillStyle = gr; actx.fillRect(0, hz, w, h - hz);
    // hazard cells ahead, drawn as patches on the ground
    for (var f = 28; f >= 1; f--) for (var l = -16; l <= 16; l++) {
      var gxc = Math.round(px + fx * f - fy * l), gyc = Math.round(py + fy * f + fx * l);
      if (gxc < 0 || gyc < 0 || gxc >= W || gyc >= H) continue;
      var cl = T.cls[gyc * W + gxc]; if (!cl) continue;
      var f0 = f - 0.5, f1 = f + 0.5;
      actx.fillStyle = cl === 2 ? 'rgba(239,90,69,.6)' : 'rgba(242,193,78,.38)';
      actx.beginPath(); actx.moveTo(sx(l - 0.5, f0), sy(f0)); actx.lineTo(sx(l + 0.5, f0), sy(f0)); actx.lineTo(sx(l + 0.5, f1), sy(f1)); actx.lineTo(sx(l - 0.5, f1), sy(f1)); actx.closePath(); actx.fill();
    }
    // the route as a line on the ground
    actx.save(); actx.lineJoin = 'round'; actx.lineCap = 'round'; actx.beginPath(); actx.moveTo(w / 2, h);
    var last = null;
    for (var i = 1; i < path.length && i < 60; i++) {
      var rx = (path[i] % W) - px, ry = ((path[i] / W) | 0) - py, pf = rx * fx + ry * fy, pl = -rx * fy + ry * fx;
      if (pf < 0.6) continue; last = [sx(pl, pf), sy(pf)]; actx.lineTo(last[0], last[1]);
    }
    actx.strokeStyle = 'rgba(20,10,6,.7)'; actx.lineWidth = 12; actx.stroke();
    actx.strokeStyle = toB ? '#ffffff' : '#f0b24a'; actx.lineWidth = 6; actx.stroke(); actx.restore();
    if (last && path.length > 1) { // beacon at the far end of the visible route
      actx.save(); actx.strokeStyle = toB ? '#f0b24a' : '#f0b24a'; actx.lineWidth = 3; actx.beginPath(); actx.moveTo(last[0], last[1]); actx.lineTo(last[0], last[1] - 46); actx.stroke();
      actx.fillStyle = '#f0b24a'; actx.beginPath(); actx.arc(last[0], last[1] - 54, 9, 0, 6.3); actx.fill(); actx.restore();
    }
    // visor frame and readouts
    actx.save(); actx.strokeStyle = 'rgba(127,224,162,.8)'; actx.lineWidth = 2;
    [[14, 14, 1, 1], [w - 14, 14, -1, 1], [14, h - 14, 1, -1], [w - 14, h - 14, -1, -1]].forEach(function (c) { actx.beginPath(); actx.moveTo(c[0] + 26 * c[2], c[1]); actx.lineTo(c[0], c[1]); actx.lineTo(c[0], c[1] + 26 * c[3]); actx.stroke(); });
    actx.font = '600 17px "IBM Plex Mono", ui-monospace, monospace'; actx.textBaseline = 'top'; actx.shadowColor = '#000'; actx.shadowBlur = 4;
    actx.fillStyle = '#f3e6d6'; actx.textAlign = 'left';
    var line1 = toB ? 'B ' + destName() : (st.pos === T.hab ? 'AT THE AIRLOCK' : 'AIRLOCK'), line2 = '';
    if (toB) line2 = fmtM(D.tgt.len) + ' · ' + Math.round(D.tgt.time / 60) + ' min · ' + D.tgt.steps.toLocaleString('en-US') + ' steps';
    else if (st.pos !== T.hab) line2 = fmtM(pathLen(D.egress));
    actx.fillText(line1, 28, 26); actx.fillStyle = '#bba591'; actx.fillText(line2, 28, 50);
    actx.textAlign = 'right'; actx.fillStyle = D.margin > 80 ? '#7fe0a2' : D.margin > 0 ? '#f2c14e' : '#ef5a45';
    actx.fillText('O₂ margin ' + (isFinite(D.margin) ? D.margin.toFixed(0) : '—') + ' L', w - 28, 26);
    var ban = 'PATH CLEAR', col = '#7fe0a2';
    if (D.breach) { ban = 'RETURN TO AIRLOCK NOW'; col = '#ef5a45'; }
    else if (toB && !D.tgt.ok) { ban = 'BEYOND POINT OF NO RETURN'; col = '#ef5a45'; }
    else if (toB && D.tgt.warns.length) { ban = warnText(D.tgt.warns[0]).replace(/^After /, 'In ').toUpperCase(); col = D.tgt.warns[0].k.charAt(0) === 'n' ? '#ef5a45' : '#f2c14e'; }
    else if (T.cls[st.pos] === 1) { ban = 'CAUTION GROUND UNDER FOOT'; col = '#f2c14e'; }
    actx.font = '600 15px "IBM Plex Mono", ui-monospace, monospace'; actx.textAlign = 'center'; actx.textBaseline = 'middle';
    while (ban.length > 62) ban = ban.slice(0, ban.lastIndexOf(' ', 60)) + '…';
    var bw = Math.min(w - 40, actx.measureText(ban).width + 28);
    actx.shadowBlur = 0; actx.fillStyle = 'rgba(20,12,8,.82)'; actx.fillRect(w / 2 - bw / 2, h - 52, bw, 30);
    actx.strokeStyle = col; actx.lineWidth = 1.5; actx.strokeRect(w / 2 - bw / 2, h - 52, bw, 30);
    actx.fillStyle = col; actx.fillText(ban, w / 2, h - 36, w - 60); actx.restore();
  }

  /* ---------- panel ---------- */
  function set(id, v) { var e = $(id); if (e && e.textContent !== v) e.textContent = v; }
  function el(tag, cls, txt) { var e = document.createElement(tag); if (cls) e.className = cls; if (txt != null) e.textContent = txt; return e; }
  var LAND = [
    ['Open canyon floor', 'esp_060148_1610.jpg', 'ESP_060148_1610'], ['Canyon wall', 'psp_004054_1675.jpg', 'PSP_004054_1675'],
    ['Talus apron under the wall', 'psp_004054_1675.jpg', 'PSP_004054_1675'], ['Layered sulfate mound', 'esp_037494_1685.jpg', 'ESP_037494_1685'],
    ['Block field', 'psp_004384_1705.jpg', 'PSP_004384_1705'], ['Fault trough (graben)', 'esp_034883_1645.jpg', 'ESP_034883_1645'],
    ['Dark dunes, loose sand', 'esp_044892_1695.jpg', 'ESP_044892_1695']
  ];
  var shotsKey = '', warnsKey = '', zlistKey = '';
  function photos(D) {
    var box = $('photos'), list = [], forced = D.shed || D.breach;
    if (D.tgt && !forced) {
      var cnt = {}, order = [], dl = T.land[st.target];
      D.tgt.route.forEach(function (c) { var k = T.land[c]; if (!cnt[k]) { cnt[k] = 0; order.push(k); } cnt[k]++; });
      list.push([dl, 'Destination']);
      order.forEach(function (k) { if (k !== dl && cnt[k] >= 6 && list.length < 3) list.push([k, 'On the way']); });
    }
    var key = list.join('|'); if (key === shotsKey) return; shotsKey = key;
    box.textContent = '';
    if (!list.length) { box.appendChild(el('p', 'note', 'Build a route to see what the ground looks like at the destination and on the way.')); return; }
    list.forEach(function (it) {
      var L = LAND[it[0]], fig = el('figure'), im = el('img'), cap = el('figcaption');
      im.src = 'img/' + L[1]; im.alt = L[0] + ', seen from orbit'; im.loading = 'lazy';
      cap.appendChild(el('b', '', it[1] + ': ' + L[0])); cap.appendChild(el('span', '', 'HiRISE ' + L[2] + ', orbital analogue of this landform'));
      fig.appendChild(im); fig.appendChild(cap); box.appendChild(fig);
    });
  }
  function warnsList(D) {
    var ul = $('warns'), forced = D.shed || D.breach, items = [];
    if (D.tgt && !forced) {
      if (!D.tgt.ok) items.push(['bad', 'Destination is beyond the point of no return.']);
      D.tgt.warns.slice(0, 6).forEach(function (w) { items.push([w.k.charAt(0) === 'n' ? 'bad' : 'warn', warnText(w)]); });
      if (D.tgt.warns.length > 6) items.push(['warn', 'And ' + (D.tgt.warns.length - 6) + ' more further along.']);
      if (!items.length) items.push(['ok', 'No obstacles or hazard zones on this route.']);
    } else items.push(['', forced ? 'Only the return route is shown in this HUD state.' : 'No route yet. Choose a destination B.']);
    var key = JSON.stringify(items); if (key === warnsKey) return; warnsKey = key;
    ul.textContent = '';
    items.forEach(function (it) { ul.appendChild(el('li', it[0] ? 'c-' + it[0] : 'muted', it[1])); });
  }
  function ago(t) {
    if (!t) return 'sample';
    var m = Math.max(0, Math.round((Date.now() - t) / 60000));
    return m < 1 ? 'just now' : m < 60 ? m + ' min ago' : m < 1440 ? Math.round(m / 60) + ' h ago' : Math.round(m / 1440) + ' d ago';
  }
  function zoneList() {
    var note = dbState === 'local' ? 'Saved in this browser and shared live between open tabs. Set AEROMARGIN_DB_URL in index.html to share between devices.'
      : dbState === 'online' ? 'Shared database connected. Reports from every crew member appear here.'
      : dbState === 'offline' ? 'Shared database unreachable. Showing the local copy; new reports are sent when it is back.' : 'Connecting to the shared database…';
    set('db-note', note);
    var key = zver + ':' + zones.map(function (z) { return z.id + (z.local ? 'L' : ''); }).join(',') + ':' + Math.floor(Date.now() / 60000);
    if (key === zlistKey) return; zlistKey = key;
    var box = $('zlist'); box.textContent = '';
    if (!zones.length) { box.appendChild(el('p', 'note', 'No hazard reports yet.')); return; }
    zones.slice().reverse().forEach(function (z) {
      var row = el('div', 'zrow'), main = el('div', 'zmain');
      row.appendChild(el('i', 'dot ' + (z.level === 3 ? 'b-bad' : 'b-warn')));
      main.appendChild(el('b', '', HZ_TYPES[z.type][0] + ' · ' + LEVELS[z.level]));
      main.appendChild(el('span', 'muted', z.by + ' · ' + ago(z.t) + ' · radius ' + z.r * CELL + ' m' + (z.local && DB_URL ? ' · not synced yet' : '')));
      if (z.desc) main.appendChild(el('p', '', z.desc));
      row.appendChild(main);
      if (z.photo) { var im = el('img'); im.src = z.photo; im.alt = 'Photo attached to the report'; row.appendChild(im); }
      var bt = el('div', 'zbtn'), loc = el('button', '', 'Locate'), rm = el('button', '', 'Remove');
      loc.type = rm.type = 'button';
      loc.addEventListener('click', function () { ping(z.y * W + z.x); cv.scrollIntoView({ block: 'center', behavior: reduce ? 'auto' : 'smooth' }); });
      rm.addEventListener('click', function () {
        if (st.walking) return;
        if (DB_URL && !z.local) remote('DELETE', '/hazards/' + z.id).catch(function () { dbState = 'offline'; zoneList(); });
        setZones(zones.filter(function (q) { return q !== z; }));
      });
      bt.appendChild(loc); bt.appendChild(rm); row.appendChild(bt); box.appendChild(row);
    });
  }
  function ping(cell) { st.ping = cell; st.pingT = Date.now(); update(); setTimeout(function () { if (!st.walking) update(); }, 6100); }
  function press(sel, attr, val) {
    Array.prototype.forEach.call(document.querySelectorAll(sel), function (bn) { bn.setAttribute('aria-pressed', bn.getAttribute(attr) === val ? 'true' : 'false'); });
  }
  function sciAt(cell) { for (var i = 0; i < SCI.length; i++) if (SCI[i].cell === cell) return SCI[i]; return null; }

  function panel(D) {
    var stt = STATE_TEXT[D.code];
    var vis = $('visor'); vis.className = 'visor s-' + stt[0];
    set('st-code', D.code); set('st-text', stt[1]);
    set('g-psi', Math.round(D.p).toLocaleString('en-US')); set('g-use', st.V.toFixed(0));
    var band = D.p >= 1200 ? 'ok' : D.p >= 650 ? 'warn' : 'bad';
    var bar = $('tankbar'); bar.style.width = Math.max(2, (D.p - 0) / 2450 * 100) + '%'; bar.className = 'fill b-' + band;
    set('g-flow', st.ema.toFixed(2)); set('g-time', fmtT((st.V - P.reserve) / st.ema * 60));
    set('g-ret', isFinite(D.retHere) ? D.retHere.toFixed(0) : '—'); set('g-margin', (D.margin >= 0 ? '+' : '−') + Math.abs(D.margin).toFixed(0));
    $('g-margin').className = 'num ' + (D.margin > 80 ? 'c-ok' : D.margin > 0 ? 'c-warn' : 'c-bad');
    set('g-beta', D.beta.toFixed(2)); set('g-fat', st.F.toFixed(2));
    set('g-el', fmtT(st.elapsed)); set('g-dose', st.dose.toFixed(1)); set('g-walk', (st.walked / 1000).toFixed(2));
    set('tank-out', Math.round(D.p).toLocaleString('en-US') + ' psia'); set('fat-out', st.F.toFixed(2) + ' → β ' + D.beta.toFixed(2));
    // visor arrow: bearing along the egress line
    var eg = D.egress, aim = eg[Math.min(eg.length - 1, 8)], dx = (aim % W) - (st.pos % W), dy = ((aim / W) | 0) - ((st.pos / W) | 0);
    var brg = (Math.atan2(dx, -dy) * 180 / Math.PI + 360) % 360, home = pathLen(eg);
    $('arrow').style.transform = 'rotate(' + (st.pos === T.hab ? 0 : brg) + 'deg)';
    set('v-home', st.pos === T.hab ? 'At the airlock' : 'Airlock ' + (home / 1000).toFixed(2) + ' km · bearing ' + ('00' + Math.round(brg)).slice(-3) + '°');
    set('v-o2', st.pos === T.hab ? 'Survival contour shown in green on the map' : 'O₂ home ' + D.retHere.toFixed(0) + ' L · margin ' + D.margin.toFixed(0) + ' L');
    // layers
    var forced = D.shed || D.breach, layer = forced ? 'relief' : userLayer, L = LEG[layer];
    Array.prototype.forEach.call(document.querySelectorAll('[data-layer]'), function (bn) {
      bn.disabled = forced; bn.setAttribute('aria-pressed', bn.getAttribute('data-layer') === layer ? 'true' : 'false');
    });
    set('leg-name', L[0] + (forced ? ' (other layers switched off by the HUD)' : '')); set('leg-lo', L[1]); set('leg-hi', L[2]);
    $('leg-bar').style.background = 'linear-gradient(90deg,' + L[3].map(function (c) { return 'rgb(' + c.join(',') + ')'; }).join(',') + ')';
    // route planner controls
    var ps = sciAt(st.pos);
    $('selA').value = st.pos === T.hab ? 'hab' : ps ? ps.id : 'pos';
    $('selB').value = st.target < 0 ? '' : st.target === T.hab ? 'hab' : st.sci ? st.sci.id : 'pt';
    $('selA').disabled = $('selB').disabled = $('swap').disabled = st.walking;
    press('[data-click]', 'data-click', clickMode); press('[data-mode]', 'data-mode', mode);
    $('hz-add').setAttribute('aria-pressed', clickMode === 'haz' ? 'true' : 'false');
    // target card
    var tc = $('tgt');
    if (D.tgt && !forced) {
      var t = D.tgt, nm = destName();
      set('t-name', nm); set('t-note', st.sci ? st.sci.note + (st.sci.done ? ' (sampled)' : '') : st.target === T.hab ? 'The habitat airlock, the research base of this EVA' : 'Free waypoint picked on the map');
      set('t-len', (t.len / 1000).toFixed(2) + ' km'); set('t-out', t.eo.toFixed(0) + ' L'); set('t-back', t.er.toFixed(0) + ' L');
      set('t-time', fmtT(t.time)); set('t-steps', '≈ ' + t.steps.toLocaleString('en-US'));
      set('t-caut', fmtM(t.caution)); $('t-caut').className = 'num ' + (t.caution > 0 ? 'c-warn' : 'c-ok');
      set('t-after', (t.after >= 0 ? '+' : '−') + Math.abs(t.after).toFixed(0) + ' L');
      $('t-after').className = 'num ' + (t.after > 80 ? 'c-ok' : t.after >= 0 ? 'c-warn' : 'c-bad');
      var detour = t.len - t.straight;
      set('t-verdict', !t.ok ? 'Beyond the point of no return. Walking there is blocked.' :
        (t.haz > 0 ? 'Straight line crosses ' + t.haz + ' hazard cells. Detour adds ' + Math.max(0, detour).toFixed(0) + ' m and keeps to firmer ground.' : 'Inside the survival contour. Clear to go.'));
      $('t-verdict').className = 'verdict ' + (t.ok ? (t.haz ? 'c-warn' : 'c-ok') : 'c-bad');
      var a = t.altStats, cmp;
      if (st.walking) cmp = '';
      else if (!a) cmp = 'The safest route and the shortest route are the same here.';
      else if (mode === 'safe') cmp = 'Safest route shown. It is ' + fmtM(Math.max(0, t.len - a.len)) + ' longer and costs ' + Math.max(0, t.eo - a.o2).toFixed(0) + ' L more O₂ than the shortest (dotted), but has ' + fmtM(Math.max(0, a.caution - t.caution)) + ' less caution ground.';
      else cmp = 'Shortest route shown. The safest (dotted) is ' + fmtM(Math.max(0, a.len - t.len)) + ' longer and avoids ' + fmtM(Math.max(0, t.caution - a.caution)) + ' of caution ground.';
      set('t-alt', cmp);
      tc.hidden = false; $('t-empty').hidden = true;
    } else {
      tc.hidden = true; $('t-empty').hidden = false;
      set('t-empty', forced ? (D.breach ? 'Science targets are locked. Only the return route is available.' : 'Science targets are hidden while fatigue is above 0.55.') : 'Choose a destination B: tap the map, use the list above, pick a target below, or say it.');
    }
    // readout strip under the map
    var rt = D.tgt && !forced ? D.tgt : null;
    set('r-len', rt ? fmtM(rt.len) : '—'); set('r-time', rt ? Math.max(1, Math.round(rt.time / 60)) + ' min' : '—');
    set('r-steps', rt ? rt.steps.toLocaleString('en-US') : '—'); set('r-o2', rt ? rt.eo.toFixed(0) + ' L' : '—');
    set('r-after', rt ? (rt.after >= 0 ? '+' : '−') + Math.abs(rt.after).toFixed(0) + ' L' : '—');
    $('r-after').className = rt ? (rt.after > 80 ? 'c-ok' : rt.after >= 0 ? 'c-warn' : 'c-bad') : '';
    $('go').disabled = st.walking || forced || !D.tgt || !D.tgt.ok || st.target === st.pos;
    $('home').disabled = st.walking || st.pos === T.hab;
    // science ranking
    var list = $('sci'); list.innerHTML = '';
    if (forced) { list.hidden = true; } else {
      list.hidden = false;
      var rows = SCI.map(function (s) {
        var eo = out.dist[s.cell], er = R().dist[s.cell], dI = 0.5 * Math.log(1 + s.s2 / P.sigmaNoise), dr = P.doseRef * T.psi[s.cell];
        var merj = (s.done || !isFinite(eo) || eo < 0.5) ? 0 : dI / (eo * (1 + P.lambdaRad * dr));
        return { s: s, eo: eo, merj: merj, ok: isFinite(eo) && D.beta * (st.V - eo - P.reserve) - er >= 0 };
      }).sort(function (a, b) { return b.merj - a.merj; });
      rows.forEach(function (r) {
        var bt = document.createElement('button'); bt.type = 'button'; bt.className = 'sci-row' + (st.sci === r.s ? ' on' : '');
        var tag = r.s.done ? 'sampled' : (st.pos === r.s.cell ? 'here' : (r.ok ? r.eo.toFixed(0) + ' L out' : 'locked'));
        bt.innerHTML = '<span class="badge">' + r.s.id + '</span><span class="nm">' + r.s.name + '</span><span class="mj">' + (r.merj ? (r.merj * 1000).toFixed(1) : '—') + '</span><span class="tg ' + (r.s.done ? '' : r.ok ? 'c-ok' : 'c-bad') + '">' + tag + '</span>';
        bt.addEventListener('click', function () { if (st.walking) return; st.target = r.s.cell; st.sci = r.s; st.msg = ''; update(); });
        list.appendChild(bt);
      });
    }
    set('msg', st.msg);
    warnsList(D); photos(D); zoneList();
    // telemetry packet
    var px = st.pos % W, py = (st.pos / W) | 0;
    var pkt = {
      packet_id: 'DEMO-EVA-' + ('00000' + Math.round(st.elapsed)).slice(-5), crew_id: 'EV-1',
      location: { frame: 'local grid, metres from airlock (synthetic terrain)', east_m: (px - C.HAB.x) * CELL, north_m: (C.HAB.y - py) * CELL, elevation_rel_m: +T.z[st.pos].toFixed(1), local_slope_deg: +T.slope[st.pos].toFixed(1), soil_thermal_inertia_tiu: +T.I[st.pos].toFixed(0), sky_view_factor: +T.psi[st.pos].toFixed(2), safety_class: ['SAFE', 'CAUTION', 'NO_ENTRY'][T.cls[st.pos]] },
      plss_life_support: { o2_primary_pressure_psia: +D.p.toFixed(0), o2_primary_usable_liters_stp: +st.V.toFixed(1), o2_flow_rate_lpm_ema180: +st.ema.toFixed(2), time_remaining_o2_sec: Math.max(0, Math.round((st.V - P.reserve) / st.ema * 60)) },
      biometrics_and_risk: { metabolic_power_watts: +st.W.toFixed(0), fatigue_index_f: +st.F.toFixed(2), present_bias_beta: +D.beta.toFixed(3), radiation_dose_rate_usv_hr: +(P.doseRef * T.psi[st.pos]).toFixed(1), cumulative_eva_dose_msv: +(st.dose / 1000).toFixed(4) },
      navigation_pnr: { distance_to_airlock_m: +home.toFixed(0), optimal_return_path_o2_liters: +D.retHere.toFixed(1), pnr_margin_liters: +D.margin.toFixed(1), hud_state: D.code, hud_cognitive_mode: forced ? 'EGRESS_ONLY' : 'STANDARD_5_LAYER' },
      navigation_route: D.tgt ? { mode: mode === 'safe' ? 'SAFEST' : 'SHORTEST', distance_m: +D.tgt.len.toFixed(0), eta_sec: Math.round(D.tgt.time), steps_estimate: D.tgt.steps, o2_liters: +D.tgt.eo.toFixed(1), caution_ground_m: +D.tgt.caution.toFixed(0), warnings: D.tgt.warns.length } : null,
      shared_hazard_reports: zones.length
    };
    set('packet', JSON.stringify(pkt, null, 2));
  }
  function update() { var D = derive(); draw(D); panel(D); drawAR(D); return D; }

  /* ---------- walking ---------- */
  function walk(cells, onDone, track) {
    if (cells.length < 2) return;
    st.walking = true; st.msg = ''; st.walkCells = track ? cells : null; st.walkIdx = 1;
    var v = speed(), i = 1, per = reduce ? cells.length : Math.max(1, Math.ceil(cells.length / 90));
    function step() {
      for (var n = 0; n < per && i < cells.length; n++, i++) {
        var a = cells[i - 1], b = cells[i], dd = stepLen(a, b);
        var e = C.edge(T, a, b, dd, v), al = Math.min(1, e.t / 180), fl = C.flow(e.W);
        st.V = Math.max(0, st.V - e.o2); st.elapsed += e.t; st.walked += e.len; st.W = e.W;
        st.ema = al * fl + (1 - al) * st.ema; st.dose += P.doseRef * T.psi[b] * e.t / 3600; st.pos = b;
      }
      st.walkIdx = i;
      $('tank').value = Math.round(st.V * P.psiPerL);
      if (i < cells.length) { update(); requestAnimationFrame(step); }
      else {
        st.walking = false; st.walkCells = null; onDone();
        if (pendingZones) { var pz = pendingZones; pendingZones = null; setZones(pz); } else update();
      }
    }
    step();
  }
  function goWalk() {
    var D = derive(); if (st.walking || D.shed || D.breach || !D.tgt || !D.tgt.ok || st.target === st.pos) return false;
    var sci = st.sci, home = st.target === T.hab;
    walk(D.tgt.route, function () {
      if (sci && !sci.done) {
        sci.done = true; var n = 0;
        SCI.forEach(function (o) { if (o !== sci && o.cls === sci.cls && !o.done) { o.s2 *= 0.25; n++; } });
        st.msg = 'Sample logged at ' + sci.name + '.' + (n ? ' Uncertainty at ' + n + ' similar target' + (n > 1 ? 's' : '') + ' collapsed, so the ranking changed.' : '');
      } else st.msg = home ? 'Back at the airlock with ' + st.V.toFixed(0) + ' L usable.' : 'Arrived at the waypoint.';
      st.target = -1; st.sci = null;
    }, true);
    return true;
  }
  $('go').addEventListener('click', goWalk);
  $('home').addEventListener('click', function () {
    var D = derive(); st.target = -1; st.sci = null;
    walk(D.egress, function () { st.msg = 'Back at the airlock with ' + st.V.toFixed(0) + ' L usable.'; });
  });
  $('reset').addEventListener('click', function () { if (st.walking) return; reset(); update(); });
  $('tank').addEventListener('input', function () { if (st.walking) return; st.V = +this.value / P.psiPerL; update(); });
  $('fat').addEventListener('input', function () { if (st.walking) return; st.F = +this.value / 100; update(); });
  Array.prototype.forEach.call(document.querySelectorAll('[data-layer]'), function (bn) {
    bn.addEventListener('click', function () { userLayer = bn.getAttribute('data-layer'); update(); });
  });

  /* ---------- route planner: point A, point B, route mode ---------- */
  function zoneMsg(c) { var z = zones[T.zoneAt[c]]; return z ? HZ_TYPES[z.type][0] + ', reported by ' + z.by : 'reported hazard'; }
  function pickDest(c, sci) { // set destination B; returns true when a route can be planned
    if (st.walking) return false;
    var D = derive();
    if (D.shed || D.breach) { st.msg = 'Targets are locked in this HUD state. Lower fatigue or use Return to airlock.'; return false; }
    if (sci) { st.target = sci.cell; st.sci = sci; st.msg = ''; return true; }
    if (c === st.pos) { st.msg = 'That is where you already are.'; return false; }
    if (T.rsl[c]) { st.msg = 'Restricted: slope with seasonal streaks (planetary protection).'; return false; }
    if (T.uz[c]) { st.msg = 'Red zone (' + zoneMsg(c) + '). Pick a destination outside it.'; return false; }
    if (T.block[c]) { st.msg = 'Too steep to walk: slope over 30°.'; return false; }
    if (!isFinite(out.dist[c])) { st.msg = 'No walkable route to that point.'; return false; }
    st.target = c; st.sci = sciAt(c); st.msg = ''; return true;
  }
  function setStart(c) { // move the start point A (the astronaut's planning position)
    if (st.walking) return false;
    if (T.rsl[c] || T.block[c]) { st.msg = 'Start point A cannot be inside a red zone.'; return false; }
    if (!isFinite(R().dist[c])) { st.msg = 'No walkable way home from that point.'; return false; }
    st.pos = c; if (st.target === c) { st.target = -1; st.sci = null; }
    st.msg = 'Start point A set. EV-1 is placed there for planning.'; return true;
  }
  function fillSelects() {
    var a = $('selA'), b = $('selB');
    function opt(sel, v, t) { var o = el('option', '', t); o.value = v; sel.appendChild(o); }
    opt(a, 'hab', 'Base airlock'); opt(b, '', 'Choose…'); opt(b, 'hab', 'Base airlock');
    SCI0.forEach(function (s) { opt(a, s.id, s.id + ' ' + s.name); opt(b, s.id, s.id + ' ' + s.name); });
    opt(a, 'pos', 'Point on the map'); opt(b, 'pt', 'Point on the map');
    var ty = $('hz-type'); Object.keys(HZ_TYPES).forEach(function (k) { opt(ty, k, HZ_TYPES[k][0]); });
  }
  function sciById(id) { for (var i = 0; i < SCI.length; i++) if (SCI[i].id === id) return SCI[i]; return null; }
  $('selA').addEventListener('change', function () {
    var v = this.value;
    if (v === 'pos') { clickMode = 'A'; st.msg = 'Tap the map to set start point A.'; }
    else setStart(v === 'hab' ? T.hab : sciById(v).cell);
    update();
  });
  $('selB').addEventListener('change', function () {
    var v = this.value;
    if (v === '') { st.target = -1; st.sci = null; st.msg = ''; }
    else if (v === 'pt') { clickMode = 'B'; st.msg = 'Tap the map to set destination B.'; }
    else if (v === 'hab') pickDest(T.hab, null);
    else pickDest(-1, sciById(v));
    update();
  });
  $('swap').addEventListener('click', function () {
    if (st.walking || st.target < 0) { st.msg = 'Set a destination B first, then swap.'; update(); return; }
    var a = st.pos, b = st.target;
    if (setStart(b)) { st.target = a; st.sci = sciAt(a); st.msg = 'A and B swapped.'; }
    update();
  });
  Array.prototype.forEach.call(document.querySelectorAll('[data-click]'), function (bn) {
    bn.addEventListener('click', function () { closeForm(); clickMode = bn.getAttribute('data-click'); st.msg = 'Tap the map to set ' + (clickMode === 'A' ? 'start point A.' : 'destination B.'); update(); });
  });
  Array.prototype.forEach.call(document.querySelectorAll('[data-mode]'), function (bn) {
    bn.addEventListener('click', function () { if (st.walking) return; mode = bn.getAttribute('data-mode'); update(); });
  });

  /* ---------- reporting a hazard zone ---------- */
  var draftPhoto = '';
  function openForm(c) {
    draft = { x: c % W, y: (c / W) | 0 };
    $('hzform').hidden = false;
    set('hz-where', 'Location: ' + ((draft.x - C.HAB.x) * CELL) + ' m east, ' + ((C.HAB.y - draft.y) * CELL) + ' m north of the airlock. Tap the map again to move it.');
    set('hz-r-out', +$('hz-r').value * CELL + ' m');
  }
  function closeForm() {
    draft = null; draftPhoto = ''; $('hzform').hidden = true; $('hzform').reset(); $('hz-prev').hidden = true;
    if (clickMode === 'haz') clickMode = 'B';
  }
  function startAdd(here) {
    clickMode = 'haz'; st.msg = 'Tap the map where the hazard is.';
    if (here) openForm(st.pos);
    update();
    if (here) $('hzform').scrollIntoView({ block: 'nearest', behavior: reduce ? 'auto' : 'smooth' });
  }
  $('hz-add').addEventListener('click', function () { if (clickMode === 'haz') { closeForm(); st.msg = ''; update(); } else startAdd(false); });
  $('hz-here').addEventListener('click', function () { openForm(st.pos); update(); });
  $('hz-cancel').addEventListener('click', function () { closeForm(); st.msg = ''; update(); });
  $('hz-r').addEventListener('input', function () { set('hz-r-out', +this.value * CELL + ' m'); update(); });
  $('hz-level').addEventListener('change', function () { update(); });
  $('hz-photo').addEventListener('change', function () {
    var f = this.files && this.files[0], prev = $('hz-prev'); draftPhoto = ''; prev.hidden = true;
    if (!f || !/^image\//.test(f.type)) return;
    var rd = new FileReader();
    rd.onload = function () {
      var im = new Image();
      im.onload = function () { // shrink before storing so that reports stay small
        var k = Math.min(1, 360 / Math.max(im.width, im.height)), c = document.createElement('canvas');
        c.width = Math.max(1, Math.round(im.width * k)); c.height = Math.max(1, Math.round(im.height * k));
        c.getContext('2d').drawImage(im, 0, 0, c.width, c.height);
        draftPhoto = c.toDataURL('image/jpeg', 0.72); prev.src = draftPhoto; prev.hidden = false;
      };
      im.src = rd.result;
    };
    rd.readAsDataURL(f);
  });
  $('hzform').addEventListener('submit', function (ev) {
    ev.preventDefault();
    if (!draft) { st.msg = 'Tap the map to mark where the hazard is.'; update(); return; }
    if (st.walking) { st.msg = 'Wait until the walk is finished, then save.'; update(); return; }
    var z = cleanZone({ id: uid(), x: draft.x, y: draft.y, r: +$('hz-r').value, type: $('hz-type').value, level: +$('hz-level').value,
      desc: $('hz-desc').value.trim(), by: $('hz-by').value.trim() || 'EV-1', t: Date.now(), photo: draftPhoto, local: !!DB_URL });
    closeForm();
    st.msg = 'Hazard zone saved: ' + HZ_TYPES[z.type][0] + ', ' + LEVELS[z.level].toLowerCase() + '. Routes now avoid it.';
    setZones(zones.concat([z])); pushZone(z); ping(z.y * W + z.x);
  });

  /* ---------- map pointer ---------- */
  function cellAt(ev) {
    var r = cv.getBoundingClientRect(), x = Math.floor((ev.clientX - r.left) / r.width * W), y = Math.floor((ev.clientY - r.top) / r.height * H);
    return (x < 0 || y < 0 || x >= W || y >= H) ? -1 : y * W + x;
  }
  cv.addEventListener('click', function (ev) {
    var c = cellAt(ev); if (c < 0) return;
    if (clickMode === 'haz') { openForm(c); st.msg = ''; update(); return; }
    if (st.walking) return;
    if (clickMode === 'A') { setStart(c); update(); return; }
    var near = null, best = 6;
    SCI.forEach(function (s) { var d = Math.hypot((s.cell % W) - (c % W), ((s.cell / W) | 0) - ((c / W) | 0)); if (d < best) { best = d; near = s; } });
    pickDest(c, near);
    update();
  });
  cv.addEventListener('pointermove', function (ev) {
    var c = cellAt(ev); if (c < 0 || !out) return;
    var eo = out.dist[c], er = R().dist[c], z = zones[T.zoneAt[c]];
    set('hover', ['Safe', 'Caution', 'No entry'][T.cls[c]] + (z ? ' (' + HZ_TYPES[z.type][0] + ', ' + z.by + ')' : '') + ' · elevation ' + (T.z[c] >= 0 ? '+' : '−') + Math.abs(T.z[c]).toFixed(0) + ' m · slope ' + T.slope[c].toFixed(0) + '° · thermal inertia ' + T.I[c].toFixed(0) + ' · dose ' + (P.doseRef * T.psi[c]).toFixed(1) + ' µSv/h · O₂ there ' + (isFinite(eo) ? eo.toFixed(0) + ' L' : 'n/a') + ' · O₂ home ' + (isFinite(er) ? er.toFixed(0) + ' L' : 'n/a'));
  });

  /* ---------- voice control (Web Speech API) with a typed fallback ---------- */
  function Lx(en, ru) { return vlang === 'ru' ? ru : en; }
  function say(text) {
    set('vreply', text);
    if (!$('speak').checked || !window.speechSynthesis || !window.SpeechSynthesisUtterance) return;
    try { window.speechSynthesis.cancel(); var u = new SpeechSynthesisUtterance(text); u.lang = vlang === 'ru' ? 'ru-RU' : 'en-US'; window.speechSynthesis.speak(u); } catch (e) {}
  }
  var COMPASS = [['north', 'north-east', 'east', 'south-east', 'south', 'south-west', 'west', 'north-west'], ['север', 'северо-восток', 'восток', 'юго-восток', 'юг', 'юго-запад', 'запад', 'северо-запад']];
  function nearestHazard() {
    var px = st.pos % W, py = (st.pos / W) | 0, best = -1, bd = 1e9;
    for (var i = 0; i < N; i++) if (T.cls[i] === 2) { var d = Math.hypot((i % W) - px, ((i / W) | 0) - py); if (d < bd) { bd = d; best = i; } }
    if (best < 0) return null;
    var brg = (Math.atan2((best % W) - px, -(((best / W) | 0) - py)) * 180 / Math.PI + 360) % 360;
    return { cell: best, dist: bd * CELL, dir: Math.round(brg / 45) % 8, zone: T.uz[best] ? zones[T.zoneAt[best]] : null, rsl: !!T.rsl[best] };
  }
  function findPlace(t) {
    var m = t.match(/(?:^|[^a-zа-я])s\s?([1-5])(?![a-zа-я0-9])/) || t.match(/(?:цел[ьи]|target|объект[ау]?)\s*(?:номер\s*|number\s*)?([1-5])/);
    if (m) return SCI[+m[1] - 1];
    if (/сульфат|sulfate|sulphate/.test(t)) return SCI[0];
    if (/контакт|contact/.test(t)) return SCI[1];
    if (/блок|валун|block|boulder/.test(t)) return SCI[2];
    if (/rsl|рсл|полос|streak|standoff/.test(t)) return SCI[3];
    if (/грабен|разлом|уступ|graben|scarp|fault/.test(t)) return SCI[4];
    if (/баз|шлюз|домой|base|airlock|habitat|(?:^|\s)home/.test(t)) return 'hab';
    return null;
  }
  function summary(D) {
    var t = D.tgt, ru = vlang === 'ru', min = Math.max(1, Math.round(t.time / 60));
    var s = ru ? 'Маршрут до точки «' + destName(1) + '», ' + (mode === 'safe' ? 'самый безопасный' : 'кратчайший') + ': ' + fmtM(t.len, 1) + ', около ' + min + ' мин, примерно ' + t.steps + ' шагов, кислород ' + t.eo.toFixed(0) + ' л.'
      : 'Route to ' + destName() + ', ' + (mode === 'safe' ? 'safest' : 'shortest') + ': ' + fmtM(t.len) + ', about ' + min + ' minutes, roughly ' + t.steps + ' steps, ' + t.eo.toFixed(0) + ' litres of oxygen.';
    if (!t.ok) s += Lx(' Warning: it is beyond the point of no return.', ' Внимание: точка за пределами зоны возврата.');
    else if (t.warns.length) s += Lx(' Warnings: ' + t.warns.length + '. ', ' Предупреждений: ' + t.warns.length + '. ') + warnText(t.warns[0], ru);
    else s += Lx(' No hazards on the way.', ' Опасностей на пути нет.');
    return s;
  }
  function runCmd(raw) {
    var t = String(raw || '').toLowerCase().replace(/ё/g, 'е').trim(); if (!t) return;
    var ru = vlang === 'ru', D, hz;
    if (/(добав|отмет|сообщ|запиш).*(опасн|зон|трещин|провал)|(add|report|mark|new).*(hazard|danger|zone)/.test(t)) {
      startAdd(/здесь|тут|here|my position/.test(t));
      return say(Lx('Hazard report started. Tap the map where the hazard is, then fill in the form.', 'Добавление опасной зоны. Отметьте место на карте и заполните форму.'));
    }
    if (/(ближайш|рядом|nearest|closest).*(опасн|угроз|hazard|danger)|(где|where).*(опасн|hazard|danger)/.test(t)) {
      hz = nearestHazard(); if (!hz) return say(Lx('There are no red zones on the map.', 'Красных зон на карте нет.'));
      ping(hz.cell);
      var what = hz.zone ? Lx(HZ_TYPES[hz.zone.type][0].toLowerCase() + ', reported by ' + hz.zone.by, HZ_TYPES[hz.zone.type][1] + ', сообщил ' + hz.zone.by)
        : hz.rsl ? Lx('restricted slope with seasonal streaks', 'закрытый склон с сезонными полосами') : Lx('slope steeper than 30 degrees', 'склон круче 30 градусов');
      return say(Lx('Nearest red zone: ' + what + '. ' + fmtM(hz.dist) + ' to the ' + COMPASS[0][hz.dir] + '. It is circled on the map.',
        'Ближайшая опасная зона: ' + what + '. ' + fmtM(hz.dist, 1) + ', направление ' + COMPASS[1][hz.dir] + '. Она обведена на карте.'));
    }
    if (/(покажи|показать|show|display).*(опасн|зон|hazard|danger|safety)/.test(t)) {
      userLayer = 'safety'; update(); return say(Lx('Safety zones are on the map: green safe, yellow caution, red no entry.', 'Зоны показаны: зелёная безопасно, жёлтая осторожно, красная вход запрещён.'));
    }
    if (/сколько.*(врем|займ|идти|минут|длит|шаг|кислород|метр|километр)|как долго|как далеко|расстояни|кислород|шагов|how long|how much|how far|how many|travel time|distance|oxygen|steps/.test(t)) {
      D = derive(); if (!D.tgt) return say(Lx('There is no route yet. Choose a destination B first.', 'Маршрута пока нет. Сначала выберите точку Б.'));
      return say(summary(D));
    }
    if (/сброс|заново|reset/.test(t)) { if (!st.walking) { reset(); update(); } return say(Lx('EVA reset. You are at the airlock.', 'Выход сброшен. Вы у шлюза.')); }
    var wantSafe = /безопасн|safe/.test(t), wantFast = /кратчайш|коротк|быстр|shortest|fastest|quickest|direct/.test(t);
    var place = findPlace(t);
    if (!place && /вернись|вернуться|возвращ|return|go back|head back/.test(t)) place = 'hab';
    if (st.walking) return say(Lx('Walk in progress. Wait until it is finished.', 'Идёт движение. Дождитесь окончания.'));
    if (wantSafe) mode = 'safe'; else if (wantFast) mode = 'fast';
    if (place) {
      var ok = place === 'hab' ? pickDest(T.hab, null) : pickDest(-1, place);
      D = update();
      if (!ok || !D.tgt) return say(place === 'hab' && st.pos === T.hab ? Lx('You are already at the base.', 'Вы уже на базе.') : Lx('Route not built. ' + st.msg, 'Маршрут не построен. ' + st.msg));
      return say(summary(D));
    }
    if (wantSafe || wantFast) {
      D = update();
      var head = wantSafe ? Lx('Safest route mode is on.', 'Включён режим самого безопасного маршрута.') : Lx('Shortest route mode is on.', 'Включён режим кратчайшего маршрута.');
      return say(head + ' ' + (D.tgt ? summary(D) : Lx('Choose a destination B.', 'Выберите точку Б.')));
    }
    if (/поехали|старт|начать|начни|вперед|пошли|иди|go|start|walk|begin/.test(t)) {
      return say(goWalk() ? Lx('Walking the route.', 'Начинаю движение по маршруту.') : Lx('Cannot start: there is no allowed route. Choose a destination B.', 'Не могу начать: нет разрешённого маршрута. Выберите точку Б.'));
    }
    say(Lx('Command not recognised. Try: route to base, find the safe path, nearest danger zone, how long will it take, add hazard zone.',
      'Команда не распознана. Попробуйте: построй маршрут до базы, найди безопасный путь, где ближайшая опасная зона, сколько времени займёт дорога, добавить опасную зону.'));
  }
  var CHIPS = {
    en: ['Route to the layered sulfates', 'Find the safe path', 'Where is the nearest danger zone?', 'How long will it take?', 'Route to base', 'Add hazard zone'],
    ru: ['Построй маршрут до сульфатов', 'Найди безопасный путь', 'Где находится ближайшая опасная зона?', 'Сколько времени займёт дорога?', 'Построй маршрут до базы', 'Добавить опасную зону']
  };
  var SR = window.SpeechRecognition || window.webkitSpeechRecognition, rec = null, listening = false;
  function voiceUI() {
    press('[data-vlang]', 'data-vlang', vlang);
    var box = $('chips'); box.textContent = '';
    CHIPS[vlang].forEach(function (c) { var b = el('button', '', c); b.type = 'button'; b.addEventListener('click', function () { $('cmd').value = c; runCmd(c); }); box.appendChild(b); });
    $('cmd').placeholder = Lx('Say or type a command', 'Скажите или введите команду');
    $('mic').textContent = listening ? Lx('Listening… tap to stop', 'Слушаю… нажмите, чтобы остановить') : Lx('Speak', 'Говорить');
    $('mic').setAttribute('aria-pressed', listening ? 'true' : 'false');
    if (!SR) set('vnote', Lx('This browser has no speech recognition (Chrome, Edge and Safari do). Typed commands and the buttons below work everywhere.', 'В этом браузере нет распознавания речи (есть в Chrome, Edge и Safari). Ввод текстом и кнопки ниже работают везде.'));
  }
  $('mic').addEventListener('click', function () {
    if (!SR) { say(Lx('Voice input is not available in this browser. Type the command instead.', 'Голосовой ввод недоступен в этом браузере. Введите команду текстом.')); return; }
    if (listening && rec) { rec.stop(); return; }
    rec = new SR(); rec.lang = vlang === 'ru' ? 'ru-RU' : 'en-US'; rec.interimResults = false; rec.maxAlternatives = 1;
    rec.onresult = function (e) { var txt = e.results[0][0].transcript; $('cmd').value = txt; runCmd(txt); };
    rec.onerror = function (e) { set('vreply', Lx('Microphone problem: ', 'Проблема с микрофоном: ') + (e.error || 'unknown') + Lx('. You can type the command.', '. Команду можно ввести текстом.')); };
    rec.onend = function () { listening = false; voiceUI(); };
    try { rec.start(); listening = true; voiceUI(); } catch (err) { listening = false; voiceUI(); }
  });
  $('cmdform').addEventListener('submit', function (ev) { ev.preventDefault(); runCmd($('cmd').value); });
  Array.prototype.forEach.call(document.querySelectorAll('[data-vlang]'), function (bn) {
    bn.addEventListener('click', function () { vlang = bn.getAttribute('data-vlang'); voiceUI(); });
  });

  /* ---------- start ---------- */
  fillSelects();
  if (window.innerWidth < 700 && $('keybox')) $('keybox').open = false;
  zones = loadLocal() || SEEDS.map(cleanZone);
  applyZones(); reset(); voiceUI(); update();
  window.addEventListener('storage', function (e) { if (e.key === LS_KEY) { var l = loadLocal(); if (l) setZones(l, true); } }); // other tabs
  if (DB_URL) { pull(); setInterval(pull, 8000); }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(function () { update(); });
})();
