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
    var block = new Uint8Array(N), hazard = new Uint8Array(N), rsl = new Uint8Array(N), dune = new Float32Array(N);
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
    return { z: z, I: I, slope: slope, psi: psi, block: block, hazard: hazard, rsl: rsl, dune: dune, hab: HAB.y * W + HAB.x };
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

  function dijkstra(T, src, v, reverse) {
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
        if (T.block[n]) continue;
        var dd = (k & 1) ? CELL * 1.41421356 : CELL;
        var c = reverse ? edge(T, n, u, dd, v).o2 : edge(T, u, n, dd, v).o2;
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
  var C = AeroCore, W = C.W, H = C.H, CELL = C.CELL, P = C.P, S = 5;
  var $ = function (id) { return document.getElementById(id); };
  var cv = $('map'); if (!cv) return;
  cv.width = W * S; cv.height = H * S;
  var ctx = cv.getContext('2d');
  var T = C.build();
  var baseCv = document.createElement('canvas'); baseCv.width = W; baseCv.height = H;
  var overCv = document.createElement('canvas'); overCv.width = W; overCv.height = H;
  var bctx = baseCv.getContext('2d'), octx = overCv.getContext('2d');
  var reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  function findSpot(x, y) { // nearest walkable, gentle cell
    for (var r = 0; r < 30; r++) for (var dy = -r; dy <= r; dy++) for (var dx = -r; dx <= r; dx++) {
      if (Math.max(Math.abs(dx), Math.abs(dy)) !== r) continue;
      var px = x + dx, py = y + dy; if (px < 1 || py < 1 || px >= W - 1 || py >= H - 1) continue;
      var i = py * W + px; if (!T.block[i] && T.slope[i] < 12 && T.I[i] > 220) return i;
    }
    return y * W + x;
  }
  var SCI0 = [
    { id: 'A', name: 'Layered sulfates', note: 'Hydrated sulfate layers, a sign of past water', cls: 'sulfate', s2: 4.0, x: 131, y: 96 },
    { id: 'B', name: 'Layer contact', note: 'Where layered deposits meet older rock', cls: 'sulfate', s2: 3.0, x: 171, y: 84 },
    { id: 'C', name: 'Block field', note: 'Bright 100–500 m blocks in dark material', cls: 'bedrock', s2: 2.0, x: 196, y: 122 },
    { id: 'D', name: 'RSL standoff', note: 'Viewpoint below seasonal slope streaks', cls: 'rsl', s2: 5.0, x: 108, y: 62 },
    { id: 'E', name: 'Graben scarp', note: 'Fault trough cutting the canyon floor', cls: 'tectonic', s2: 1.5, x: 104, y: 136 }
  ];
  var SCI = [];
  var st = {};
  var ret = { nom: C.dijkstra(T, T.hab, P.vNom, true), shed: C.dijkstra(T, T.hab, P.vShed, true) };
  var out = null, outKey = '';
  var userLayer = 'relief';

  function reset() {
    SCI = SCI0.map(function (s) { return { id: s.id, name: s.name, note: s.note, cls: s.cls, s2: s.s2, cell: findSpot(s.x, s.y), done: false }; });
    st = { V: 1500 / P.psiPerL, F: 0.2, pos: T.hab, target: -1, sci: null, elapsed: 0, dose: 0, walked: 0, walking: false, msg: '' };
    st.ema = C.flow(C.metW(T, T.hab, T.hab + 1, CELL, P.vNom));
    st.W = 0;
    $('tank').value = 1500; $('fat').value = 20;
    outKey = '';
  }
  function shedMode() { return st.F > 0.55; }
  function speed() { return shedMode() ? P.vShed : P.vNom; }
  function R() { return shedMode() ? ret.shed : ret.nom; }
  function ensureOut() {
    var key = st.pos + ':' + (shedMode() ? 1 : 0);
    if (key !== outKey) { out = C.dijkstra(T, st.pos, speed(), false); outKey = key; }
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
  var LEG = {
    relief: ['Synthetic shaded relief', 'low', 'high', RAMPS.relief],
    slope: ['Slope, degrees', '0°', '25°+', [[70, 150, 96], [222, 196, 84], [226, 128, 52], [200, 52, 44]]],
    thermal: ['Thermal inertia, TIU', '150 loose sand', '452 bedrock', RAMPS.thermal],
    rad: ['Dose rate, µSv/h', '12 sheltered', '26.7 open sky', RAMPS.rad],
    o2: ['O₂ to walk home, litres', '0', '200+', RAMPS.o2]
  };
  function paintBase(layer) {
    var img = bctx.createImageData(W, H), d = img.data, Rr = R().dist;
    for (var y = 0; y < H; y++) for (var x = 0; x < W; x++) {
      var i = y * W + x, xl = Math.max(0, x - 1), xr = Math.min(W - 1, x + 1), yu = Math.max(0, y - 1), yd = Math.min(H - 1, y + 1);
      var gx = (T.z[y * W + xr] - T.z[y * W + xl]) / ((xr - xl) * CELL), gy = (T.z[yd * W + x] - T.z[yu * W + x]) / ((yd - yu) * CELL);
      var sh = (gx * 0.62 - gy * 0.62 + 0.48) / Math.sqrt(gx * gx + gy * gy + 1) / 0.48; // light from south-west, so the south-facing wall is lit
      sh = Math.min(1.3, Math.max(0.5, sh));
      var c, k = 1;
      if (layer === 'relief') {
        var zz = T.z[i]; c = ramp(RAMPS.relief, zz < 320 ? (zz + 90) / 410 * 0.66 : 0.66 + Math.min(1, (zz - 320) / 1700) * 0.34);
        var dn = T.dune[i]; c = [c[0] + (62 - c[0]) * dn * 0.7, c[1] + (52 - c[1]) * dn * 0.7, c[2] + (56 - c[2]) * dn * 0.7];
        k = sh;
      } else if (layer === 'slope') { c = ramp(LEG.slope[3], T.slope[i] / 25); k = 0.55 + 0.45 * sh; }
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
  function straightHaz(a, b) {
    var x0 = a % W, y0 = (a / W) | 0, x1 = b % W, y1 = (b / W) | 0, n = Math.max(Math.abs(x1 - x0), Math.abs(y1 - y0)), c = 0;
    for (var i = 1; i < n; i++) { var x = Math.round(x0 + (x1 - x0) * i / n), y = Math.round(y0 + (y1 - y0) * i / n), q = y * W + x; if (T.hazard[q] || T.block[q]) c++; }
    return { haz: c, dist: Math.hypot(x1 - x0, y1 - y0) * CELL };
  }
  function pathLen(cells) { var L = 0; for (var i = 1; i < cells.length; i++) { var a = cells[i - 1], b = cells[i]; var dd = (a % W !== b % W && ((a / W) | 0) !== ((b / W) | 0)) ? CELL * 1.41421356 : CELL; L += Math.hypot(dd, T.z[b] - T.z[a]); } return L; }

  function derive() {
    ensureOut();
    var Rf = R(), b = C.beta(st.F), p = st.V * P.psiPerL, budget = b * (st.V - P.reserve);
    var retHere = Rf.dist[st.pos], margin = budget - retHere;
    var breach = (st.V - retHere <= P.reserve) || p < 650;
    var D = { beta: b, p: p, budget: budget, retHere: retHere, margin: margin, breach: breach, shed: shedMode(), tgt: null };
    if (st.target >= 0 && isFinite(out.dist[st.target])) {
      var eo = out.dist[st.target], er = Rf.dist[st.target], route = C.chain(out.prev, st.target, st.pos).reverse();
      var sh = straightHaz(st.pos, st.target), len = pathLen(route);
      D.tgt = { eo: eo, er: er, after: b * (st.V - eo - P.reserve) - er, route: route, back: C.chain(Rf.prev, st.target, T.hab), len: len, haz: sh.haz, straight: sh.dist };
      D.tgt.ok = D.tgt.after >= 0;
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
    SAND_TRAP_WARN: ['warn', 'The straight course crosses loose sand or slopes over 18°. Route redrawn over firmer ground.'],
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
  function label(txt, x, y, color, align) {
    ctx.save(); ctx.font = '600 22px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = align || 'left'; ctx.textBaseline = 'middle';
    ctx.lineWidth = 5; ctx.strokeStyle = 'rgba(20,10,6,.85)'; ctx.strokeText(txt, x, y); ctx.fillStyle = color; ctx.fillText(txt, x, y); ctx.restore();
  }
  function draw(D) {
    var forced = D.shed || D.breach, layer = forced ? 'relief' : userLayer;
    var key = layer + (layer === 'o2' ? (D.shed ? 1 : 0) : '');
    if (key !== paintedKey) { paintBase(layer); paintedKey = key; }
    ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(baseCv, 0, 0, W * S, H * S);
    // overlay: outside-contour shading + RSL zones
    var Rf = R().dist, img = octx.createImageData(W, H), d = img.data, bud = D.budget, b = D.beta, g = new Float32Array(W * H);
    for (var i = 0; i < W * H; i++) {
      var f = Rf[i] + b * out.dist[i] - bud; if (!isFinite(f)) f = 1e3; g[i] = f;
      var o = i * 4;
      if (T.rsl[i]) { var hx = ((i % W) + ((i / W) | 0)) % 3 === 0; d[o] = 105; d[o + 1] = 200; d[o + 2] = 232; d[o + 3] = hx ? 200 : 70; }
      else if (f > 0) { d[o] = 30; d[o + 1] = 6; d[o + 2] = 4; d[o + 3] = 150; }
    }
    octx.putImageData(img, 0, 0);
    ctx.drawImage(overCv, 0, 0, W * S, H * S);
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
      strokePath(D.tgt.route, D.tgt.ok ? '#ffffff' : '#ef5a45', 3.5);
    }
    if (st.pos !== T.hab) strokePath(D.egress, '#f0b24a', forced ? 6 : 3.5);
    // science targets
    if (!forced) SCI.forEach(function (s) {
      var X = cx(s.cell), Y = cy(s.cell);
      ctx.save(); ctx.translate(X, Y); ctx.rotate(Math.PI / 4);
      ctx.fillStyle = s.done ? 'rgba(60,48,42,.9)' : '#1c1410'; ctx.strokeStyle = s.done ? '#8a7667' : '#f3e6d6'; ctx.lineWidth = 2.5;
      ctx.fillRect(-15, -15, 30, 30); ctx.strokeRect(-15, -15, 30, 30); ctx.restore();
      ctx.save(); ctx.font = '700 22px "IBM Plex Mono", ui-monospace, monospace'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
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
    // map furniture: scale bar and north arrow
    var bx = 28, by = H * S - 30, km = 1000 / CELL * S;
    ctx.save(); ctx.strokeStyle = '#f3e6d6'; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(bx, by - 7); ctx.lineTo(bx, by); ctx.lineTo(bx + km, by); ctx.lineTo(bx + km, by - 7); ctx.stroke(); ctx.restore();
    label('1 km', bx + km + 10, by - 3, '#f3e6d6');
    label('N ↑', W * S - 28, 32, '#f3e6d6', 'right');
    label('CANYON WALL', W * S / 2, 32, 'rgba(243,230,214,.8)', 'center');
  }

  /* ---------- panel ---------- */
  function set(id, v) { var e = $(id); if (e && e.textContent !== v) e.textContent = v; }
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
    // target card
    var tc = $('tgt');
    if (D.tgt && !forced) {
      var t = D.tgt, nm = st.sci ? st.sci.id + ' · ' + st.sci.name : 'Waypoint';
      set('t-name', nm); set('t-note', st.sci ? st.sci.note + (st.sci.done ? ' (sampled)' : '') : 'Free waypoint picked on the map');
      set('t-len', (t.len / 1000).toFixed(2) + ' km'); set('t-out', t.eo.toFixed(0) + ' L'); set('t-back', t.er.toFixed(0) + ' L');
      set('t-after', (t.after >= 0 ? '+' : '−') + Math.abs(t.after).toFixed(0) + ' L');
      $('t-after').className = 'num ' + (t.after > 80 ? 'c-ok' : t.after >= 0 ? 'c-warn' : 'c-bad');
      var detour = t.len - t.straight;
      set('t-verdict', !t.ok ? 'Beyond the point of no return. Walking there is blocked.' :
        (t.haz > 0 ? 'Straight line crosses ' + t.haz + ' hazard cells. Detour adds ' + Math.max(0, detour).toFixed(0) + ' m and keeps to firmer ground.' : 'Inside the survival contour. Clear to go.'));
      $('t-verdict').className = 'verdict ' + (t.ok ? (t.haz ? 'c-warn' : 'c-ok') : 'c-bad');
      tc.hidden = false; $('t-empty').hidden = true;
    } else {
      tc.hidden = true; $('t-empty').hidden = false;
      set('t-empty', forced ? (D.breach ? 'Science targets are locked. Only the return route is available.' : 'Science targets are hidden while fatigue is above 0.55.') : 'Click the map or pick a target below to plan a route.');
    }
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
    // telemetry packet
    var px = st.pos % W, py = (st.pos / W) | 0;
    var pkt = {
      packet_id: 'DEMO-EVA-' + ('00000' + Math.round(st.elapsed)).slice(-5), crew_id: 'EV-1',
      location: { frame: 'local grid, metres from airlock (synthetic terrain)', east_m: (px - C.HAB.x) * CELL, north_m: (C.HAB.y - py) * CELL, elevation_rel_m: +T.z[st.pos].toFixed(1), local_slope_deg: +T.slope[st.pos].toFixed(1), soil_thermal_inertia_tiu: +T.I[st.pos].toFixed(0), sky_view_factor: +T.psi[st.pos].toFixed(2) },
      plss_life_support: { o2_primary_pressure_psia: +D.p.toFixed(0), o2_primary_usable_liters_stp: +st.V.toFixed(1), o2_flow_rate_lpm_ema180: +st.ema.toFixed(2), time_remaining_o2_sec: Math.max(0, Math.round((st.V - P.reserve) / st.ema * 60)) },
      biometrics_and_risk: { metabolic_power_watts: +st.W.toFixed(0), fatigue_index_f: +st.F.toFixed(2), present_bias_beta: +D.beta.toFixed(3), radiation_dose_rate_usv_hr: +(P.doseRef * T.psi[st.pos]).toFixed(1), cumulative_eva_dose_msv: +(st.dose / 1000).toFixed(4) },
      navigation_pnr: { distance_to_airlock_m: +home.toFixed(0), optimal_return_path_o2_liters: +D.retHere.toFixed(1), pnr_margin_liters: +D.margin.toFixed(1), hud_state: D.code, hud_cognitive_mode: forced ? 'EGRESS_ONLY' : 'STANDARD_5_LAYER' }
    };
    set('packet', JSON.stringify(pkt, null, 2));
  }
  function update() { var D = derive(); draw(D); panel(D); return D; }

  /* ---------- walking ---------- */
  function walk(cells, onDone) {
    if (cells.length < 2) return;
    st.walking = true; st.msg = '';
    var v = speed(), i = 1, per = reduce ? cells.length : Math.max(1, Math.ceil(cells.length / 90));
    function step() {
      for (var n = 0; n < per && i < cells.length; n++, i++) {
        var a = cells[i - 1], b = cells[i], dd = (a % W !== b % W && ((a / W) | 0) !== ((b / W) | 0)) ? CELL * 1.41421356 : CELL;
        var e = C.edge(T, a, b, dd, v), al = Math.min(1, e.t / 180), fl = C.flow(e.W);
        st.V = Math.max(0, st.V - e.o2); st.elapsed += e.t; st.walked += e.len; st.W = e.W;
        st.ema = al * fl + (1 - al) * st.ema; st.dose += P.doseRef * T.psi[b] * e.t / 3600; st.pos = b;
      }
      $('tank').value = Math.round(st.V * P.psiPerL);
      if (i < cells.length) { update(); requestAnimationFrame(step); }
      else { st.walking = false; onDone(); update(); }
    }
    step();
  }
  $('go').addEventListener('click', function () {
    var D = derive(); if (!D.tgt || !D.tgt.ok) return;
    var sci = st.sci;
    walk(D.tgt.route, function () {
      if (sci && !sci.done) {
        sci.done = true; var n = 0;
        SCI.forEach(function (o) { if (o !== sci && o.cls === sci.cls && !o.done) { o.s2 *= 0.25; n++; } });
        st.msg = 'Sample logged at ' + sci.name + '.' + (n ? ' Uncertainty at ' + n + ' similar target' + (n > 1 ? 's' : '') + ' collapsed, so the ranking changed.' : '');
      } else st.msg = 'Arrived at the waypoint.';
      st.target = -1; st.sci = null;
    });
  });
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
  function cellAt(ev) {
    var r = cv.getBoundingClientRect(), x = Math.floor((ev.clientX - r.left) / r.width * W), y = Math.floor((ev.clientY - r.top) / r.height * H);
    return (x < 0 || y < 0 || x >= W || y >= H) ? -1 : y * W + x;
  }
  cv.addEventListener('click', function (ev) {
    if (st.walking) return; var c = cellAt(ev); if (c < 0) return;
    var D = derive(); if (D.shed || D.breach) { st.msg = 'Targets are locked in this HUD state. Lower fatigue or use Return to airlock.'; update(); return; }
    var near = null, best = 6;
    SCI.forEach(function (s) { var d = Math.hypot((s.cell % W) - (c % W), ((s.cell / W) | 0) - ((c / W) | 0)); if (d < best) { best = d; near = s; } });
    if (near) { st.target = near.cell; st.sci = near; st.msg = ''; }
    else if (T.rsl[c]) { st.msg = 'Restricted: slope with seasonal streaks (planetary protection).'; }
    else if (T.block[c]) { st.msg = 'Too steep to walk: slope over 30°.'; }
    else if (!isFinite(out.dist[c])) { st.msg = 'No walkable route to that point.'; }
    else { st.target = c; st.sci = null; st.msg = ''; }
    update();
  });
  cv.addEventListener('pointermove', function (ev) {
    var c = cellAt(ev); if (c < 0 || !out) return;
    var eo = out.dist[c], er = R().dist[c];
    set('hover', 'Elevation ' + (T.z[c] >= 0 ? '+' : '−') + Math.abs(T.z[c]).toFixed(0) + ' m · slope ' + T.slope[c].toFixed(0) + '° · thermal inertia ' + T.I[c].toFixed(0) + ' · dose ' + (P.doseRef * T.psi[c]).toFixed(1) + ' µSv/h · O₂ there ' + (isFinite(eo) ? eo.toFixed(0) + ' L' : 'n/a') + ' · O₂ home ' + (isFinite(er) ? er.toFixed(0) + ' L' : 'n/a'));
  });

  reset(); update();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(function () { update(); });
})();
