import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { validateManifest } from '../tools/validate_manifest.mjs';
const { AeroCore: C, AeroData: D } = createRequire(import.meta.url)('../app.js');
const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const man = JSON.parse(readFileSync(join(root, 'data', 'manifest.json'), 'utf8'));
const tile = man.tiles.find((t) => t.layer === 'elevation');
const file = readFileSync(join(root, 'data', tile.file));
const ab = (b) => b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);

test('manifest.json is complete and matches the files on disk', () => {
  assert.deepEqual(validateManifest(man), []);
});

test('the loader accepts the tile that the manifest describes', async () => {
  const a = D.validateTile(ab(file), tile);
  assert.equal(a.length, tile.width * tile.height);
  assert.equal(await D.sha256(ab(file)), tile.sha256);
});

test('the loader rejects tiles that do not match the manifest', () => {
  assert.throws(() => D.validateTile(ab(file.subarray(0, file.length - 4)), tile), /bytes/);
  assert.throws(() => D.validateTile(ab(file), { ...tile, width: tile.width + 1 }), /bytes/);
  const nan = Buffer.from(file); nan.writeFloatLE(NaN, 40);
  assert.throws(() => D.validateTile(ab(nan), tile), /NaN/);
  const high = Buffer.from(file); high.writeFloatLE(50000, 40);
  assert.throws(() => D.validateTile(ab(high), tile), /outside the manifest range/);
  assert.throws(() => D.validateTile(ab(file), { id: 'x' }), /incomplete/);
});

test('a changed file fails the manifest checksum', () => {
  const errs = validateManifest(man, () => { const b = Buffer.from(file); b.writeFloatLE(-4000, 0); return b; });
  assert.ok(errs.some((e) => /sha256/.test(e)));
});

test('the walking grid is built from the tile: real elevations, an open airlock, a reachable floor', () => {
  const T = C.build(D.raster(D.validateTile(ab(file), tile), tile), null);
  assert.ok(T.zLo >= tile.min - 60 && T.zHi <= tile.max + 60, 'interpolated elevations stay close to the tile range');
  assert.ok(T.slope[T.hab] < 5 && !T.block[T.hab], 'the airlock stands on gentle open ground');
  const home = C.dijkstra(T, T.hab, C.P.vNom, true);
  let reach = 0; for (let i = 0; i < C.N; i++) if (isFinite(home.dist[i])) reach++;
  assert.ok(reach / C.N > 0.6, 'most of the map can be walked back from');
  assert.equal(T.hasSoil, false, 'no soil tile is bundled, and the grid says so');
});
