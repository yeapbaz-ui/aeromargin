// Checks data/manifest.json: every dataset has the required fields, and every tile on disk matches its entry.
// Run with: node tools/validate_manifest.mjs   (also runs inside `npm test`)
import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');

export function validateManifest(man, readFile = (f) => readFileSync(join(root, 'data', f))) {
  const errors = [];
  const need = (obj, keys, where) => keys.forEach((k) => { if (obj[k] === undefined || obj[k] === '') errors.push(`${where}: missing "${k}"`); });
  need(man, ['region', 'tiles', 'datasets', 'last_checked'], 'manifest');
  need(man.region || {}, ['name', 'center_lat', 'center_lon_east', 'width_m', 'height_m', 'cell_m'], 'region');
  const ids = new Set();
  for (const d of man.datasets || []) {
    const where = `dataset ${d.id || '?'}`;
    need(d, ['id', 'role', 'name', 'provider', 'url', 'units', 'license', 'last_checked'], where);
    if (!('product_id' in d)) errors.push(`${where}: missing "product_id"`);
    if (d.product_id === null && !d.verify) errors.push(`${where}: product_id is null and there is no "verify" note`);
    if (d.url && !/^https:\/\//.test(d.url)) errors.push(`${where}: url must be https`);
    if (d.last_checked && !/^\d{4}-\d{2}-\d{2}$/.test(d.last_checked)) errors.push(`${where}: last_checked must be YYYY-MM-DD`);
    if (ids.has(d.id)) errors.push(`${where}: duplicate id`);
    ids.add(d.id);
  }
  for (const t of man.tiles || []) {
    const where = `tile ${t.id || '?'}`;
    need(t, ['id', 'layer', 'dataset', 'file', 'width', 'height', 'bytes', 'resolution_m', 'margin_px', 'valid_min', 'valid_max', 'sha256'], where);
    if (!ids.has(t.dataset)) errors.push(`${where}: dataset "${t.dataset}" is not listed`);
    let buf;
    try { buf = readFile(t.file); } catch (e) { errors.push(`${where}: cannot read ${t.file}`); continue; }
    if (buf.length !== t.bytes || buf.length !== t.width * t.height * 4) errors.push(`${where}: file has ${buf.length} bytes, expected ${t.width * t.height * 4}`);
    const sha = createHash('sha256').update(buf).digest('hex');
    if (sha !== t.sha256) errors.push(`${where}: sha256 is ${sha}, manifest says ${t.sha256}`);
  }
  return errors;
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const errors = validateManifest(JSON.parse(readFileSync(join(root, 'data', 'manifest.json'), 'utf8')));
  if (errors.length) { console.error('manifest.json is not valid:\n - ' + errors.join('\n - ')); process.exit(1); }
  console.log('manifest.json is valid');
}
