// Model tests. Run with `npm test` (Node 20+, no dependencies).
import test from 'node:test';
import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
const { AeroCore: C } = createRequire(import.meta.url)('../app.js');

test('oxygen flow rises with metabolic rate and stays in a plausible suit range', () => {
  let prev = 0;
  for (let met = 1; met <= 8; met += 0.5) {
    const f = C.flowFromMet(met);
    assert.ok(f > prev, `flow must rise with MET (at ${met})`);
    assert.ok(f > 0.25 && f < 3, `flow ${f.toFixed(2)} L/min at ${met} MET is outside 0.25-3`);
    prev = f;
  }
  // From 2 MET up the flow is inside 0.5-3 L/min. At rest (1 MET, 93 W for 80 kg) physics gives about 0.31 L/min.
  for (let met = 2; met <= 8; met++) assert.ok(C.flowFromMet(met) >= 0.5 && C.flowFromMet(met) <= 3);
  assert.ok(Math.abs(C.flowFromMet(1) - 0.311) < 0.01);
});

test('MET and watts convert both ways', () => {
  assert.ok(Math.abs(C.met(C.P.wPerMetKg * C.P.mBody) - 1) < 1e-9);
  assert.ok(Math.abs(C.met(C.P.Wbase) - 1.236) < 0.01, 'standing in the suit is about 1.2 MET');
  assert.ok(Math.abs(C.flow(280) - 0.866) < 0.002, '280 W gives 0.866 L/min, as stated on the site');
});

test('RAD dose rate: mSv/day and microSv/h agree', () => {
  assert.ok(Math.abs(C.mSvDayToUSvHour(0.64) - 26.667) < 0.001);
  assert.ok(Math.abs(C.uSvHourToMSvDay(C.mSvDayToUSvHour(0.64)) - 0.64) < 1e-12);
  assert.ok(Math.abs(C.P.doseRef - C.mSvDayToUSvHour(C.P.doseDaily)) < 1e-12, 'the model uses the same figure in both units');
});

test('walking uphill costs more than flat, and soft ground costs more than firm ground', () => {
  const flat = C.metWraw(0, 452, 0.8), up = C.metWraw(Math.sin(10 * Math.PI / 180), 452, 0.8), soft = C.metWraw(0, 180, 0.8);
  assert.ok(up > flat && soft > flat);
  assert.ok(C.metWraw(-0.05, 452, 0.8) >= C.P.Wbase, 'power never drops below standing');
});
