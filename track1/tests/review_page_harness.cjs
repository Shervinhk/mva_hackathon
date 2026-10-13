/* Execute the generated review-page JavaScript with fabricated DOM inputs.
 * This checks selection and exported bytes, not browser layout or rendering.
 */
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');

const html = fs.readFileSync(process.argv[2], 'utf8');
const payload = html.match(/<script id="payload" type="application\/json">([\s\S]*?)<\/script>/)[1];
const code = html.match(/<\/script><script>([\s\S]*?)<\/script>/)[1];
const downloads = [];
const blobs = new Map();
class Element {
  constructor(tag) { this.tag = tag; this.children = []; this.value = ''; this.checked = false; this.textContent = ''; }
  append(...items) {
    this.children.push(...items);
    if (this.tag === 'select' && !this.value && items[0]) this.value = items[0].value;
  }
  setAttribute() {}
  remove() {}
  click() { if (this.download) downloads.push({name: this.download, text: blobs.get(this.href)}); }
}
const ids = Object.fromEntries(['payload', 'candidates', 'reviewed', 'manifest', 'csv', 'message'].map(id => [id, new Element('div')]));
ids.payload.textContent = payload;
const context = vm.createContext({
  document: {getElementById: id => ids[id], createElement: tag => new Element(tag),
    createTextNode: text => ({textContent: text}), body: new Element('body')},
  structuredClone, setTimeout: fn => fn(),
  Blob: class { constructor(parts) { this.text = parts.join(''); } },
  URL: {createObjectURL: blob => { const id = 'blob:' + blobs.size; blobs.set(id, blob.text); return id; }, revokeObjectURL() {}},
});
vm.runInContext(code, context);
ids.csv.onclick();
assert.equal(downloads.length, 0, 'Unreviewed export must not produce a file');
assert.match(ids.message.textContent, /Confirm local review/);
ids.reviewed.checked = true;
ids.csv.onclick();
assert.equal(downloads.length, 0, 'Empty selection must not produce a file');
vm.runInContext(`elements[0].include.checked=true; elements[0].reviewed.checked=true; elements[0].reason.value='Fabricated evidence, with a comma and a "quotation".';`, context);
ids.csv.onclick();
assert.match(ids.message.textContent, /Supply EPCR/);
assert.equal(downloads.length, 0, 'Missing confidence must not be replaced by a ranking score');
vm.runInContext("elements[0].confidence.value='0.42'", context);
ids.manifest.onclick();
ids.csv.onclick();
assert.equal(downloads.length, 2);
const selection = JSON.parse(downloads[0].text);
assert.equal(selection.reviewed_locally, true);
assert.equal(selection.choices[0].epcr, .42);
assert.equal(selection.choices[0].evidence_reviewed, true);
assert.equal(downloads[1].name, 'predictions.csv');
assert.ok(downloads[1].text.includes('\r\n'), 'CSV requires actual row separators');
assert.ok(!downloads[1].text.includes('\\r\\n'), 'CSV must not contain escaped row-separator text');
process.stdout.write(JSON.stringify({selection, csv: downloads[1].text}));
