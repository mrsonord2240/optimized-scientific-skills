// Executes selected, inspected, offline demos from isolated pinned packages.
// This is bounded smoke coverage, not a scientific audit of every workflow.
import assert from 'node:assert/strict';
import { execFileSync, spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import os from 'node:os';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { parseArgs } from 'node:util';

const { values } = parseArgs({ options: {
  python: { type: 'string', default: 'python' },
  rscript: { type: 'string', default: 'Rscript' },
  only: { type: 'string', multiple: true },
  'export-only': { type: 'boolean', default: false },
  'scratch-parent': { type: 'string', default: os.tmpdir() },
} });
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const evidence = JSON.parse(readFileSync(path.join(root, 'authoring/submissions/mrsonord2240/pilot-review-evidence.json')));
const commit = evidence.submitted_source.commit;
assert.match(commit, /^[0-9a-f]{40}$/);
const scratch = mkdtempSync(path.join(path.resolve(values['scratch-parent']), 'scientific-pilot-smoke-'));
const git = (...args) => execFileSync('git', ['-C', root, ...args]);
const cases = [
  ['bio-alignment-msa-statistics', 'examples/selftest.py', ['All checks passed.']],
  ['bio-alignment-pairwise', 'examples/empirical_pvalue.py', ['Empirical p-value:']],
  ['bio-proteomics-data-import', 'examples/load_maxquant.py', ['read=14, after bookkeeping=11, quantified=11', 'MNAR']],
  ['bio-proteomics-peptide-identification', 'examples/fdr_filtering.py', ['Target PSMs at q <= 0.01:']],
  ['bio-proteomics-protein-inference', 'examples/protein_groups.py', ['Parsimony kept 7 groups out of 8', 'Protein groups passing 1% picked-group FDR: 2']],
  ['bio-proteomics-spectral-libraries', 'examples/build_library.py', ['R^2 =', 'both retained: True']],
  ['bio-proteomics-quantification', 'examples/lfq_normalization.py', ['interactors 5/5, sticky binders called 0/10']],
];
const selected = new Set(values.only ?? [...cases.map(test => test[0]),
  'bio-single-cell-preprocessing', 'bio-single-cell-clustering']);
const runnable = new Set([...cases.map(test => test[0]),
  'bio-single-cell-preprocessing', 'bio-single-cell-clustering']);
if (!values['export-only']) {
  for (const id of selected)
    assert(runnable.has(id), `No smoke test defined for ${id}; use --export-only for snapshot export`);
}
const exportManifest = { source_commit: commit, files: {} };
for (const id of selected) {
  assert(evidence.skills.some(skill => skill.id === id), `Unknown pilot Skill: ${id}`);
  const prefix = `skills/${id}/`;
  const files = git('ls-tree', '-r', '--name-only', commit, '--', `skills/${id}`).toString().trim().split('\n');
  for (const file of files) {
    const output = path.join(scratch, id, file.slice(prefix.length));
    mkdirSync(path.dirname(output), { recursive: true });
    const bytes = git('show', `${commit}:${file}`);
    writeFileSync(output, bytes);
    exportManifest.files[`${id}/${file.slice(prefix.length)}`] = createHash('sha256').update(bytes).digest('hex');
  }
}
writeFileSync(path.join(scratch, 'source-manifest.json'), JSON.stringify(exportManifest, null, 2) + '\n');
if (values['export-only']) {
  process.stdout.write(JSON.stringify({ source_commit: commit, scratch_directory: scratch }) + '\n');
  process.exit(0);
}
const results = [];
function run(id, script, command, args, expected) {
  if (values.only && !values.only.includes(id)) return;
  const directory = path.join(scratch, id);
  const result = spawnSync(command, args, {
    cwd: directory, encoding: 'utf8', timeout: 90000,
    env: { ...process.env, PYTHONDONTWRITEBYTECODE: '1', PYTHONPATH: '', MPLBACKEND: 'Agg', LC_ALL: 'C', LANG: 'C' },
    maxBuffer: 4 * 1024 * 1024,
  });
  const stdout = result.stdout ?? '';
  const stderr = result.stderr ?? '';
  const checks = expected.map(text => ({ expected: text, passed: stdout.includes(text) }));
  const passed = result.status === 0 && checks.every(check => check.passed);
  const ordinal = results.length + 1;
  writeFileSync(path.join(directory, `smoke-${ordinal}.stdout.log`), stdout);
  writeFileSync(path.join(directory, `smoke-${ordinal}.stderr.log`), stderr);
  results.push({ id, script, command, args, status: passed ? 'passed' : 'failed',
    exit_code: result.status, error: result.error?.message ?? null, checks,
    source_sha256: createHash('sha256').update(git('show', `${commit}:skills/${id}/${script}`)).digest('hex'),
    stdout, stderr });
  process.stderr.write(`${passed ? 'PASS' : 'FAIL'} ${id}/${script}\n`);
}
for (const [id, script, expected] of cases)
  run(id, script, values.python, ['-B', script], expected);

// Deterministic synthetic counts: no patient data, downloads, or model services.
const scanpyFixture = `
import gzip, pathlib, runpy
import numpy as np
from scipy import io, sparse
import scanpy as sc
rng = np.random.default_rng(29)
counts = rng.poisson(rng.lognormal(0, 0.15, (180, 1)) * rng.gamma(0.5, 1.2, (1, 3000)))
folder = pathlib.Path('filtered_feature_bc_matrix')
folder.mkdir()
with gzip.open(folder / 'matrix.mtx.gz', 'wb') as f: io.mmwrite(f, sparse.coo_matrix(counts.T))
with gzip.open(folder / 'barcodes.tsv.gz', 'wt') as f: f.write(''.join(f'cell{i}\\n' for i in range(180)))
with gzip.open(folder / 'features.tsv.gz', 'wt') as f:
    f.write(''.join(f'gene{i}\\t' + ('MT-' if i < 30 else '') + f'G{i}\\tGene Expression\\n' for i in range(3000)))
runpy.run_path('examples/preprocess_scanpy.py', run_name='__main__')
result = sc.read_h5ad('preprocessed.h5ad')
assert result.n_obs >= 144 and result.n_vars > 50
assert 'counts' in result.layers and 'X_pca' in result.obsm
assert not np.allclose(result.X.data, result.layers['counts'].data)
print('SYNTHETIC_PREPROCESS_ASSERTIONS_OK')
`;
run('bio-single-cell-preprocessing', 'examples/preprocess_scanpy.py', values.python,
  ['-B', '-c', scanpyFixture], ['SYNTHETIC_PREPROCESS_ASSERTIONS_OK']);
const clusterFixture = `
import runpy
import numpy as np
import scanpy as sc
from scipy import sparse
rng = np.random.default_rng(41)
adata = sc.AnnData(sparse.csr_matrix(rng.poisson(1.0, (180, 120)).astype(float)))
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad('preprocessed.h5ad')
runpy.run_path('examples/cluster_scanpy.py', run_name='__main__')
result = sc.read_h5ad('clustered.h5ad')
assert result.n_obs == 180 and result.n_vars == 120
assert 'leiden' in result.obs and result.obsm['X_umap'].shape == (180, 2)
assert np.isfinite(result.obsm['X_umap']).all()
print('SYNTHETIC_CLUSTER_ASSERTIONS_OK')
`;
run('bio-single-cell-clustering', 'examples/cluster_scanpy.py', values.python,
  ['-B', '-c', clusterFixture], ['SYNTHETIC_CLUSTER_ASSERTIONS_OK']);

if (selected.has('bio-proteomics-quantification')) {
const quant = path.join(scratch, 'bio-proteomics-quantification');
const rows = ['protein,ion,run,intensity'];
for (let protein = 1; protein <= 4; protein++)
  for (let ion = 1; ion <= 3; ion++)
    for (let sample = 1; sample <= 3; sample++)
      rows.push(`P${protein},I${ion},S${sample},${protein * ion * sample * 10000}`);
writeFileSync(path.join(quant, 'smoke-peptides.csv'), rows.join('\n') + '\n');
run('bio-proteomics-quantification', 'scripts/maxlfq_iq.R', values.rscript,
  ['--vanilla', 'scripts/maxlfq_iq.R', 'smoke-peptides.csv', 'smoke-maxlfq.csv'],
  ['proteins: 4 | samples: 3 | disconnected: 0']);
}

if (selected.has('bio-proteomics-data-import')) {
const dataImport = path.join(scratch, 'bio-proteomics-data-import');
writeFileSync(path.join(dataImport, 'smoke-proteinGroups.txt'),
  'Protein IDs\tReverse\tPotential contaminant\tOnly identified by site\tLFQ intensity S1\tLFQ intensity S2\n' +
  'P1\t\t\t\t1024\t2048\nP2\t\t\t\t0\t4096\nREV__P3\t+\t\t\t100\t200\nCON__P4\t\t+\t\t100\t200\nP5\t\t\t+\t100\t200\n');
run('bio-proteomics-data-import', 'examples/load_maxquant_qfeatures.R', values.rscript,
  ['--vanilla', 'examples/load_maxquant_qfeatures.R', 'smoke-proteinGroups.txt'],
  ['QFeatures log2 LFQ: 2 protein groups x 2 samples | -Inf: FALSE']);
}

const versions = spawnSync(values.python, ['-c',
  'import sys,importlib.metadata as m; print(sys.version); print({p:m.version(p) for p in ["biopython","numpy","pandas","scipy"]})'], { encoding: 'utf8' });
process.stdout.write(JSON.stringify({ schema_version: 1,
  scope: 'provider_offline_smoke_only_not_aipoch_approval', source_commit: commit,
  isolated_packages: selected.size, scratch_directory: scratch,
  python_versions: versions.stdout?.trim(),
  rscript_version: spawnSync(values.rscript, ['--version'], { encoding: 'utf8' }).stdout?.trim(),
  tests: results, passed: results.filter(test => test.status === 'passed').length,
  failed: results.filter(test => test.status !== 'passed').length,
}, null, 2) + '\n');
if (results.length === 0 || results.some(test => test.status !== 'passed')) process.exitCode = 1;
