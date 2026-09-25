// Provider preflight; never grants Marketplace approval or executes Skill code.
// Usage: node authoring/verify-pilot.mjs --marketplace <clone> --audit <clone> [--intake <generated-json>]
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { existsSync, readFileSync, readdirSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { parseArgs } from 'node:util';

const { values } = parseArgs({ options: {
  marketplace: { type: 'string' }, audit: { type: 'string' }, intake: { type: 'string' },
} });
assert(values.marketplace && values.audit, '--marketplace and --audit are required');
const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const base = path.join(root, 'authoring/submissions/mrsonord2240');
const evidenceBytes = readFileSync(path.join(base, 'pilot-review-evidence.json'));
const evidence = JSON.parse(evidenceBytes);
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const git = (repo, ...args) => execFileSync('git', ['-C', repo, ...args], { maxBuffer: 32 * 1024 * 1024 });
const blob = (repo, commit, file) => git(repo, 'show', `${commit}:${file}`);
const audit = evidence.external_audit;
const source = evidence.submitted_source;
const provenanceBytes = blob(root, source.commit, 'PROVENANCE.json');
const provenance = JSON.parse(provenanceBytes);
assert.equal(evidence.original_source.repository, provenance.sources['gptomics-bioskills'].upstream_repository);
assert.equal(evidence.original_source.commit, provenance.sources['gptomics-bioskills'].upstream_commit);
const configs = readdirSync(base, { withFileTypes: true })
  .filter(entry => entry.isDirectory())
  .map(entry => path.join(base, entry.name, 'release.config.json'));
assert.equal(configs.length, 20);
assert.equal(evidence.skills.length, 20);
assert.equal(new Set(evidence.skills.map(skill => skill.id)).size, 20);
for (const commit of [source.commit, audit.commit, evidence.original_source.commit])
  assert.match(commit, /^[0-9a-f]{40}$/);
const auditOrigin = git(values.audit, 'config', '--get', 'remote.origin.url').toString().trim().replace(/\.git$/, '');
assert.equal(auditOrigin, audit.repository);
const intakeArgs = configs.flatMap(config => ['--manifest', config, '--source', root]);
const generated = JSON.parse(execFileSync(process.execPath,
  [path.resolve(values.marketplace, 'scripts/intake-skill.mjs'), ...intakeArgs],
  { cwd: values.marketplace, maxBuffer: 8 * 1024 * 1024 }));
const localIntake = values.intake ? path.resolve(values.intake) : path.join(base, 'intake-review-input.json');
if (existsSync(localIntake)) {
  const saved = JSON.parse(readFileSync(localIntake));
  assert.deepEqual(generated, saved, 'Saved intake evidence differs from current official intake');
}
const allAuditPaths = git(values.audit, 'ls-tree', '-r', '--name-only', audit.commit)
  .toString().trim().split('\n');
const skills = [];
for (const skill of evidence.skills) {
  const config = JSON.parse(readFileSync(path.join(base, skill.id, 'release.config.json')));
  assert.equal(config.id, skill.id);
  assert.deepEqual(config.source, { repository: source.repository, commit: source.commit, path: `skills/${skill.id}` });
  assert.equal(skill.intake_record, `${skill.id}@${config.version}`);
  assert(generated[skill.intake_record]);
  const upstream = provenance.skills.find(entry => entry.id === skill.id);
  assert(upstream);
  assert.equal(upstream.upstream_path, skill.provenance.upstream_path);
  assert.equal(upstream.relative_to_upstream, skill.provenance.modification_status);
  assert.equal(upstream.score, skill.upstream_audit.score);
  assert.equal(upstream.audited_on, skill.upstream_audit.audited_on);
  assert.equal(upstream.fix_log, skill.upstream_audit.fix_log);
  assert.equal(skill.aipoch_review.status, 'pending_independent_review');
  for (const field of ['reviewed_by', 'reviewed_on', 'score'])
    assert.equal(skill.aipoch_review[field], null, 'Provider preflight must not fabricate approval');
  for (const scope of skill.aipoch_review.review_scope)
    assert(Object.hasOwn(evidence.review_scope_definitions, scope));
  assert(skill.aipoch_review.unresolved_limitations.length > 0);
  const prefix = `${config.source.path}/`;
  const files = git(root, 'ls-tree', '-r', '--name-only', source.commit, '--', config.source.path)
    .toString().trim().split('\n');
  const executables = files.filter(file => /\.(py|sh|r)$/i.test(file));
  assert.deepEqual([...skill.contract.bundled_executables].sort(),
    executables.map(file => file.slice(prefix.length)).sort(), 'Incomplete executable inventory');
  // Working-tree syntax/runtime checks can only be bound to the release if every
  // packaged local file equals its Git blob. No newline normalization is allowed.
  for (const file of files)
    assert(readFileSync(path.join(root, file)).equals(blob(root, source.commit, file)), `Working file differs: ${file}`);
  const frontmatter = blob(root, source.commit, `${prefix}SKILL.md`).toString().split(/^---\s*$/m)[1];
  assert.match(frontmatter, /^\s*tool_type:/m);
  assert.match(frontmatter, /^\s*primary_tool:/m);
  assert.match(frontmatter, /^\s*author: GPTomics\s*$/m);
  assert.doesNotMatch(frontmatter, /^\s*category:/m);
  const reports = allAuditPaths.filter(file => file.startsWith(`audits/skills/${skill.id}/`) && file.endsWith('/report.json'));
  const matches = reports.map(file => ({ file, bytes: blob(values.audit, audit.commit, file) }))
    .map(item => ({ ...item, report: JSON.parse(item.bytes) }))
    .filter(item => item.report.final?.score === upstream.score && item.report.meta?.evaluated_on === upstream.audited_on);
  assert.equal(matches.length, 1, `Expected one matching upstream report: ${skill.id}`);
  const match = matches[0];
  const recordPath = match.file.replace(/report\.json$/, 'record.json');
  const recordBytes = blob(values.audit, audit.commit, recordPath);
  const record = JSON.parse(recordBytes);
  assert.equal(record.source.author, 'GPTomics');
  assert.match(record.source.commit, /^[0-9a-f]{40}$/);
  const auditEvidence = {
    report_path: match.file, report_sha256: hash(match.bytes),
    report_url: `${audit.repository}/blob/${audit.commit}/${match.file}`,
    record_path: recordPath, record_sha256: hash(recordBytes),
    record_url: `${audit.repository}/blob/${audit.commit}/${recordPath}`,
    fix_log_sha256: hash(blob(values.audit, audit.commit, upstream.fix_log)),
    fix_log_url: `${audit.repository}/blob/${audit.commit}/${upstream.fix_log}`,
    audited_source: record.source,
    auditor_independent: match.report.meta?.auditor_independent ?? null,
    submitted_byte_equivalence: 'unverified_different_repository_and_commit',
  };
  if (skill.upstream_audit.evidence)
    assert.deepEqual(skill.upstream_audit.evidence, auditEvidence);
  assert.deepEqual(skill.package_evidence, generated[skill.intake_record]);
  for (const test of skill.provider_review?.runtime_checks ?? []) {
    assert.equal(test.id, skill.id);
    assert.equal(test.source_commit, source.commit);
    assert(skill.contract.bundled_executables.includes(test.script));
    assert.equal(test.source_sha256, hash(blob(root, source.commit, `${prefix}${test.script}`)));
    if (test.status === 'passed') {
      assert.equal(test.exit_code, 0);
      assert.equal(test.error, null);
      assert((test.checks ?? []).every(check => check.passed));
    }
  }
  skills.push({ id: skill.id, package_evidence: generated[skill.intake_record],
    upstream_audit_evidence: auditEvidence,
    executable_sha256: Object.fromEntries(executables.map(file => [file.slice(prefix.length), hash(blob(root, source.commit, file))])),
    packaged_files_equal_worktree: true });
}
if (evidence.provider_review) {
  const tests = evidence.skills.flatMap(skill => skill.provider_review.runtime_checks);
  assert.deepEqual(evidence.provider_review.functional_checks, {
    total: tests.length, passed: tests.filter(test => test.status === 'passed').length,
    failed: tests.filter(test => test.status !== 'passed').length,
    skills_exercised: new Set(tests.map(test => test.id)).size,
  });
}
process.stdout.write(JSON.stringify({
  schema_version: 1, scope: 'provider_preflight_not_marketplace_approval',
  source_commit: source.commit, audit_commit: audit.commit,
  marketplace_tool_commit: git(values.marketplace, 'rev-parse', 'HEAD').toString().trim(),
  evidence_sha256: hash(evidenceBytes), provenance_sha256: hash(provenanceBytes),
  intake_records_verified: skills.length, skills,
}, null, 2) + '\n');
