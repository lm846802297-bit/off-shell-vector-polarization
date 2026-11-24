# Off-shell vector polarization storage audit — 2026-11-24

**Scope:** read-only current-tree audit and adoption of the central three-tier scientific data storage policy.  
**Repository:** `lm846802297-bit/off-shell-vector-polarization`  
**Audited branch:** `main`

## 1. Repository baseline

- GitHub repository metadata size: **686 KiB**.
- Complete recursive main tree: **4 blobs**, approximately **0.0138 MiB** summed file bytes.
- Recursive tree response was complete (`truncated=false`).
- No file exceeds 1 MiB.
- The current tree contains project instructions, README, and standards lock; no numerical result bundle is present.

## 2. Storage-tier decision

- **Ordinary Git:** appropriate for all current theory notes, project instructions, derivations, manifests, and future small/medium native numerical results.
- **Git LFS:** not enabled; no current data object justifies it.
- **External object storage:** not configured; no bulk or long-term dataset currently requires it.

For future numerical work, retain native output, exact source/input identities, run environment, parameters, convergence and independent physical checks. Use immutable result directories. Select LFS or external storage only when actual data size, access pattern, licensing, or retention needs justify it.

## 3. Integrity and migration boundary

- No existing file was modified, deleted, or moved.
- No Git history was rewritten.
- No scientific calculation was run or rerun.
- No LFS or external storage was configured.

## 4. Standards adoption

The downstream lock pins central commit `14155487afc075d301fa7f6cf452c7bff137e869`, which includes central PR #42's three-tier Git / Git LFS / external object storage policy. Off-shell current definitions, gauge/BRST conventions, Slavnov-Taylor identities, pole residues, and polarization-density-matrix requirements remain project-local.

## 5. Conclusion

This is a very small theory repository with no current numerical data bundle. Ordinary Git is appropriate. **No migration and no LFS activation are required.** Apply the three-tier policy when substantial native numerical outputs are added.
