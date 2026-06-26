# Before & After Code Comparison: ICA Section

## Overview
This document shows the exact differences between the sequential ICA loop and the new parallel version.

---

## Original Code (Sequential, lines ~280-350)

```matlab
if RUN_ICA

    for s_idx = 4:numel(rr_subjects)

        subjID   = rr_subjects(s_idx).nc_label;   % clean 'Nc##' used for all outputs
        subjPath = rr_subjects(s_idx).subjPath;

        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));
        % -----------------------------------------------------------------
        % 9. ICA on 1 Hz copy, transfer weights to 0.5 Hz data.
        % -----------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);

        EEG = pop_iclabel(EEG, 'default');

        % -----------------------------------------------------------------
        % ICLabel auto-rejection requested:
        % Commented-out/manual route removed. Reject all components with
        % Brain probability < 0.50.
        % -----------------------------------------------------------------

        % OLD MANUAL ROUTE:
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other
        % EEG = review_iclabel_interactive(EEG, 'nCols', 10);
        EEG = pop_subcomp(EEG, []);

        if ~isfield(EEG.etc, 'ic_classification') || ...
           ~isfield(EEG.etc.ic_classification, 'ICLabel') || ...
           ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            error('%s: ICLabel classifications not found after pop_iclabel.', subjID);
        end

        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);

        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);

        fprintf('  ICLabel auto-rejection: rejecting %d/%d ICs with Brain probability < %.2f.\n', ...
            numel(reject_comps), numel(brain_prob), ICLABEL_MIN_BRAIN_PROB);

        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);

        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');

        % --- REST OF LOOP (theta filtering, epoching, etc.) ---
        % [~400 more lines of code, processing continues below]

    end
end
```

---

## Parallelized Code (parfor version)

```matlab
if RUN_ICA

    % -------------------------------------------------------------------------
    % NEW: Setup Parallel Pool
    % -------------------------------------------------------------------------
    if isempty(gcp('nocreate'))
        pool = parpool('local');
    end

    % -------------------------------------------------------------------------
    % KEY CHANGE: Use parfor instead of for
    % -------------------------------------------------------------------------
    parfor s_idx = 4:numel(rr_subjects)

        subjID   = rr_subjects(s_idx).nc_label;
        subjPath = rr_subjects(s_idx).subjPath;

        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));
        % -----------------------------------------------------------------
        % 9. ICA on 1 Hz copy, transfer weights to 0.5 Hz data.
        % -----------------------------------------------------------------
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);

        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);

        EEG = pop_iclabel(EEG, 'default');

        % -----------------------------------------------------------------
        % ICLabel auto-rejection
        % -----------------------------------------------------------------
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other
        EEG = pop_subcomp(EEG, []);

        if ~isfield(EEG.etc, 'ic_classification') || ...
           ~isfield(EEG.etc.ic_classification, 'ICLabel') || ...
           ~isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            error('%s: ICLabel classifications not found after pop_iclabel.', subjID);
        end

        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);

        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);

        % NOTE: fprintf removed - NOT allowed inside parfor loop
        % Progress can be monitored via output files or pool window

        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);

        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);

        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');

        % --- REST OF LOOP (theta filtering, epoching, etc.) ---
        % [Identical to sequential version, runs in parallel per worker]

    end  % END parfor (not 'end' for regular for loop)
end
```

---

## Detailed Change Summary

### Change 1: Pool Setup (New, ~3 lines)
```matlab
% ADDED: Before parfor loop
if isempty(gcp('nocreate'))
    pool = parpool('local');
end
```
- `gcp('nocreate')`: Get current pool (returns [] if none)
- `parpool('local')`: Create pool with default worker count
- **Effect**: All subsequent iterations run in parallel across workers

### Change 2: Loop Declaration (1 line)
```matlab
% BEFORE
for s_idx = 4:numel(rr_subjects)

% AFTER
parfor s_idx = 4:numel(rr_subjects)
```
- **Effect**: Loop iterations run in parallel instead of sequentially
- Loop index `s_idx` is automatically distributed across workers

### Change 3: Remove fprintf (1 line removed)
```matlab
% BEFORE
fprintf('  ICLabel auto-rejection: rejecting %d/%d ICs with Brain probability < %.2f.\n', ...
    numel(reject_comps), numel(brain_prob), ICLABEL_MIN_BRAIN_PROB);

% AFTER
% (Line removed - fprintf not allowed in parfor)
```
- **Why**: Display functions (fprintf, disp, disp) can't run inside parfor
- **Alternative**: Monitor via file writes or pool window

### Everything Else
- Identical to original code
- All processing logic remains unchanged
- File I/O patterns remain unchanged (no collisions because subjID is unique)

---

## Performance Impact

### Timing Estimate
Assuming:
- 12 subjects (s_idx = 4:15)
- 2-3 minutes per subject for ICA + filtering + epoching
- 4-core system

| Config | Time | Notes |
|--------|------|-------|
| Sequential (original) | 24-36 min | Each subject waits for previous |
| Parallel, 2 workers | 12-20 min | ~1.8× speedup |
| Parallel, 4 workers | 6-10 min | ~3.5× speedup (overhead cost) |
| Parallel, 8 workers | 4-7 min | ~5× speedup (I/O bottleneck) |

**Speedup Formula**: `Expected ≈ Ideal / (1 + (K * overhead_fraction))`
- Ideal = number of workers
- K = tuning factor (typically 1-2)
- Most overhead from I/O, not computation

---

## Migration Checklist

- [ ] Backup original script
- [ ] Replace `for s_idx = 4:numel(rr_subjects)` with `parfor s_idx = 4:numel(rr_subjects)`
- [ ] Add pool setup code before parfor:
  ```matlab
  if isempty(gcp('nocreate'))
      pool = parpool('local');
  end
  ```
- [ ] Remove fprintf statements inside parfor (or comment out)
- [ ] Test with subset first (change 4:numel to 4:6 for testing)
- [ ] Monitor output files to verify parallel execution
- [ ] Compare runtime with original version
- [ ] Optional: Add error handling (try-catch inside parfor)

---

## Verification Steps

### Step 1: Check Parallel Execution
Look for multiple .set files being saved simultaneously:
```bash
# In a terminal, while script runs:
watch -n 1 'ls -l study_filepath/*_ICA_pruned_June26.set | tail -5'
```
If files appear rapidly and out-of-order (not 4,5,6,7...), parallel execution is working.

### Step 2: Monitor Pool
```matlab
% In MATLAB, in another window
poolobj = gcp;
p = poolobj.plot();  % Shows worker activity
```

### Step 3: Compare Outputs
Sequential and parallel versions should produce **identical** .set files:
```matlab
% In MATLAB, after both runs
s1 = load([filename '_S1.mat']);
s2 = load([filename '_S2.mat']);
whos(s1, s2)  % Compare sizes/types
```

---

## Reverting to Sequential

If you need to go back to sequential:

```matlab
% Change
parfor s_idx = 4:numel(rr_subjects)

% Back to
for s_idx = 4:numel(rr_subjects)

% Remove pool setup
% if isempty(gcp('nocreate'))
%     pool = parpool('local');
% end

% Restore fprintf if desired
fprintf('Processing %s...\n', subjID);
```

---

## Advanced: Custom Worker Count

```matlab
% Default (all logical cores)
parpool('local');

% Specific number
parpool('local', 4);

% Specific number, with custom settings
pool = parpool('local', 4);
pool.IdleTimeout = Inf;  % Don't shut down idle pool
pool.Description = 'ICA Processing';

% Check actual worker count
p = gcp;
fprintf('Using %d workers\n', p.NumWorkers);

% Gracefully shutdown
delete(gcp('nocreate'));
```

---

## Summary of Changes

| Aspect | Before | After |
|--------|--------|-------|
| Loop type | `for` | `parfor` |
| Pool setup | None | `parpool('local')` |
| Worker count | 1 (serial) | 2-8 (system dependent) |
| fprintf inside loop | Yes | No (removed) |
| Execution pattern | Sequential | Parallel |
| Output files | Saved serially | Saved in parallel |
| Total runtime | 24-36 min | 6-10 min (4 workers) |
| Code readability | Simple | Slightly complex |
| Debugging | Easy | Harder (distributed) |

---

## Files for Implementation

1. **S2_RR_preprocessing_ICA_PARALLEL_SECTION.m**
   - Complete parfor code block
   - Drop-in replacement for lines ~280-350 in original script

2. **ICA_PARALLELIZATION_GUIDE.md**
   - High-level overview
   - Feature extraction discussion
   - Performance expectations

3. **PARFOR_IMPLEMENTATION_DETAILS.md**
   - Detailed technical reference
   - Parfor restrictions explained
   - Error handling patterns

4. **BEFORE_AFTER_CODE_COMPARISON.md** (this file)
   - Side-by-side comparison
   - Migration checklist
   - Verification steps
