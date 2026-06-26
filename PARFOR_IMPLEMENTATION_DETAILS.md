# Detailed Parfor Implementation for ICA Section

## Executive Summary
The ICA processing loop (currently `for s_idx = 4:numel(rr_subjects)`) has been refactored to use `parfor` for **parallel execution across multiple workers**. Expected speedup: **3-7× faster** on typical workstations with 4-8 cores.

---

## Section 1: Why Parfor Works for ICA

### ✅ Independence Criterion
Each subject's ICA processing is **completely independent**:
- No data is shared between subjects
- Each subject's ICA decomposition is unique
- File writes use unique subject identifiers
- No global state modifications

### ✅ Computational Intensity
ICA operations are CPU-intensive:
- `pop_runica()`: Eigenvalue decomposition (high CPU)
- `pop_iclabel()`: Deep learning inference (high CPU)
- Theta filtering & Hilbert transform: High memory/CPU

**Result**: Benefits from parallel execution justify overhead.

### ✅ No Side Effects
The ICA loop doesn't:
- Modify shared variables
- Call external scripts that change state
- Have conditional branching based on other subjects
- Accumulate results (writes individual files)

---

## Section 2: Parfor Restrictions & Solutions

### Restriction 1: No fprintf/Display Inside parfor Loop

**❌ NOT Allowed Inside parfor:**
```matlab
parfor s_idx = 1:10
    fprintf('Processing subject %d\n', s_idx);  % ❌ ERROR!
end
```

**✅ Solution Applied:**
```matlab
% Remove all fprintf/disp statements from inside parfor
% OR use ParforProgressbar (advanced, not used here)

% Progress can be monitored via:
% 1. Checking output files on disk
% 2. Using Parallel Pool monitoring window
% 3. Post-loop summary
```

### Restriction 2: Slicing vs Broadcast Variables

**Broadcast Variables** (read-only on all workers):
```matlab
% These are copied to each worker - OK to use in parfor
study_filepath         % Used to load/save files
ICLABEL_MIN_BRAIN_PROB % Parameter threshold
baseline_win_ms        % Baseline window
do_theta_epochs        % Flag
do_trimming            % Flag
protect_labels         % Cell array of channel labels
```

**Sliced Variables** (distributed across workers):
```matlab
% This is divided into chunks, each chunk processed by one worker
parfor s_idx = 4:numel(rr_subjects)
    % s_idx is automatically sliced
end
```

**❌ NOT Allowed (implicit broadcasting):**
```matlab
parfor s_idx = 1:10
    EEG = EEG_array(s_idx);  % ❌ Trying to modify EEG_array
end
```

**✅ OK (element assignment from struct array):**
```matlab
parfor s_idx = 1:10
    subjID = rr_subjects(s_idx).nc_label;  % ✅ Reading from array element
end
```

### Restriction 3: Output Variables

**❌ Cannot modify variables outside parfor inside the loop:**
```matlab
output_table = [];
parfor s_idx = 1:10
    output_table = [output_table; results];  % ❌ ERROR!
end
```

**✅ Use cell array collection (if needed):**
```matlab
output_cells = cell(10, 1);
parfor s_idx = 1:10
    output_cells{s_idx} = compute_result(s_idx);
end
output_table = vertcat(output_cells{:});
```

**Note**: This ICA implementation uses **file-based output** (no shared table collection inside parfor), so this isn't needed.

---

## Section 3: Implementation Details

### 3.1 Parallel Pool Setup

```matlab
% Check if pool exists; create if not
if isempty(gcp('nocreate'))
    pool = parpool('local');  % Default: uses all available cores
end
```

**What this does:**
- `gcp('nocreate')`: Get current pool without creating one (returns empty if none exists)
- `parpool('local')`: Create local pool (same machine, not cluster)
- Default worker count = logical core count

**Options:**

```matlab
% Explicit 4-worker pool
parpool('local', 4);

% Use all cores
parpool('local', feature('numcores'));

% Use logical cores (including hyper-threading)
parpool('local', feature('numcores', 'logical'));

% Restart pool (close + recreate)
delete(gcp('nocreate'));
parpool('local');
```

### 3.2 Parfor Loop Structure

```matlab
parfor s_idx = 4:numel(rr_subjects)
    % 1. Extract subject metadata
    subjID   = rr_subjects(s_idx).nc_label;
    
    % 2. Load subject data
    EEG = pop_loadset(...);  % Broadcast read from disk
    
    % 3. Process (independent operations)
    EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
    EEG_ica = pop_runica(EEG_ica, 'extended', 1);
    
    % 4. Save results (unique filenames)
    pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...);
    
    % 5. Clear local variables (good practice)
    clear EEG EEG_ica ...
end
```

**Why this structure works:**
- Line 1-2: Extracting metadata is trivial
- Line 4-5: Load is broadcast (same data sent to each worker)
- Line 7-9: Processing is local to worker (no shared state)
- Line 11: Save with unique ID (no file collisions)
- Line 13: Clean up worker memory

### 3.3 Variables in Current Implementation

**Variables accessed inside parfor loop:**

| Variable | Type | Access | Source |
|----------|------|--------|--------|
| `rr_subjects` | struct array | read `s_idx` element | broadcast (main workspace) |
| `study_filepath` | char/string | read full path | broadcast |
| `ICLABEL_MIN_BRAIN_PROB` | double | read value | broadcast |
| `baseline_win_ms` | vector | read | broadcast |
| `do_theta_epochs` | logical | read | broadcast |
| `do_trimming` | logical | read | broadcast |
| `protect_labels` | cell array | read | broadcast |
| `epoch_outpath` | char/string | read | broadcast |
| `outcome_codes` | cell array | read | broadcast |
| `outcome_window` | vector | read | broadcast |
| `outcome_window_epoch` | vector | computed locally | OK |
| `THETA_LO`, `THETA_HI`, `THETA_ORD` | double | read | broadcast |
| `all_trial_data` | struct | read field `(subjID)` | broadcast |
| `misaligned_trigger_IDs` | cell array | read | broadcast |
| `OUTCOME_DELAY_C2_S` | double | read | broadcast |

✅ **All valid** for parfor (read-only, no modifications)

---

## Section 4: Comparison - Before & After

### Before (Sequential)
```matlab
if RUN_ICA
    for s_idx = 4:numel(rr_subjects)
        subjID = rr_subjects(s_idx).nc_label;
        EEG = pop_loadset(...);
        % Process 1 subject
        EEG_ica = pop_runica(...);  % 2-5 minutes per subject
        pop_saveset(...);
        % Next subject waits for this one to finish
    end
end
% For 15 subjects: ~30-75 minutes total (sequential)
```

### After (Parallel with 4 Workers)
```matlab
if RUN_ICA
    if isempty(gcp('nocreate'))
        parpool('local', 4);
    end
    
    parfor s_idx = 4:numel(rr_subjects)
        subjID = rr_subjects(s_idx).nc_label;
        EEG = pop_loadset(...);
        % Process 4 subjects simultaneously
        EEG_ica = pop_runica(...);  % Runs in parallel
        pop_saveset(...);
    end
end
% For 15 subjects: ~8-20 minutes total (4 workers)
% Speedup: ~3.5-4× (not perfect 4× due to I/O & overhead)
```

---

## Section 5: Error Handling

### Current Implementation (as provided)
The parfor loop will stop at first error. For robustness, add try-catch:

```matlab
parfor s_idx = 4:numel(rr_subjects)
    try
        subjID = rr_subjects(s_idx).nc_label;
        % ... processing ...
    catch ME
        % Error handling (limited in parfor)
        % Save error info to file for post-processing
        fid = fopen(fullfile(study_filepath, [subjID '_ERROR.txt']), 'w');
        fprintf(fid, 'Error: %s\n', ME.message);
        fclose(fid);
    end
end
```

**Caveat**: fprintf inside parfor doesn't work, so we write to file instead.

---

## Section 6: Feature Extraction - Why Not Parallelized (Yet)

### Current Feature Extraction Loop
```matlab
for s_idx = 1:numel(rr_subjects_feat)
    % Load data
    % Extract features
    % Append to all_trials_table
end
```

### Problems with Direct Parallelization

**1. Table Concatenation**
```matlab
% Inside parfor - NOT ALLOWED:
all_trials_table = [all_trials_table; subj_features];  % ❌
```

**Solution**: Use cell array pattern
```matlab
subj_tables = cell(numel(rr_subjects_feat), 1);
parfor s_idx = 1:numel(rr_subjects_feat)
    subj_tables{s_idx} = extract_features(s_idx);  % ✅
end
all_trials_table = vertcat(subj_tables{:});
```

**2. Stage Assignment Logic**
Requires sequential reading of `revTrial` vectors and complex masking across subjects.

**3. PLV Computation**
Sliding window requires contextual knowledge of neighboring trials within subject.

**4. Memory Requirements**
Loading all theta/phase EEGs for parallelization would exceed typical RAM.

### Recommendation
**Only parallelize feature extraction if:**
- You have >30 subjects
- Each subject has >1000 trials
- Your system has 8+ cores and >32 GB RAM
- Current runtime is >2 hours per batch

For typical use (15 subjects, 400 trials each), sequential is simpler and nearly as fast.

---

## Section 7: Advanced: Progress Monitoring

### Option 1: Post-Loop Summary (Current)
```matlab
parfor s_idx = 4:numel(rr_subjects)
    % Processing
end
fprintf('ICA complete for %d subjects\n', numel(rr_subjects) - 3);
```

### Option 2: File-Based Progress (Robust)
```matlab
parfor s_idx = 4:numel(rr_subjects)
    % Processing
    % Create marker file when done
    fid = fopen(fullfile(study_filepath, [subjID '_ICA_DONE.txt']), 'w');
    fprintf(fid, 'Done\n');
    fclose(fid);
end

% Check progress
done_files = dir(fullfile(study_filepath, '*_ICA_DONE.txt'));
fprintf('Progress: %d/%d subjects\n', numel(done_files), numel(rr_subjects) - 3);
```

### Option 3: Data Queue (Advanced)
```matlab
if isempty(gcp('nocreate'))
    parpool('local');
end

D = parallel.pool.DataQueue;
afterEach(D, @(x) fprintf('Subject %s ICA done\n', x));

parfor s_idx = 4:numel(rr_subjects)
    subjID = rr_subjects(s_idx).nc_label;
    % ... processing ...
    send(D, subjID);
end
```

---

## Section 8: Troubleshooting

### Issue: "Too many output arguments"
**Cause**: Trying to return values from parfor
**Solution**: Use file I/O instead (already done in this implementation)

### Issue: "Cannot slice variable"
**Cause**: Complex indexing in parfor
**Solution**: Use only simple slicing like `rr_subjects(s_idx)`

### Issue: "Workspace variables are not available"
**Cause**: Forgot to declare variables as broadcast
**Solution**: Ensure all read-only variables are in main workspace before parfor

### Issue: Slow I/O performance
**Cause**: Multiple workers writing to same disk simultaneously
**Solution**: Use unique subject identifiers (already done)

---

## Summary Checklist

- [x] ICA loop is independent (one subject per iteration)
- [x] No shared table/array modifications inside loop
- [x] All file I/O uses unique identifiers
- [x] Broadcast variables are read-only
- [x] No fprintf/disp inside parfor loop
- [x] Parallel pool setup is present
- [x] Error handling considered
- [x] Feature extraction left sequential (not worth parallelizing)

---

## Files Provided

1. **ICA_PARALLELIZATION_GUIDE.md** - High-level overview and recommendations
2. **S2_RR_preprocessing_ICA_PARALLEL_SECTION.m** - Drop-in replacement code
3. **PARFOR_IMPLEMENTATION_DETAILS.md** - This file (detailed technical reference)
