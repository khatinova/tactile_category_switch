# Feature Extraction Parallelization Analysis

## TL;DR
**Yes, it's possible, but NOT recommended.** Here's why:

| Aspect | Difficulty | Benefit | Risk |
|--------|-----------|---------|------|
| **ICA Parallelization** | Easy | 3-5× speedup | Low |
| **Feature Extraction Parallelization** | Hard | 1.5-2× speedup | Medium |

---

## Current Feature Extraction Loop Structure

```matlab
if RUN_FEATURE_EXTRACTION
    for s_idx = 1:numel(rr_subjects_feat)
        % 1. Load EEG data
        % 2. Extract per-trial features
        % 3. Append to all_trials_table
    end
end
```

### What Happens in Each Iteration

**Per Subject (~2-3 minutes):**
1. Load broadband EEG epochs
2. Load theta-filtered epochs
3. Load theta-phase epochs
4. Extract per-trial features (N2, P300, FRN, RewP, theta, PLV)
5. Compute sliding-window PLV (context-dependent)
6. **Concatenate to shared table** ← THE PROBLEM

---

## Why Parallelization is Difficult

### Problem 1: Shared Table Concatenation ⚠️

**Sequential (works fine):**
```matlab
for s_idx = 1:numel(rr_subjects_feat)
    subj_features = extract_features(s_idx);
    all_trials_table = [all_trials_table; subj_features];  % ✅ Direct append
end
```

**Parallel (doesn't work directly):**
```matlab
parfor s_idx = 1:numel(rr_subjects_feat)
    subj_features = extract_features(s_idx);
    all_trials_table = [all_trials_table; subj_features];  % ❌ ERROR!
    % parfor can't modify shared variables
end
```

**Workaround (cell array pattern):**
```matlab
subj_tables = cell(numel(rr_subjects_feat), 1);
parfor s_idx = 1:numel(rr_subjects_feat)
    subj_tables{s_idx} = extract_features(s_idx);  % ✅ OK (cell indexing allowed)
end
all_trials_table = table();
for s_idx = 1:numel(subj_tables)
    all_trials_table = [all_trials_table; subj_tables{s_idx}];
end
```

**Cost:** Additional memory (whole table in cell array before concatenating)

### Problem 2: Stage Assignment Logic

```matlab
% This sequential logic is complex and context-dependent
for b = 1:height(beh.correct)
    beh_correct = [beh_correct, beh.correct(b, :)];
end

rev_trials_vec = beh.revTrial(:);  % Per-subject reversal point

T = assign_stages_preserve_LE_RN(...
    T, block_col, trial_col, rev_trials_vec, ...);
```

**Issue:** Each subject's stage assignment depends on their reversal trial
- Can't compute in parallel (need sequential context)
- Could move outside parfor, but adds complexity

### Problem 3: Sliding Window PLV

```matlab
for bi2 = 1:n_bucket
    win_lo = max(1, bi2 - PLV_WINDOW_HALF);
    win_hi = min(n_bucket, bi2 + PLV_WINDOW_HALF);
    window_rows = bucket_rows(win_lo:win_hi);  % ← Context-dependent
```

**Issue:** PLV calculation requires neighboring trials in same subject
- Can't safely parallelize across subjects
- Could parallelize within-subject, but complex

### Problem 4: Memory Requirements

Loading all epochs for all subjects simultaneously:
```
Per subject: 
  - Broadband EEG: ~200 MB (1000 trials × 1500 timepoints × 128 channels)
  - Theta EEG: ~200 MB
  - Phase data: ~200 MB
  - Total per subject: ~600 MB

For 15 subjects with 4 workers:
  - 4 workers × 600 MB = ~2.4 GB just for data
  - Plus MATLAB overhead: ~3.5-4 GB
  - Typical workspace: 8-16 GB (tight!)
```

---

## Current Feature Extraction Runtime

**Typical Timing (15 subjects, 400 trials each):**

```
Sequential: 30-40 minutes
  ├─ Load EEG files: 5 min (disk I/O)
  ├─ Extract features: 25-30 min (computation)
  └─ Stage/PLV/table building: 5 min

Parallel (if optimized): 20-30 minutes
  ├─ Load EEG files: 2 min (4 workers, parallel I/O)
  ├─ Extract features: 10-15 min (4× compute)
  └─ Stage/PLV/table building: 8-12 min (sequential bottleneck)
```

**Speedup: ~1.5-2× (not impressive given complexity)**

---

## Implementation Path (If You Want to Try)

### Step 1: Refactor Feature Extraction Function

Extract into standalone function:
```matlab
function subj_table = extract_subject_features(subj, rr_subjects_feat, ...
    epoch_outpath, all_trial_data, RR_behav_table, ...)
    % Self-contained: loads, processes, returns single table
    % No side effects, no global modifications
    subj_table = table();
    % ... all feature extraction logic ...
end
```

### Step 2: Pre-compute Stage Assignment

Before parfor:
```matlab
% Compute stages for all subjects sequentially
stage_tables = cell(numel(rr_subjects_feat), 1);
for s_idx = 1:numel(rr_subjects_feat)
    stage_tables{s_idx} = assign_stages(...);
end
```

### Step 3: Parallelize with Cell Array Collection

```matlab
subj_tables = cell(numel(rr_subjects_feat), 1);

parfor s_idx = 1:numel(rr_subjects_feat)
    subj_table = extract_subject_features(...);
    subj_tables{s_idx} = subj_table;  % ✅ Cell indexing allowed
end

% Concatenate after loop
all_trials_table = table();
for s_idx = 1:numel(subj_tables)
    all_trials_table = [all_trials_table; subj_tables{s_idx}];
end
```

### Step 4: Handle FRN/RewP Computation

These are difference waves (must be computed per stage after all subjects):
```matlab
% After all_trials_table is built
frn_rewp_stage_table = kh_compute_frn_rewp_by_stage(...);
```

---

## Detailed Comparison: Sequential vs Parallel

### Sequential Version (Current)
```
Pros:
  ✅ Simple, readable code
  ✅ Straightforward table concatenation
  ✅ Easy debugging
  ✅ Minimal memory overhead
  ✅ No synchronization issues

Cons:
  ❌ 30-40 minutes runtime
  ❌ CPU underutilized (1 core)
```

### Parallel Version (Proposed)
```
Pros:
  ✅ 1.5-2× speedup (20-30 minutes)
  ✅ Better CPU utilization
  ✅ Still feasible with 8+ GB RAM

Cons:
  ❌ Complex refactoring required (200+ lines changed)
  ❌ Cell array pattern (less elegant)
  ❌ Stage assignment logic moved (adds complexity)
  ❌ Memory usage increases 3-4×
  ❌ Harder to debug (distributed execution)
  ❌ Limited speedup (1.5-2× only)
  ❌ Risk of subtle bugs in table concatenation
```

---

## Reality Check: Is It Worth It?

### Current Situation
- Sequential feature extraction: **30-40 minutes** ✓ Acceptable
- Sequential ICA: **15-30 minutes** ✓ Acceptable
- **Total pipeline: 50-70 minutes** ✓ Reasonable for 15 subjects

### After ICA Parallelization (Just parfor ICA)
- **ICA parallelized: 4-8 minutes** (3-5× speedup)
- Feature extraction: 30-40 minutes (unchanged)
- **Total pipeline: 35-50 minutes** ⚡ Good improvement!

### After Both Parallelized
- ICA parallelized: 4-8 minutes
- Feature extraction parallelized: 20-30 minutes (1.5-2× speedup)
- **Total pipeline: 25-40 minutes** ⚡⚡ Marginal improvement

**Real-world question:** Is saving 10-15 minutes worth 200+ lines of complex refactoring, increased risk, and harder debugging?

---

## Recommendation

### ✅ DO Parallelize ICA
- Easy to implement (3 changes)
- Massive speedup (3-5×)
- Low risk (independent subjects)
- Already provided in PR #4

### ⚠️ DON'T Parallelize Feature Extraction
**Unless you have:**
1. >30 subjects (not 15)
2. Each with >1000 trials (you have 400)
3. Tight time constraints
4. 16+ GB RAM available
5. Experienced team for maintenance

---

## If You Really Want to Try It Anyway

I can provide:
1. **Feature extraction refactoring code** (standalone function)
2. **Parallel cell-array wrapper** (parfor boilerplate)
3. **Memory optimization guide** (reduce per-subject footprint)
4. **Debugging toolkit** (monitor parallel workers)

**Effort estimate:** 2-3 hours implementation + 1-2 hours testing/debugging

---

## Alternative: Optimize Within Sequential Loop

Faster improvements without parallelization complexity:

```matlab
% 1. Cache all_trial_data lookups
subj_data_cache = all_trial_data.(subj);

% 2. Pre-allocate tables
subj_features = table();
subj_features.N2_amp = nan(n_rows, 1);  % Pre-allocate instead of growing

% 3. Vectorize operations where possible
% Instead of: for ti = 1:n_rows ... compute feature ...
% Use: vectorized operations on batch

% Expected gain: 10-15% speedup, zero risk
```

**This alone might get you 30-40 min → 25-35 min without parallelization complexity!**

---

## Summary Table

| Optimization | Speedup | Difficulty | Risk | Effort | Recommended |
|--------------|---------|-----------|------|--------|-------------|
| ICA parfor | 3-5× | Easy | Low | 15 min | ✅ YES |
| Feature extract parfor | 1.5-2× | Hard | Medium | 2-3 hrs | ⚠️ Maybe |
| Sequential optimizations | 1.1-1.2× | Easy | Low | 1 hr | ✅ YES |
| Do nothing | 1× | - | - | 0 | ⚠️ Current |

---

## My Honest Opinion

**For your use case (15 subjects):**

1. ✅ **Implement ICA parfor** (already done in PR #4)
   - Expected total time: 35-50 min (vs 50-70 min)
   - Effort: 15 minutes
   - Risk: Low

2. ✓ **Consider sequential optimizations** (caching, pre-allocation)
   - Expected: 25-35 min total
   - Effort: 1 hour
   - Risk: None

3. ❌ **Skip feature extraction parfor** (for now)
   - Marginal benefit (1.5-2×)
   - High implementation complexity
   - Maintenance burden
   - Better to wait until you have 30+ subjects

---

## Decision Tree

```
Do you have 30+ subjects?
  └─ YES → Consider feature extraction parfor (worth the effort)
  └─ NO  → Skip it (see below)

Do you have tight time constraints?
  └─ YES → Consider sequential optimizations (1 hour, 1.1-1.2× speedup)
  └─ NO  → You're fine with 35-50 min after ICA parfor

Is your current pipeline acceptable?
  └─ YES → Just do ICA parfor (PR #4), you're done
  └─ NO  → Add sequential optimizations
```

---

## If You Change Your Mind Later

I can provide feature extraction parallelization when you:
- Scale to 30+ subjects
- Have time for refactoring
- Need the extra speedup

Just ask! 🚀
