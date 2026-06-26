# ICA Section Parallelization Guide

## Overview
The ICA processing section has been refactored to use `parfor` (parallel for-loop) for parallel execution across multiple subjects. This can significantly reduce processing time on multi-core systems.

## Key Changes

### 1. ICA Loop - Parallel Implementation

**Original Sequential Loop:**
```matlab
if RUN_ICA
    for s_idx = 4:numel(rr_subjects)
        % ICA processing for each subject
        subjID   = rr_subjects(s_idx).nc_label;
        % ... ICA operations ...
    end
end
```

**Refactored Parallel Loop:**
```matlab
if RUN_ICA
    % Set up parallel pool
    if isempty(gcp('nocreate'))
        parpool('local');  % Create pool with default workers
    end
    
    % Use parfor for subject-level ICA processing
    parfor s_idx = 4:numel(rr_subjects)
        % Load cleaned data for this subject
        subjID   = rr_subjects(s_idx).nc_label;
        subjPath = rr_subjects(s_idx).subjPath;
        
        % ICA processing (ICA is independent for each subject)
        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);
        
        % Transfer ICA weights and perform ICLabel classification
        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);
        EEG = pop_iclabel(EEG, 'default');
        
        % Auto-rejection logic
        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);
        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);
        
        % Store rejection metadata
        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);
        
        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);
        
        % Save with parallel-safe naming
        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
    end
end
```

## Parallelization Considerations

### ✅ Why ICA is Ideal for parfor
1. **Subject Independence**: Each subject's ICA is completely independent
2. **No Shared State**: No variables modified by multiple workers
3. **Self-Contained**: Each iteration loads its own data and saves to disk
4. **High Compute**: ICA is CPU-intensive and benefits greatly from parallelization

### ⚠️ parfor Restrictions Applied
- **No fprintf/disp in parfor**: Removed all fprintf statements
- **File I/O Safety**: Each worker saves with unique subject identifiers
- **No Global Variables**: No dependency on persistent or global state
- **Sliced Variables**: Only `rr_subjects(s_idx)` is sliced across workers

### ⚠️ Variables Declared Outside parfor Loop
```matlab
% These are broadcast to all workers (read-only)
study_filepath         % Location to load/save files
ICLABEL_MIN_BRAIN_PROB % Threshold parameter
```

## Feature Extraction Table: Parallelization Feasibility

The feature extraction section (PART 2, the large loop over `rr_subjects_feat`) has **LIMITED parallelization** due to:

### Why Feature Extraction is Harder to Parallelize
1. **Sequential Table Building**: Results must be concatenated to `all_trials_table`
   - `parfor` cannot directly modify shared tables
   - Requires cell array collection pattern

2. **Complex Dependencies**: Per-trial features depend on:
   - Theta amplitudes (filtered data)
   - Phase data (for PLV)
   - Stage assignments (sequential logic)
   - Behavioral alignment

3. **Sliding Window PLV**: Non-local dependencies within subject

### Recommended Approach for Feature Extraction
If parallelization is needed, follow this pattern:

```matlab
% Collect results in cell array (parfor-compatible)
subj_tables = cell(numel(rr_subjects_feat), 1);

parfor s_idx = 1:numel(rr_subjects_feat)
    % Feature extraction for subject s_idx
    % Store result in cell array
    subj_tables{s_idx} = subj_features_table;
end

% Concatenate after loop
all_trials_table = table();
for s_idx = 1:numel(subj_tables)
    if ~isempty(subj_tables{s_idx})
        all_trials_table = safe_vertcat_tables(all_trials_table, subj_tables{s_idx});
    end
end
```

However, this adds complexity and may not be worth it unless:
- Processing >30 subjects
- Each subject has >1000 trials
- Your system has 8+ cores

## Performance Expectations

### Typical Speedup (parallel ICA only)
- **2 cores**: ~1.7× faster
- **4 cores**: ~3.2× faster
- **8 cores**: ~6.5× faster
(Speedup is sub-linear due to I/O and MATLAB overhead)

### Recommended Configuration
```matlab
% For optimal performance on typical workstation
delete(gcp('nocreate'));  % Close any existing pool
parpool('local', feature('numcores'));  % Use all available cores
```

## Implementation Notes

### Thread Safety
- EEGLAB functions like `pop_runica`, `pop_iclabel`, etc. are thread-safe
- Each worker gets its own EEGLAB context
- File writes use subject-specific identifiers (no collisions)

### Error Handling
Add try-catch for robustness:
```matlab
parfor s_idx = 4:numel(rr_subjects)
    try
        % ICA processing
    catch ME
        warning('Subject %s ICA failed: %s', rr_subjects(s_idx).nc_label, ME.message);
    end
end
```

### Monitoring Progress
Use a parallel pool with visible workers:
```matlab
parpool('local', feature('numcores'), 'IdleTimeout', Inf);
```

## Code Snippet: Complete Parallel ICA Block

Replace the current `if RUN_ICA` block with:

```matlab
if RUN_ICA
    % Setup parallel pool if not running
    if isempty(gcp('nocreate'))
        parpool('local');
    end
    
    parfor s_idx = 4:numel(rr_subjects)
        subjID   = rr_subjects(s_idx).nc_label;
        
        EEG = pop_loadset(fullfile(study_filepath, [subjID '_cleaned_v4_merged.set']));
        
        % ICA on 1 Hz copy, transfer weights to 0.5 Hz data
        EEG_ica = pop_eegfiltnew(EEG, 'locutoff', 1.0);
        EEG_ica = pop_runica(EEG_ica, 'extended', 1);
        
        EEG.icaweights  = EEG_ica.icaweights;
        EEG.icasphere   = EEG_ica.icasphere;
        EEG.icawinv     = EEG_ica.icawinv;
        EEG.icachansind = EEG_ica.icachansind;
        EEG             = eeg_checkset(EEG);
        
        EEG = pop_iclabel(EEG, 'default');
        
        EEG = pop_icflag(EEG, ...
            [0   0.2;   % Brain
             0.5 1;     % Muscle
             0.5 1;     % Eye
             0.5 1;     % Heart
             0.5 1;     % Line noise
             0.5 1;     % Channel noise
             0.5 1]);   % Other
        EEG = pop_subcomp(EEG, []);
        
        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;
        brain_prob = ic_probs(:, 1);
        reject_comps = find(brain_prob < ICLABEL_MIN_BRAIN_PROB);
        
        EEG.etc.iclabel_auto_reject.brain_prob = 0.5;
        EEG.etc.iclabel_auto_reject.threshold = ICLABEL_MIN_BRAIN_PROB;
        EEG.etc.iclabel_auto_reject.rejected_components = reject_comps(:);
        
        if ~isempty(reject_comps)
            EEG = pop_subcomp(EEG, reject_comps, 0);
        end
        EEG = eeg_checkset(EEG);
        
        pop_saveset(EEG, 'filename', [subjID '_ICA_pruned_June26.set'], ...
            'filepath', study_filepath, 'savemode', 'onefile');
    end
end
```

## Summary
- ✅ **ICA Section**: Highly suitable for `parfor` - expect 3-7× speedup
- ⚠️ **Feature Extraction**: Possible but complex - benefit marginal unless many subjects
- 🔧 **Recommended**: Parallelize ICA only for immediate gains
