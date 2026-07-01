# RR EEG Preprocessing and Feature Extraction Pipeline

## Overview

I've created `S2_RR_preprocess_epoch_eeg.m` - a dedicated script for preprocessing and feature extraction of your RR (EGI) data. This script mirrors the KH pipeline but uses RR-specific configurations.

## Key Adaptations Made

### 1. Channel Mapping (RR vs KH)
- **KH**: FCz, Cz, Pz, P1, P2
- **RR**: E11 (FCz equivalent), E7 (Cz equivalent), E62, E67, E72 (parietal)
- **Somatosensory**: E36, E104, E41, E103

### 2. Event Codes
- **Current**: `{'rewa','puni','corr','inco'}` (example)
- **⚠️ VERIFY**: Check your actual RR event codes and update line ~290

### 3. Reference Scheme
- **RR**: Cz (E7/channel 61) first, then common average reference
- **KH**: Direct average reference

### 4. No Mastoid Regression
- RR data doesn't have LM/RM channels, so this step is skipped

## Files Created/Modified

1. **`pipeline/S2_RR_preprocess_epoch_eeg.m`** - New RR preprocessing script
2. **`pipeline/S4_merge_feature_tables.m`** - Updated to load from new locations

## Steps to Run

### 1. First, configure the RR script:

```matlab
cd /path/to/tactile_category_switch/pipeline
open S2_RR_preprocess_epoch_eeg.m
```

**Edit these settings:**
- Line 44: `RR_VALID_PARTICIPANTS = [1:10];` - Update with your actual RR participant numbers
- Line 290: `outcome_codes = {'rewa','puni','corr','inco'};` - **VERIFY with your RR event structure**
- Line 134: Check subject folder pattern (`Nc*` vs `RR*`)

### 2. Run the RR pipeline:

```matlab
% Step 2a: Basic preprocessing and feature extraction
S2_RR_preprocess_epoch_eeg

% Step 2b: Enhanced feature extraction (peak latencies, advanced PLV, z-scores)
S3_RR_extract_eeg_features
```

**S2_RR** will:
- Preprocess all RR subjects
- Generate epoched data in: `Salient mod switch RR/Results/EEG analysis/Epoched_data_noisefiltering/`
- Create basic feature tables in: `Salient mod switch RR/Results/EEG analysis/Outcome_feature_tables_v4_merged/`

**S3_RR** will:
- Add peak latencies (E11_neg_peak_lat, P300_peak_lat)
- Compute enhanced cross-trial PLV
- Add within-subject z-scores
- Generate FRN/RewP difference waves
- Create enhanced feature tables and grand average plots

### 3. Merge KH and RR data:

```matlab
S4_merge_feature_tables
```

This will combine KH and RR feature tables into `group_feature_table_combined.mat`

## Output Structure

```
Salient mod switch RR/Results/EEG analysis/
├── Epoched_data_noisefiltering/
│   ├── Nc01_outcome.set, Nc01_outcome_trimmed.set
│   ├── Nc01_outcome_theta.set, Nc01_outcome_theta_trimmed.set
│   └── Nc01_trial2epoch.mat
├── Outcome_feature_tables_v4_merged/
│   ├── RR_outcome_all_trials_v4_merged.mat/.csv (from S2_RR)
│   ├── RR_outcome_all_trials_v4_merged_enhanced.mat/.csv (from S3_RR)
│   ├── RR_outcome_stage_features_v4_merged.mat/.csv (from S2_RR)
│   └── RR_frn_rewp_by_stage.mat/.csv (from S3_RR)
└── Figures/
    └── outcome_v4_merged_QC/
        └── Nc01_RR_merged_v4_outcome_QC.pdf
```

## Important Notes

### ⚠️ Before Running:
1. **Verify RR event codes** - Check your .mff files for actual event markers
2. **Check subject naming** - Confirm if your RR subjects use "Nc01" or "RR01" pattern
3. **Validate channel locations** - Ensure E11, E7, E62, etc. exist in your montage

### 🔧 If You Get Errors:
1. **"Channel E11 not found"**: Check your EGI montage - you might need different channel numbers
2. **"No .mff file found"**: Verify the subject folder structure matches the script expectations
3. **"Missing all_trial_data"**: Ensure S1 (behavior script) was run first and includes RR subjects

### 🎯 Key Features Extracted:
**S2_RR (Basic)**:
- **ERP**: N2, P300, E11/E7 combined waveform
- **Theta**: 4-8 Hz amplitude at E11
- **PLV**: Basic fronto-parietal and fronto-somatosensory phase-locking
- **FRN/RewP**: Stage-level difference waves

**S3_RR (Enhanced)**:
- **Peak latencies**: E11_neg_peak_lat, P300_peak_lat
- **Advanced PLV**: Condition-specific cross-trial computation
- **Z-scores**: Within-subject normalized features
- **Grand averages**: ERP plots and difference wave visualizations

The output format matches your KH pipeline exactly, so all downstream analysis scripts should work with the merged data.