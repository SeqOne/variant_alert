# coding=utf-8
# =========================================================================
# Modules
# =========================================================================

ACMG_CLASS_BREAKING_CHANGE = {
    "other": {
        "..": "warning",
        "other": "null",
        "Benign": "null",
        "Benign/Likely_benign": "null",
        "Likely_benign": "null",
        "Uncertain_significance": "unknown",
        "Conflicting_interpretations_of_pathogenicity": "unknown",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "..": {
        "..": "null",
        "other": "null",
        "Benign": "null",
        "Benign/Likely_benign": "null",
        "Likely_benign": "null",
        "Uncertain_significance": "unknown",
        "Conflicting_interpretations_of_pathogenicity": "unknown",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Benign": {
        "..": "warning",
        "other": "null",
        "Benign": "null",
        "Benign/Likely_benign": "null",
        "Likely_benign": "null",
        "Uncertain_significance": "unknown",
        "Conflicting_interpretations_of_pathogenicity": "unknown",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Benign/Likely_benign": {
        "..": "warning",
        "other": "null",
        "Benign": "null",
        "Benign/Likely_benign": "null",
        "Likely_benign": "null",
        "Uncertain_significance": "unknown",
        "Conflicting_interpretations_of_pathogenicity": "unknown",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Likely_benign": {
        "..": "warning",
        "other": "null",
        "Benign": "null",
        "Benign/Likely_benign": "null",
        "Likely_benign": "null",
        "Uncertain_significance": "unknown",
        "Conflicting_interpretations_of_pathogenicity": "unknown",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Uncertain_significance": {
        "..": "warning",
        "other": "minor",
        "Benign": "minor",
        "Benign/Likely_benign": "minor",
        "Likely_benign": "minor",
        "Uncertain_significance": "null",
        "Conflicting_interpretations_of_pathogenicity": "null",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Conflicting_interpretations_of_pathogenicity": {
        "..": "warning",
        "other": "minor",
        "Benign": "minor",
        "Benign/Likely_benign": "minor",
        "Likely_benign": "minor",
        "Uncertain_significance": "null",
        "Conflicting_interpretations_of_pathogenicity": "null",
        "Likely_pathogenic": "major",
        "Pathogenic/Likely_pathogenic": "major",
        "Pathogenic": "major",
    },
    "Likely_pathogenic": {
        "..": "warning",
        "other": "major",
        "Benign": "major",
        "Benign/Likely_benign": "major",
        "Likely_benign": "major",
        "Uncertain_significance": "major",
        "Conflicting_interpretations_of_pathogenicity": "major",
        "Likely_pathogenic": "null",
        "Pathogenic/Likely_pathogenic": "null",
        "Pathogenic": "minor",
    },
    "Pathogenic/Likely_pathogenic": {
        "..": "warning",
        "other": "major",
        "Benign": "major",
        "Benign/Likely_benign": "major",
        "Likely_benign": "major",
        "Uncertain_significance": "major",
        "Conflicting_interpretations_of_pathogenicity": "major",
        "Likely_pathogenic": "null",
        "Pathogenic/Likely_pathogenic": "null",
        "Pathogenic": "minor",
    },
    "Pathogenic": {
        "..": "warning",
        "other": "major",
        "Benign": "major",
        "Benign/Likely_benign": "major",
        "Likely_benign": "major",
        "Uncertain_significance": "major",
        "Conflicting_interpretations_of_pathogenicity": "major",
        "Likely_pathogenic": "minor",
        "Pathogenic/Likely_pathogenic": "minor",
        "Pathogenic": "null",
    },
}

CLINVAR_REVIEW_STATUS_RANK = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_interpretations": 0.5,
    "no_assertion_provided": 0,
    "no_assertion_criteria_provided": 0,
}

CLINVAR_PATHO_STATUS_RANK = {
    "Pathogenic": 5,
    "Pathogenic/Likely_pathogenic": 4.5,
    "Likely_pathogenic": 4,
}