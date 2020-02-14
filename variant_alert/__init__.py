# coding=utf-8
# =========================================================================
# Modules
# =========================================================================


class VariantAlertParams(object):
    """
    Params shared through SamQ group commands
    """

    def __init__(self):
        self.clinvar_vcf_comparator = ""
        self.output_directory = ""
