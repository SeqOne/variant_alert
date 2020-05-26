# coding=utf-8
# =========================================================================
# Modules
# =========================================================================
import csv
import os
import re

import pysam

from const.clinvar import (
    ACMG_CLASS_BREAKING_CHANGE,
    CLINVAR_REVIEW_STATUS_RANK,
    CLINVAR_PATHO_STATUS_RANK,
    EXCEPTION_GENES,
)
from logger import logger
from parsers.vcf import VCFParser, MalformedVCFError


class ClinVarComparisonError(Exception):
    pass


class ClinvarVCFParser(VCFParser):
    def __init__(self, vcf_file_path):
        super().__init__(vcf_file_path)
        self.logger = logger.get_logger(__name__)

        self.vcf_type_name = "Clinvar VCF"
        self._get_fields_dict_methods = {
            **self._get_fields_dict_methods,
            **{
                "ALLELEID": self._get_variant_record_alleleid,
                "GENEINFO": self._get_variant_record_geneinfo,
                "CLNREVSTAT": self._get_variant_record_clnrevstat,
                "CLNRECSTAT": self._get_variant_record_clnrecstat,
                "CLNSIG": self._get_variant_record_clnsig,
                "OLD_CLNSIG": self._get_variant_record_old_clnsig,
                "MC": self._get_variant_record_mc,
                "RS": self._get_variant_record_rs,
            },
        }

        self._variant_buffer = []
        self._gene_buffer = []

    @staticmethod
    def _get_variant_record_alleleid(variant_record: pysam.VCFRecord):
        return {"alleleid": variant_record.info.get("ALLELEID")}

    @staticmethod
    def _get_variant_record_geneinfo(variant_record: pysam.VCFRecord):
        geneinfo_field = variant_record.info.get("GENEINFO")
        if geneinfo_field is None:
            geneinfo = "NA"
            geneinfo_id = "NA"
        else:
            # split on | to get only the first gene if there is many
            geneinfo, geneinfo_id = geneinfo_field.split("|")[0].split(":")
        return {"gene": geneinfo, "gene_id": geneinfo_id}

    @staticmethod
    def _get_variant_record_clnrevstat(variant_record: pysam.VCFRecord):
        return {"clnrevstat": ",".join(variant_record.info.get("CLNREVSTAT"))}

    @staticmethod
    def _get_variant_record_clnrecstat(variant_record: pysam.VCFRecord):
        return {"clnrecstat": variant_record.info.get("CLNRECSTAT")}

    @staticmethod
    def _get_variant_record_clnsig(variant_record: pysam.VCFRecord):
        if variant_record.info.get("CLNSIG") is None:
            return {"clnsig": ".."}
        else:
            return {"clnsig": variant_record.info.get("CLNSIG")[0]}

    @staticmethod
    def _get_variant_record_old_clnsig(variant_record: pysam.VCFRecord):
        return {"old_clnsig": variant_record.info.get("OLD_CLNSIG")}

    def _get_variant_record_mc(self, variant_record: pysam.VCFRecord):
        variant_record_mc = variant_record.info.get("MC")
        if variant_record_mc is None:
            return {"mc": []}
        try:
            variant_record_mc_values = self._split_piped_annotation(
                "MC", variant_record_mc
            )
        except AttributeError as e:
            raise MalformedVCFError(e)

        variant_record_mc_annotations = []
        for info_subvalues in variant_record_mc_values:
            so_id = info_subvalues[0].split(":")[-1]
            molecular_consequence = info_subvalues[1]
            variant_record_mc_annotations.append(
                {"so_id": so_id, "molecular_consequence": molecular_consequence}
            )
        return {"mc": variant_record_mc_annotations}

    @staticmethod
    def _get_variant_record_rs(variant_record: pysam.VCFRecord):
        return {"rs": variant_record.info.get("RS")}

    @property
    def variant_buffer(self):
        if not self._variant_buffer:
            self.set_variant_buffer(self.get_cln_variant_buffer())
        return self._variant_buffer

    def set_variant_buffer(self, value):
        self._variant_buffer = value

    @property
    def gene_buffer(self):
        if not self._gene_buffer:
            self.set_gene_buffer(self.get_pathogenic_gene_buffer())
        return self._gene_buffer

    def set_gene_buffer(self, value):
        self._gene_buffer = value

    @property
    def release(self):
        return self.file_date.strftime("%Y%m%d")

    def get_cln_variant_buffer(self) -> dict:
        """
        Put necessary variant record info into a dict to ease comparison with another vcf
        """
        return self.get_variant_buffer(
            [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "CLNSIG",
                "OLD_CLNSIG",
                "CLNREVSTAT",
                "CLNRECSTAT",
                "ALLELEID",
                "GENEINFO",
                "MC",
                "RS",
            ]
        )

    def get_pathogenic_gene_buffer(self):
        """
            Get pathogenic genes from put_variants_in_buffer with highest status
            (review confidence and ACMG classification)
            """
        variant_buffer = self.variant_buffer
        gene_buffer = {}
        for variant_id in variant_buffer:
            clnsig = variant_buffer[variant_id]["clnsig"]
            gene = variant_buffer[variant_id]["gene_id"]
            review = variant_buffer[variant_id]["clnrevstat"]
            if re.search(r"[pP]athogenic$", clnsig):
                gene_buffer.setdefault(
                    gene,
                    {
                        "best_review": {
                            "review": "",
                            "review_rank": -1,
                            "clnsig_rank": -1,
                            "clnsig": "",
                        },
                        "best_clnsig": {
                            "review": "",
                            "review_rank": -1,
                            "clnsig_rank": -1,
                            "clnsig": "",
                        },
                        "gene": variant_buffer[variant_id]["gene"],
                    },
                )
                # find the highest review confidence variant
                # if a better review status variant is detected, keep
                if (
                    gene_buffer[gene]["best_review"]["review_rank"]
                    < CLINVAR_REVIEW_STATUS_RANK[review]
                ):
                    gene_buffer[gene]["best_review"] = {
                        "review": review,
                        "review_rank": CLINVAR_REVIEW_STATUS_RANK[review],
                        "clnsig_rank": CLINVAR_PATHO_STATUS_RANK[clnsig],
                        "clnsig": clnsig,
                    }
                # if an equal better review status variant is detected, keep the highest pathogenicity classication
                elif (
                    gene_buffer[gene]["best_review"]["review_rank"]
                    == CLINVAR_REVIEW_STATUS_RANK[review]
                ):
                    if (
                        gene_buffer[gene]["best_review"]["clnsig_rank"]
                        < CLINVAR_PATHO_STATUS_RANK[clnsig]
                    ):
                        gene_buffer[gene]["best_review"] = {
                            "review": review,
                            "review_rank": CLINVAR_REVIEW_STATUS_RANK[review],
                            "clnsig_rank": CLINVAR_PATHO_STATUS_RANK[clnsig],
                            "clnsig": clnsig,
                        }
                # find the highest pathogenic variant
                # if a better pathogenic classification is detected, keep
                if (
                    gene_buffer[gene]["best_clnsig"]["clnsig_rank"]
                    < CLINVAR_PATHO_STATUS_RANK[clnsig]
                ):
                    gene_buffer[gene]["best_clnsig"] = {
                        "review": review,
                        "review_rank": CLINVAR_REVIEW_STATUS_RANK[review],
                        "clnsig_rank": CLINVAR_PATHO_STATUS_RANK[clnsig],
                        "clnsig": clnsig,
                    }
                # if an equal pathogenic classification is detected, keep the highest review status variant
                elif (
                    gene_buffer[gene]["best_clnsig"]["clnsig_rank"]
                    == CLINVAR_PATHO_STATUS_RANK[clnsig]
                ):
                    if (
                        gene_buffer[gene]["best_clnsig"]["review_rank"]
                        < CLINVAR_REVIEW_STATUS_RANK[review]
                    ):
                        gene_buffer[gene]["best_clnsig"] = {
                            "review": review,
                            "review_rank": CLINVAR_REVIEW_STATUS_RANK[review],
                            "clnsig_rank": CLINVAR_PATHO_STATUS_RANK[clnsig],
                            "clnsig": clnsig,
                        }
        return gene_buffer


class ClinvarVCFComparator(object):
    """
    Comparison object used to compare two ClinVar VCF file
    """

    @property
    def source_vcf(self):
        return self._source_vcf

    @property
    def target_vcf(self):
        return self._target_vcf

    def __init__(self, source_vcf: str, target_vcf: str):
        self._source_vcf = ClinvarVCFParser(source_vcf)
        self._target_vcf = ClinvarVCFParser(target_vcf)

        self.logger = logger.get_logger(__name__)

    def compare_references(self):
        """
        Uses VCFs to compare their genome reference, raise an error if they are not the same
        """
        if self.source_vcf.reference != self.target_vcf.reference:
            self.logger.error("Not the same genome references")
            raise ClinVarComparisonError(
                "Genome reference between the two VCF file is not the same: \n"
                f"source: {self.source_vcf.reference}\n"
                f"target: {self.target_vcf.reference}"
            )
        else:
            self.logger.info("The genome reference is the same across the VCF files")
            return self.source_vcf.reference

    def compare_variants(self):
        """
        Uses VCFs variant buffer and compare them to detect all changes and output a list of variant dict
        """
        source_variant_buffer = self.source_vcf.variant_buffer.copy()
        target_variant_buffer = self.target_vcf.variant_buffer.copy()

        compared_variants = []
        variants_have_been_lost = False
        self.logger.info("Comparing source and target vcf variants")
        # return variants from the latest vcf
        for variant_id in target_variant_buffer:
            target_variant_clnsig = target_variant_buffer[variant_id]["clnsig"]
            target_variant_clinvar_id = target_variant_buffer[variant_id]["id"]
            target_variant_clnrevstat = target_variant_buffer[variant_id]["clnrevstat"]
            target_variant_clnrecstat = target_variant_buffer[variant_id]["clnrecstat"]
            target_variant_gene_info = target_variant_buffer[variant_id]["gene"]
            target_variant_geneinfo_id = target_variant_buffer[variant_id]["gene_id"]
            # return variants with important change of classification
            if variant_id in source_variant_buffer:
                source_variant_clnsig = source_variant_buffer[variant_id]["clnsig"]
                if source_variant_clnsig not in ACMG_CLASS_BREAKING_CHANGE.keys():
                    source_variant_clnsig = "other"

                if (
                    target_variant_clnsig
                    not in ACMG_CLASS_BREAKING_CHANGE[source_variant_clnsig].keys()
                ):
                    target_variant_clnsig = "other"

                breaking_change = ACMG_CLASS_BREAKING_CHANGE[source_variant_clnsig][
                    target_variant_clnsig
                ]

                if target_variant_clnsig != source_variant_clnsig:
                    compared_variants.append(
                        {
                            "variant_id": variant_id,
                            "clinvar_id": target_variant_clinvar_id,
                            "old_classification": source_variant_clnsig,
                            "new_classification": target_variant_clnsig,
                            "breaking_change": breaking_change,
                            "confidence": target_variant_clnrevstat,
                            "reclassification_status": target_variant_clnrecstat,
                            "gene_info": target_variant_gene_info,
                            "gene_info_id": target_variant_geneinfo_id,
                            "name_clinvar_old": self.source_vcf.release,
                            "name_clinvar_new": self.target_vcf.release,
                        }
                    )

                del source_variant_buffer[variant_id]

            # return new variants that did not exist in older version
            else:
                if target_variant_clnsig not in ACMG_CLASS_BREAKING_CHANGE.keys():
                    target_variant_clnsig = "other"

                breaking_change = ACMG_CLASS_BREAKING_CHANGE[".."][
                    target_variant_clnsig
                ]

                compared_variants.append(
                    {
                        "variant_id": variant_id,
                        "clinvar_id": target_variant_clinvar_id,
                        "old_classification": "..",
                        "new_classification": target_variant_clnsig,
                        "breaking_change": breaking_change,
                        "confidence": target_variant_clnrevstat,
                        "reclassification_status": target_variant_clnrecstat,
                        "gene_info": target_variant_gene_info,
                        "gene_info_id": target_variant_geneinfo_id,
                        "name_clinvar_old": self.source_vcf.release,
                        "name_clinvar_new": self.target_vcf.release,
                    }
                )
        # return variants present in old vcf that are not available in latest vcf
        for variant_id in source_variant_buffer:
            variants_have_been_lost = True
            source_variant_clnsig = source_variant_buffer[variant_id]["clnsig"]
            source_variant_clinvar_id = source_variant_buffer[variant_id]["id"]
            source_variant_clnrevstat = source_variant_buffer[variant_id]["clnrevstat"]
            source_variant_clnrecstat = source_variant_buffer[variant_id]["clnrecstat"]
            source_variant_gene_info = source_variant_buffer[variant_id]["gene"]
            source_variant_geneinfo_id = source_variant_buffer[variant_id]["gene_id"]

            if source_variant_clnsig not in ACMG_CLASS_BREAKING_CHANGE.keys():
                source_variant_clnsig = "other"

            breaking_change = ACMG_CLASS_BREAKING_CHANGE[source_variant_clnsig][".."]

            compared_variants.append(
                {
                    "variant_id": variant_id,
                    "clinvar_id": source_variant_clinvar_id,
                    "old_classification": source_variant_clnsig,
                    "new_classification": "..",
                    "breaking_change": breaking_change,
                    "confidence": source_variant_clnrevstat,
                    "reclassification_status": source_variant_clnrecstat,
                    "gene_info": source_variant_gene_info,
                    "gene_info_id": source_variant_geneinfo_id,
                    "name_clinvar_old": self.source_vcf.release,
                    "name_clinvar_new": self.target_vcf.release,
                }
            )
        self.logger.info("Comparison done")
        return compared_variants, variants_have_been_lost

    def compare_genes(self):
        """
        Compare VCFs gene buffer to detect all changes between them and output a list of genes with
        annotations and change status. In more biological terms it gives :
            - new genes with pathogenic variants,
            - update on the highest review status classification,
            - update on the highest ACMG classfication
        """

        self.logger.info("Comparing source and target vcf genes")
        compared_genes = []
        # create pathogenic clinvar gene dictionary
        source_pathogenic_genes = self.source_vcf.gene_buffer.copy()
        target_pathogenic_genes = self.target_vcf.gene_buffer.copy()

        # append lost pathogenic genes between 2 versions
        lost_pathogenic_gene_set = (
            source_pathogenic_genes.keys() - target_pathogenic_genes.keys()
        )
        for gene_key in lost_pathogenic_gene_set:
            compared_genes.append(
                {
                    "gene_info_id": gene_key,
                    "gene_info": source_pathogenic_genes[gene_key]["gene"],
                    "pathogenic_class_status": "LOST_PATHOGENICITY",
                    "pathogenic_class_old": source_pathogenic_genes[gene_key][
                        "best_clnsig"
                    ]["clnsig"],
                    "pathogenic_class_best_review_confidence_old": source_pathogenic_genes[
                        gene_key
                    ][
                        "best_clnsig"
                    ][
                        "review"
                    ],
                    "pathogenic_class_new": "..",
                    "pathogenic_class_best_review_confidence_new": "..",
                    "review_confidence_status": "NO_REVIEW",
                    "review_confidence_best_pathogenic_class_old": source_pathogenic_genes[
                        gene_key
                    ][
                        "best_review"
                    ][
                        "clnsig"
                    ],
                    "review_confidence_old": source_pathogenic_genes[gene_key][
                        "best_review"
                    ]["review"],
                    "review_confidence_best_pathogenic_class_new": "..",
                    "review_confidence_new": "..",
                    "name_clinvar_old": self.source_vcf.release,
                    "name_clinvar_new": self.target_vcf.release,
                }
            )

        # append new pathogenic genes between 2 versions
        new_pathogenic_gene_set = (
            target_pathogenic_genes.keys() - source_pathogenic_genes.keys()
        )
        for gene_key in new_pathogenic_gene_set:
            compared_genes.append(
                {
                    "gene_info_id": gene_key,
                    "gene_info": target_pathogenic_genes[gene_key]["gene"],
                    "pathogenic_class_status": "NEW_PATHOGENICITY",
                    "pathogenic_class_old": "..",
                    "pathogenic_class_best_review_confidence_old": "..",
                    "pathogenic_class_new": target_pathogenic_genes[gene_key][
                        "best_clnsig"
                    ]["clnsig"],
                    "pathogenic_class_best_review_confidence_new": target_pathogenic_genes[
                        gene_key
                    ][
                        "best_clnsig"
                    ][
                        "review"
                    ],
                    "review_confidence_status": "NEW_REVIEW",
                    "review_confidence_best_pathogenic_class_old": "..",
                    "review_confidence_old": "..",
                    "review_confidence_best_pathogenic_class_new": target_pathogenic_genes[
                        gene_key
                    ][
                        "best_review"
                    ][
                        "clnsig"
                    ],
                    "review_confidence_new": target_pathogenic_genes[gene_key][
                        "best_review"
                    ]["review"],
                    "name_clinvar_old": self.source_vcf.release,
                    "name_clinvar_new": self.target_vcf.release,
                }
            )

        # Create a dictionary of common pathogenic genes
        # Apply conditions to select for each genes :
        # - highest ACMG classfication status (and associated review confidence)
        # - highest review confidence status (and associated ACMG classifcation)
        # Report review confidence status and ACMG classfication status changes
        common_pathogenic_gene_set = (
            source_pathogenic_genes.keys() & target_pathogenic_genes.keys()
        )
        for gene_key in common_pathogenic_gene_set:
            geneinfo = target_pathogenic_genes[gene_key]["gene"]
            pathogenic_class_status = "UNCHANGED"
            pathogenic_class_old = source_pathogenic_genes[gene_key]["best_clnsig"][
                "clnsig"
            ]
            pathogenic_class_best_review_confidence_old = source_pathogenic_genes[
                gene_key
            ]["best_clnsig"]["review"]
            pathogenic_class_new = target_pathogenic_genes[gene_key]["best_clnsig"][
                "clnsig"
            ]
            pathogenic_class_best_review_confidence_new = target_pathogenic_genes[
                gene_key
            ]["best_clnsig"]["review"]
            review_confidence_status = "UNCHANGED"
            review_confidence_best_pathogenic_class_old = source_pathogenic_genes[
                gene_key
            ]["best_review"]["clnsig"]
            review_confidence_old = source_pathogenic_genes[gene_key]["best_clnsig"][
                "review"
            ]
            review_confidence_best_pathogenic_class_new = target_pathogenic_genes[
                gene_key
            ]["best_review"]["clnsig"]
            review_confidence_new = target_pathogenic_genes[gene_key]["best_review"][
                "review"
            ]

            # Select the highest ACMG classification status
            if (
                source_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
                == target_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
            ):
                pass
            elif (
                source_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
                > target_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
            ):
                pathogenic_class_status = "DOWNGRADED_PATHOGENICITY_STATUS"
            elif (
                source_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
                < target_pathogenic_genes[gene_key]["best_clnsig"]["clnsig_rank"]
            ):
                pathogenic_class_status = "UPGRADED_PATHOGENICITY_STATUS"

            # Select the highest review confidence status
            if (
                source_pathogenic_genes[gene_key]["best_review"]["review_rank"]
                == target_pathogenic_genes[gene_key]["best_review"]["review_rank"]
            ):
                pass
            elif (
                source_pathogenic_genes[gene_key]["best_review"]["review_rank"]
                > target_pathogenic_genes[gene_key]["best_review"]["review_rank"]
            ):
                review_confidence_status = "DOWNGRADED_REVIEW_CONFIDENCE"
            elif (
                source_pathogenic_genes[gene_key]["best_review"]["review_rank"]
                < target_pathogenic_genes[gene_key]["best_review"]["review_rank"]
            ):
                review_confidence_status = "UPGRADED_REVIEW_CONFIDENCE"
            else:
                raise ValueError(
                    "This error is not supposed to happen. Something has probably changed in your "
                    "data types, please check"
                )
            if (
                review_confidence_status == "UNCHANGED"
                and pathogenic_class_status == "UNCHANGED"
            ):
                continue
            else:
                compared_genes.append(
                    {
                        "gene_info_id": gene_key,
                        "gene_info": geneinfo,
                        "pathogenic_class_status": pathogenic_class_status,
                        "pathogenic_class_old": pathogenic_class_old,
                        "pathogenic_class_best_review_confidence_old": pathogenic_class_best_review_confidence_old,
                        "pathogenic_class_new": pathogenic_class_new,
                        "pathogenic_class_best_review_confidence_new": pathogenic_class_best_review_confidence_new,
                        "review_confidence_status": review_confidence_status,
                        "review_confidence_best_pathogenic_class_old": review_confidence_best_pathogenic_class_old,
                        "review_confidence_old": review_confidence_old,
                        "review_confidence_best_pathogenic_class_new": review_confidence_best_pathogenic_class_new,
                        "review_confidence_new": review_confidence_new,
                        "name_clinvar_old": self.source_vcf.release,
                        "name_clinvar_new": self.target_vcf.release,
                    }
                )
        # This is a hot fix to deal with exception genes (until clinvcf deal with gene extraction from clinvar XML release)
        for item in compared_genes:
            if item["gene_info_id"] in EXCEPTION_GENES.keys():
                item["gene_info_id"] = EXCEPTION_GENES[item["gene_info_id"]]
                item["gene_info"] = EXCEPTION_GENES[item["gene_info"]]
            else:
                pass

        return compared_genes

    def write_variant_comparison(
        self,
        variant_comparison_result: list,
        output_directory: str = os.getcwd(),
        output_format: str = "tsv",
    ):
        variant_comparison_output_formatters = {
            "tsv": self.write_variant_comparison_as_tsv,
            "vcf": self.write_variant_comparison_as_vcf,
        }
        try:
            variant_comparison_output_formatters[output_format](
                variant_comparison_result, output_directory
            )
        except KeyError:
            self.logger.exception(
                "Selected output formatter is not available please choose between:"
                f" {', '.join(variant_comparison_output_formatters.keys())}"
            )
            raise

    def write_variant_comparison_as_tsv(
        self, variant_comparison_result: list, output_directory: str
    ):
        headers = [
            "variant_id",
            "clinvar_id",
            "old_classification",
            "new_classification",
            "breaking_change",
            "confidence",
            "reclassification_status",
            "gene_info",
            "gene_info_id",
            "name_clinvar_old",
            "name_clinvar_new",
        ]
        output = os.path.join(
            output_directory,
            f"clinvar_variant_diff_from_{self.source_vcf.release}_to_{self.target_vcf.release}.tsv",
        )
        csv.register_dialect("tsv", delimiter="\t")
        self.logger.info("Opening output variant comparison file")
        with open(output, "w", newline="") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=headers, dialect="tsv")

            writer.writeheader()
            for variant in variant_comparison_result:
                writer.writerow(variant)

        self.logger.info(f"variant comparison has been written to file {output}")

    def write_variant_comparison_as_vcf(
        self, variant_comparison_result: list, output_directory: str = os.getcwd()
    ):
        # @formatter:off
        headers = [
            '##INFO=<ID=CLNSRCSIG,Number=.,Type=String,Description="Clinical significance for this single variant in source version">',  # noqa E501
            '##INFO=<ID=CLNTGTSIG,Number=.,Type=String,Description="Clinical significance for this single variant in target version">',  # noqa E501
            '##INFO=<ID=CLNBRKCHANGE,Number=.,Type=String,Description="Clinical significance change level between source and target clinvar version">',  # noqa E501
            '##INFO=<ID=CLNGENE,Number=.,Type=String,Description="Retained gene to describe clinical significance">',  # noqa E501
            '##INFO=<ID=CLNGENEID,Number=.,Type=String,Description="Retained gened id to describe clinical significance">',  # noqa E501
            f"##clinvar_source_release={self.source_vcf.release}",
            f"##clinvar_target_release={self.target_vcf.release}",
        ]
        # @formatter:on

        output = os.path.join(
            output_directory,
            f"clinvar_variant_diff_from_{self.source_vcf.release}_to_{self.target_vcf.release}.vcf",
        )

        output_header = self.target_vcf.vcf_file.header
        source_header = (
            self.source_vcf.vcf_file.header
        )  # Only necessary to use variant record object properly

        for header in headers:
            output_header.add_line(header)
            source_header.add_line(header)

        output_vcf = pysam.VariantFile(output, mode="w", header=output_header)

        variant_record_dict = {}
        for variant_record in self.target_vcf.vcf_file.fetch():
            variant_id = self.target_vcf.set_genome_variant_id(
                self.target_vcf.reference,
                variant_record.chrom,
                variant_record.pos,
                variant_record.ref,
                variant_record.alts[0],
            )
            variant_record_dict.setdefault(variant_id, variant_record.copy())

        for variant_record in self.source_vcf.vcf_file.fetch():
            variant_id = self.source_vcf.set_genome_variant_id(
                self.source_vcf.reference,
                variant_record.chrom,
                variant_record.pos,
                variant_record.ref,
                variant_record.alts[0],
            )
            variant_record_dict.setdefault(variant_id, variant_record.copy())
        for compared_variant in variant_comparison_result:
            variant_record = variant_record_dict.get(
                compared_variant["variant_id"], None
            )
            if variant_record:
                variant_record.info.update(
                    {
                        "CLNSRCSIG": compared_variant["old_classification"],
                        "CLNTGTSIG": compared_variant["new_classification"],
                        "CLNBRKCHANGE": compared_variant["breaking_change"],
                        "CLNGENE": compared_variant["gene_info"],
                        "CLNGENEID": compared_variant["gene_info_id"],
                    }
                )
                output_vcf.write(variant_record)
        output_vcf.close()

    def write_gene_comparison(
        self, gene_comparison_result: list, output_directory: str = os.getcwd()
    ):
        headers = [
            "gene_info_id",
            "gene_info",
            "pathogenic_class_status",
            "pathogenic_class_old",
            "pathogenic_class_best_review_confidence_old",
            "pathogenic_class_new",
            "pathogenic_class_best_review_confidence_new",
            "review_confidence_status",
            "review_confidence_best_pathogenic_class_old",
            "review_confidence_old",
            "review_confidence_best_pathogenic_class_new",
            "review_confidence_new",
            "name_clinvar_old",
            "name_clinvar_new",
        ]
        # Create the gene_diff.tsv file, with changes in gene status between 2 ClinVar versions
        output = os.path.join(
            output_directory,
            f"clinvar_gene_diff_from_{self.source_vcf.release}_to_{self.target_vcf.release}.tsv",
        )
        csv.register_dialect("tsv", delimiter="\t")
        self.logger.info("Opening output gene comparison file")
        with open(output, "w", newline="") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=headers, dialect="tsv")

            writer.writeheader()
            for gene in gene_comparison_result:
                writer.writerow(gene)

        self.logger.info(f"gene comparison has been written to file {output}")

    def exception_gene_check(self, gene_key: str):
        # This is a hot fix to deal with exception genes (until clinvcf deal with gene extraction from clinvar XML release)
        if gene_key in EXCEPTION_GENES.keys():
            return EXCEPTION_GENES[gene_key]
        else:
            return gene_key

    def write_clinvarome(self, output_directory: str = os.getcwd()):
        headers = [
            "gene_info_id",
            "gene_info",
            "pathogenic_class_new",
            "pathogenic_class_best_review_confidence_new",
            "review_confidence_best_pathogenic_class_new",
            "review_confidence_new",
            "name_clinvar_new",
        ]

        output = os.path.join(
            output_directory, f"clinvarome_{self.target_vcf.release}.tsv"
        )
        csv.register_dialect("tsv", delimiter="\t")
        self.logger.info("Opening output clinvarome file")
        with open(output, "w", newline="") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=headers, dialect="tsv")
            writer.writeheader()

            pathogenic_genes = self.target_vcf.gene_buffer

            for gene_key in pathogenic_genes:
                writer.writerow(
                    {
                        "gene_info_id": self.exception_gene_check(gene_key),
                        "gene_info": self.exception_gene_check(
                            pathogenic_genes[gene_key]["gene"]
                        ),
                        "pathogenic_class_new": pathogenic_genes[gene_key][
                            "best_clnsig"
                        ]["clnsig"],
                        "pathogenic_class_best_review_confidence_new": pathogenic_genes[
                            gene_key
                        ]["best_clnsig"]["review"],
                        "review_confidence_best_pathogenic_class_new": pathogenic_genes[
                            gene_key
                        ]["best_review"]["clnsig"],
                        "review_confidence_new": pathogenic_genes[gene_key][
                            "best_review"
                        ]["review"],
                        "name_clinvar_new": self.target_vcf.release,
                    }
                )
        self.logger.info(f"clinvarome has been written to file {output}")
