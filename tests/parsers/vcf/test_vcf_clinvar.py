import os
import unittest

from parsers.vcf.clinvar import (
    ClinvarVCFParser,
    ClinvarVCFComparator,
    ClinVarComparisonError,
)
from tests import TEST_DIR


class ClinvarVCFParserTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_20191202.vcf.gz"
        )
        self.clinvar_vcf_parser = ClinvarVCFParser(self.clinvar_vcf)

    def test_clinvar_vcf_parser_init(self):
        self.assertEqual("Clinvar VCF", self.clinvar_vcf_parser.vcf_type_name)

        self.assertEqual(
            "2019-12-02", self.clinvar_vcf_parser.file_date.strftime("%Y-%m-%d")
        )
        self.assertEqual("20191202", self.clinvar_vcf_parser.release)
        self.assertIsNotNone(self.clinvar_vcf_parser.variant_buffer)
        self.assertIsNotNone(self.clinvar_vcf_parser.gene_buffer)

    def test_clinvar_vcf_parser_variant_buffer(self):
        self.assertEqual(71, len(self.clinvar_vcf_parser.variant_buffer))
        variants = []
        for variant_record in self.clinvar_vcf_parser.vcf_file.fetch():
            variants.append(
                self.clinvar_vcf_parser.set_genome_variant_id(
                    self.clinvar_vcf_parser.reference,
                    variant_record.chrom,
                    variant_record.pos,
                    variant_record.ref,
                    variant_record.alts[0],
                )
            )
        self.assertCountEqual(variants, self.clinvar_vcf_parser.variant_buffer.keys())

    def test_clinvar_vcf_parser_gene_buffer(self):
        self.assertEqual(2, len(self.clinvar_vcf_parser.gene_buffer))
        expected_gene_buffer = {
            "9636": {
                "best_conf": {
                    "conf": "no_assertion_criteria_provided",
                    "conf_rank": 0,
                    "patho_rank": 5,
                    "clnsig": "Pathogenic",
                },
                "best_clnsig": {
                    "conf": "no_assertion_criteria_provided",
                    "conf_rank": 0,
                    "patho_rank": 5,
                    "clnsig": "Pathogenic",
                },
                "geneinfo": "ISG15",
            },
            "375790": {
                "best_conf": {
                    "conf": "no_assertion_criteria_provided",
                    "conf_rank": 0,
                    "patho_rank": 5,
                    "clnsig": "Pathogenic",
                },
                "best_clnsig": {
                    "conf": "no_assertion_criteria_provided",
                    "conf_rank": 0,
                    "patho_rank": 5,
                    "clnsig": "Pathogenic",
                },
                "geneinfo": "AGRN",
            },
        }
        self.assertCountEqual(expected_gene_buffer, self.clinvar_vcf_parser.gene_buffer)


class ClinvarVCFComparatorTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_target.vcf.gz"
        )
        cls.source_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_source.vcf.gz"
        )
        cls.target_clinvar_vcf_with_missing_variants = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_target_with_lost_variants.vcf.gz"
        )
        cls.hg38_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/vcf_new_38.vcf.gz"
        )

    def test_compare_reference(self):
        source_clinvar_vcf_parser = ClinvarVCFParser(self.source_clinvar_vcf)
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        self.assertEqual(
            source_clinvar_vcf_parser.reference,
            clinvar_vcf_comparator.compare_references(),
        )

        clinvar_vcf_comparator = ClinvarVCFComparator(
            self.source_clinvar_vcf, self.hg38_clinvar_vcf
        )
        with self.assertRaises(ClinVarComparisonError):
            clinvar_vcf_comparator.compare_references()

    def test_compare_variants(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        expected_compared_variants = [
            {
                "variant_id": "GRCh37_1_949422_G_A",
                "clinvar_id": "475282",
                "old_classification": "..",
                "new_classification": "Pathogenic",
                "breaking_change": "major",
                "confidence": "criteria_provided,_single_submitter",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "variant_id": "GRCh37_3_37056026_C_T",
                "clinvar_id": "183381",
                "old_classification": "Pathogenic",
                "new_classification": "Benign",
                "breaking_change": "major",
                "confidence": "no_assertion_criteria_provided",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "variant_id": "GRCh37_3_50416389_C_T",
                "clinvar_id": "542075",
                "old_classification": "Uncertain_significance",
                "new_classification": "Likely_pathogenic",
                "breaking_change": "major",
                "confidence": "criteria_provided,_single_submitter",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "variant_id": "GRCh37_6_150655494_C_T",
                "clinvar_id": "475278",
                "old_classification": "Benign",
                "new_classification": "Uncertain_significance",
                "breaking_change": "unknown",
                "confidence": "criteria_provided,_single_submitter",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "variant_id": "GRCh37_4_84383745_T_TC",
                "clinvar_id": "128215",
                "old_classification": "Benign",
                "new_classification": "Pathogenic",
                "breaking_change": "major",
                "confidence": "criteria_provided,_multiple_submitters,_no_conflicts",
                "gene_info": "ABRAXAS1",
                "gene_info_id": "84142",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
        ]

        clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()
        self.assertFalse(variants_have_been_lost)
        self.assertCountEqual(expected_compared_variants, compared_variants)

        # With missing variants
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf,
            target_vcf=self.target_clinvar_vcf_with_missing_variants,
        )
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()
        self.assertTrue(variants_have_been_lost)

    def test_compare_gene(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        expected_compared_genes = [
            {
                "gene_info_id": "84142",
                "gene_info": "ABRAXAS1",
                "pathogenic_class_status": "NEW_PATHOGENICITY",
                "pathogenic_class_old": "..",
                "pathogenic_class_best_review_confidence_old": "..",
                "pathogenic_class_new": "Pathogenic",
                "pathogenic_class_best_review_confidence_new": "criteria_provided,_multiple_submitters,_no_conflicts",
                "review_confidence_status": "NEW_REVIEW",
                "review_confidence_best_pathogenic_class_old": "..",
                "review_confidence_old": "..",
                "review_confidence_best_pathogenic_class_new": "Pathogenic",
                "review_confidence_new": "criteria_provided,_multiple_submitters,_no_conflicts",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "gene_info_id": "9636",
                "gene_info": "ISG15",
                "pathogenic_class_status": "UNCHANGED",
                "pathogenic_class_old": "Pathogenic",
                "pathogenic_class_best_review_confidence_old": "no_assertion_criteria_provided",
                "pathogenic_class_new": "Pathogenic",
                "pathogenic_class_best_review_confidence_new": "criteria_provided,_single_submitter",
                "review_confidence_status": "UPGRADED_REVIEW_CONFIDENCE",
                "review_confidence_best_pathogenic_class_old": "Pathogenic",
                "review_confidence_old": "no_assertion_criteria_provided",
                "review_confidence_best_pathogenic_class_new": "Pathogenic",
                "review_confidence_new": "criteria_provided,_single_submitter",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
        ]

        clinvar_vcf_comparator.compare_references()

        compared_genes = clinvar_vcf_comparator.compare_genes()
        self.assertCountEqual(expected_compared_genes, compared_genes)

    def test_write_variant_comparison(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")
        output_format = "tsv"
        # @formatters:off
        expected_comparison_output_content = [
            "variant_id\tclinvar_id\told_classification\tnew_classification\tbreaking_change\tconfidence\tgene_info\tgene_info_id\tname_clinvar_old\tname_clinvar_new\n",  # noqa E501
            "GRCh37_1_949422_G_A\t475282\t..\tPathogenic\tmajor\tcriteria_provided,_single_submitter\tISG15\t9636\t20190805\t20190915\n",  # noqa E501
            "GRCh37_3_37056026_C_T\t183381\tPathogenic\tBenign\tmajor\tno_assertion_criteria_provided\tISG15\t9636\t20190805\t20190915\n",  # noqa E501
            "GRCh37_3_50416389_C_T\t542075\tUncertain_significance\tLikely_pathogenic\tmajor\tcriteria_provided,_single_submitter\tISG15\t9636\t20190805\t20190915\n",  # noqa E501
            "GRCh37_6_150655494_C_T\t475278\tBenign\tUncertain_significance\tunknown\tcriteria_provided,_single_submitter\tISG15\t9636\t20190805\t20190915\n",  # noqa E501
            "GRCh37_4_84383745_T_TC\t128215\tBenign\tPathogenic\tmajor\tcriteria_provided,_multiple_submitters,_no_conflicts\tABRAXAS1\t84142\t20190805\t20190915\n",  # noqa E501
            # "GRCh37_variant_X_153296290_G_A\t571208\tUncertain_significance\t..\twarning\tcriteria_provided,_single_submitter\tISG15\t9636\t20190805\t20190915\n",
            # # noqa E501
        ]
        # @formatters:on

        clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()

        clinvar_vcf_comparator.write_variant_comparison(
            compared_variants, output_directory, output_format
        )
        comparison_output = os.path.join(
            output_directory,
            (
                f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{clinvar_vcf_comparator.target_vcf.release}.{output_format}"
            ),
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            self.assertCountEqual(expected_comparison_output_content, co.readlines())

    def test_write_gene_comparison(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")
        # @formatters:off
        expected_comparison_output_content = [
            "gene_info_id\tgene_info\tpathogenic_class_status\tpathogenic_class_old\tpathogenic_class_best_review_confidence_old\tpathogenic_class_new\tpathogenic_class_best_review_confidence_new\treview_confidence_status\treview_confidence_best_pathogenic_class_old\treview_confidence_old\treview_confidence_best_pathogenic_class_new\treview_confidence_new\tname_clinvar_old\tname_clinvar_new\n",  # noqa E501
            "84142\tABRAXAS1\tNEW_PATHOGENICITY\t..\t..\tPathogenic\tcriteria_provided,_multiple_submitters,_no_conflicts\tNEW_REVIEW\t..\t..\tPathogenic\tcriteria_provided,_multiple_submitters,_no_conflicts\t20190805\t20190915\n",  # noqa E501
            "9636\tISG15\tUNCHANGED\tPathogenic\tno_assertion_criteria_provided\tPathogenic\tcriteria_provided,_single_submitter\tUPGRADED_REVIEW_CONFIDENCE\tPathogenic\tno_assertion_criteria_provided\tPathogenic\tcriteria_provided,_single_submitter\t20190805\t20190915\n",  # noqa E501
        ]
        # @formatters:on

        clinvar_vcf_comparator.compare_references()
        compared_genes = clinvar_vcf_comparator.compare_genes()
        clinvar_vcf_comparator.write_gene_comparison(compared_genes, output_directory)
        comparison_output = os.path.join(
            output_directory,
            (
                f"clinvar_gene_diff_from_{clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{clinvar_vcf_comparator.target_vcf.release}.tsv"
            ),
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            self.assertCountEqual(expected_comparison_output_content, co.readlines())

    def test_write_clinvarome(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")
        # @formatters:off
        expected_clinvarome_output_content = [
            "gene_info_id\tgene_info\tpathogenic_class_new\tpathogenic_class_best_review_confidence_new\treview_confidence_best_pathogenic_class_new\treview_confidence_new\tname_clinvar_new\n",  # noqa E501
            "9636\tISG15\tPathogenic\tcriteria_provided,_single_submitter\tPathogenic\tcriteria_provided,_single_submitter\t20190915\n",  # noqa E501
            "84142\tABRAXAS1\tPathogenic\tcriteria_provided,_multiple_submitters,_no_conflicts\tPathogenic\tcriteria_provided,_multiple_submitters,_no_conflicts\t20190915\n",  # noqa E501
        ]
        # @formatters:on
        clinvar_vcf_comparator.compare_references()
        clinvar_vcf_comparator.write_clinvarome(output_directory)
        comparison_output = os.path.join(
            output_directory,
            f"clinvarome_{clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            self.assertCountEqual(expected_clinvarome_output_content, co.readlines())

    def test_write_variant_comparison_as_vcf(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")

        clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()
        clinvar_vcf_comparator.write_variant_comparison_as_vcf(
            compared_variants, output_directory
        )
        comparison_output = os.path.join(
            output_directory,
            f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}_to_"
            f"{clinvar_vcf_comparator.target_vcf.release}.vcf",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    TEST_DIR,
                    "parsers/vcf/data/expected_clinvar_variant_diff_from_20190805_to_20190915.vcf",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    @classmethod
    def tearDownClass(cls) -> None:
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=cls.source_clinvar_vcf, target_vcf=cls.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")

        tsv_variant_comparison_output = os.path.join(
            output_directory,
            (
                f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{clinvar_vcf_comparator.target_vcf.release}.tsv"
            ),
        )
        if os.path.exists(tsv_variant_comparison_output):
            os.remove(tsv_variant_comparison_output)

        vcf_variant_comparison_output = os.path.join(
            output_directory,
            (
                f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{clinvar_vcf_comparator.target_vcf.release}.vcf"
            ),
        )
        if os.path.exists(vcf_variant_comparison_output):
            os.remove(vcf_variant_comparison_output)

        tsv_gene_comparison_output = os.path.join(
            output_directory,
            (
                f"clinvar_gene_diff_from_{clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{clinvar_vcf_comparator.target_vcf.release}.tsv"
            ),
        )
        if os.path.exists(tsv_gene_comparison_output):
            os.remove(tsv_gene_comparison_output)

        tsv_clinvarome_output = os.path.join(
            output_directory,
            f"clinvarome_{clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        if os.path.exists(tsv_clinvarome_output):
            os.remove(tsv_clinvarome_output)
