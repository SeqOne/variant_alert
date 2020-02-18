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
        self.vcf = os.path.join(TEST_DIR, "parsers/vcf/data/clinvar_20191202.vcf.gz")
        self.parser = ClinvarVCFParser(self.vcf)

    def test_init(self):
        self.assertEqual("Clinvar VCF", self.parser.vcf_type_name)

        self.assertEqual("2019-12-02", self.parser.file_date.strftime("%Y-%m-%d"))
        self.assertEqual("20191202", self.parser.release)
        self.assertIsNotNone(self.parser.variant_buffer)
        self.assertIsNotNone(self.parser.gene_buffer)

    def test_clinvar_vcf_parser_variant_buffer(self):
        self.assertEqual(71, len(self.parser.variant_buffer))
        variants = []
        for variant_record in self.parser.vcf_file.fetch():
            variants.append(
                self.parser.set_genome_variant_id(
                    self.parser.reference,
                    variant_record.chrom,
                    variant_record.pos,
                    variant_record.ref,
                    variant_record.alts[0],
                )
            )
        self.assertCountEqual(variants, self.parser.variant_buffer.keys())

    def test_clinvar_vcf_parser_gene_buffer(self):
        self.assertEqual(2, len(self.parser.gene_buffer))
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
        self.assertCountEqual(expected_gene_buffer, self.parser.gene_buffer)


class ClinvarVCFComparatorTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_target.vcf.gz"
        )
        cls.source_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_source.vcf.gz"
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
                "gene_info": "MLH1",
                "gene_info_id": "8636",
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

    def test_compare_gene(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        expected_compared_genes = [
            {
                "gene_info_id": "8636",
                "gene_info": "MLH1",
                "pathogenic_class_status": "LOST_PATHOGENICITY",
                "pathogenic_class_old": "Pathogenic",
                "pathogenic_class_best_review_confidence_old": "no_assertion_criteria_provided",
                "pathogenic_class_new": "..",
                "pathogenic_class_best_review_confidence_new": "..",
                "review_confidence_status": "NO_REVIEW",
                "review_confidence_best_pathogenic_class_old": "Pathogenic",
                "review_confidence_old": "no_assertion_criteria_provided",
                "review_confidence_best_pathogenic_class_new": "..",
                "review_confidence_new": "..",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
            {
                "gene_info_id": "9636",
                "gene_info": "ISG15",
                "pathogenic_class_status": "NEW_PATHOGENICITY",
                "pathogenic_class_old": "..",
                "pathogenic_class_best_review_confidence_old": "..",
                "pathogenic_class_new": "Pathogenic",
                "pathogenic_class_best_review_confidence_new": "criteria_provided,_single_submitter",
                "review_confidence_status": "NEW_REVIEW",
                "review_confidence_best_pathogenic_class_old": "..",
                "review_confidence_old": "..",
                "review_confidence_best_pathogenic_class_new": "Pathogenic",
                "review_confidence_new": "criteria_provided,_single_submitter",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190915",
            },
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
        ]

        clinvar_vcf_comparator.compare_references()

        compared_genes = clinvar_vcf_comparator.compare_genes()
        self.assertCountEqual(expected_compared_genes, compared_genes)

    def test_write_variant_comparison(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")

        clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()
        clinvar_vcf_comparator.write_variant_comparison(
            compared_variants, output_directory
        )
        comparison_output = os.path.join(
            output_directory,
            f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    TEST_DIR,
                    f"parsers/vcf/data/expected_clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{clinvar_vcf_comparator.target_vcf.release}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    def test_write_gene_comparison(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")

        clinvar_vcf_comparator.compare_references()
        compared_genes = clinvar_vcf_comparator.compare_genes()
        clinvar_vcf_comparator.write_gene_comparison(compared_genes, output_directory)
        comparison_output = os.path.join(
            output_directory,
            f"clinvar_gene_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    TEST_DIR,
                    f"parsers/vcf/data/expected_clinvar_gene_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{clinvar_vcf_comparator.target_vcf.release}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    def test_write_clinvarome(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        output_directory = os.path.join(TEST_DIR, "parsers/vcf/data")
        clinvar_vcf_comparator.compare_references()
        clinvar_vcf_comparator.write_clinvarome(output_directory)
        comparison_output = os.path.join(
            output_directory,
            f"clinvarome_{clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    TEST_DIR,
                    f"parsers/vcf/data/expected_clinvarome_{clinvar_vcf_comparator.target_vcf.release}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

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
            f"clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{clinvar_vcf_comparator.target_vcf.release}.vcf",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    TEST_DIR,
                    f"parsers/vcf/data/expected_clinvar_variant_diff_from_{clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{clinvar_vcf_comparator.target_vcf.release}.vcf",
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


class MissingVariantClinvarVCFComparatorTestCase(ClinvarVCFComparatorTestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_target_with_lost_variants.vcf.gz"
        )
        cls.source_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar_source.vcf.gz"
        )
        cls.hg38_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/vcf_new_38.vcf.gz"
        )

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
                "name_clinvar_new": "20190925",
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
                "name_clinvar_new": "20190925",
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
                "name_clinvar_new": "20190925",
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
                "name_clinvar_new": "20190925",
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
                "name_clinvar_new": "20190925",
            },
            {
                "variant_id": "GRCh37_X_153296172_G_A",
                "clinvar_id": "402986",
                "old_classification": "Benign",
                "new_classification": "..",
                "breaking_change": "warning",
                "confidence": "criteria_provided,_single_submitter",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190925",
            },
            {
                "variant_id": "GRCh37_X_153296290_G_A",
                "clinvar_id": "571208",
                "old_classification": "Uncertain_significance",
                "new_classification": "..",
                "breaking_change": "warning",
                "confidence": "criteria_provided,_single_submitter",
                "gene_info": "ISG15",
                "gene_info_id": "9636",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190925",
            },
        ]

        clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = clinvar_vcf_comparator.compare_variants()
        self.assertTrue(variants_have_been_lost)
        self.assertCountEqual(expected_compared_variants, compared_variants)

    def test_compare_gene(self):
        clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=self.source_clinvar_vcf, target_vcf=self.target_clinvar_vcf
        )
        expected_compared_genes = [
            {
                "gene_info_id": "8636",
                "gene_info": "MLH1",
                "pathogenic_class_status": "LOST_PATHOGENICITY",
                "pathogenic_class_old": "Pathogenic",
                "pathogenic_class_best_review_confidence_old": "no_assertion_criteria_provided",
                "pathogenic_class_new": "..",
                "pathogenic_class_best_review_confidence_new": "..",
                "review_confidence_status": "NO_REVIEW",
                "review_confidence_best_pathogenic_class_old": "Pathogenic",
                "review_confidence_old": "no_assertion_criteria_provided",
                "review_confidence_best_pathogenic_class_new": "..",
                "review_confidence_new": "..",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190925",
            },
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
                "name_clinvar_new": "20190925",
            },
            {
                "gene_info_id": "9636",
                "gene_info": "ISG15",
                "pathogenic_class_status": "NEW_PATHOGENICITY",
                "pathogenic_class_old": "..",
                "pathogenic_class_best_review_confidence_old": "..",
                "pathogenic_class_new": "Pathogenic",
                "pathogenic_class_best_review_confidence_new": "criteria_provided,_single_submitter",
                "review_confidence_status": "NEW_REVIEW",
                "review_confidence_best_pathogenic_class_old": "..",
                "review_confidence_old": "..",
                "review_confidence_best_pathogenic_class_new": "Pathogenic",
                "review_confidence_new": "criteria_provided,_single_submitter",
                "name_clinvar_old": "20190805",
                "name_clinvar_new": "20190925",
            },
        ]

        clinvar_vcf_comparator.compare_references()

        compared_genes = clinvar_vcf_comparator.compare_genes()
        self.assertCountEqual(expected_compared_genes, compared_genes)
