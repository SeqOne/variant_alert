import json
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
        self.vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/clinvar_20191231.vcf.gz"
        )
        self.parser = ClinvarVCFParser(self.vcf)

    def test_init(self):
        self.assertEqual("Clinvar VCF", self.parser.vcf_type_name)

        self.assertEqual("2019-12-31", self.parser.file_date.strftime("%Y-%m-%d"))
        self.assertEqual("20191231", self.parser.release)
        self.assertIsNotNone(self.parser.variant_buffer)
        self.assertIsNotNone(self.parser.gene_buffer)

    def test_clinvar_vcf_parser_variant_buffer(self):
        with open(
            os.path.join(
                TEST_DIR, "parsers/vcf/data/clinvar/expected_variant_buffer.json"
            )
        ) as gold_standard:
            self.assertCountEqual(
                json.load(gold_standard), self.parser.variant_buffer,
            )

    def test_clinvar_vcf_parser_gene_buffer(self):
        with open(
            os.path.join(TEST_DIR, "parsers/vcf/data/clinvar/expected_gene_buffer.json")
        ) as gold_standard:
            self.assertCountEqual(
                json.load(gold_standard), self.parser.gene_buffer,
            )


class ClinvarVCFComparatorTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/clinvar_target.vcf.gz"
        )
        cls.source_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/clinvar_source.vcf.gz"
        )
        cls.hg38_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/vcf_new_38.vcf.gz"
        )
        cls.output_directory = os.path.join(TEST_DIR, "parsers/vcf/data/clinvar")
        cls.output_ext_prefix = ""

        cls.source_clinvar_vcf_parser = ClinvarVCFParser(cls.source_clinvar_vcf)
        cls.clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=cls.source_clinvar_vcf, target_vcf=cls.target_clinvar_vcf
        )
        cls.clinvar_vcf_comparator_with_wrong_genome = ClinvarVCFComparator(
            cls.source_clinvar_vcf, cls.hg38_clinvar_vcf
        )

    def test_compare_reference(self):
        self.assertEqual(
            self.source_clinvar_vcf_parser.reference,
            self.clinvar_vcf_comparator.compare_references(),
        )

        with self.assertRaises(ClinVarComparisonError):
            self.clinvar_vcf_comparator_with_wrong_genome.compare_references()

    def test_compare_variants(self):
        with open(
            os.path.join(
                TEST_DIR,
                f"parsers/vcf/data/clinvar/expected_compared_variants_{self.clinvar_vcf_comparator.target_vcf.release}.json",
            )
        ) as gold_standard:
            self.clinvar_vcf_comparator.compare_references()
            (
                compared_variants,
                variants_have_been_lost,
            ) = self.clinvar_vcf_comparator.compare_variants()
            self.assertFalse(variants_have_been_lost)
            self.assertCountEqual(json.load(gold_standard), compared_variants)

    def test_compare_gene(self):
        with open(
            os.path.join(
                TEST_DIR,
                f"parsers/vcf/data/clinvar/expected_compared_genes_{self.clinvar_vcf_comparator.target_vcf.release}.json",
            )
        ) as gold_standard:
            self.clinvar_vcf_comparator.compare_references()

            compared_genes = self.clinvar_vcf_comparator.compare_genes()
            self.assertCountEqual(json.load(gold_standard), compared_genes)

    def test_write_variant_comparison(self):
        self.clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = self.clinvar_vcf_comparator.compare_variants()
        self.clinvar_vcf_comparator.write_variant_comparison(
            compared_variants, self.output_directory
        )
        comparison_output = os.path.join(
            self.output_directory,
            f"clinvar_variant_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{self.clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    self.output_directory,
                    f"expected_clinvar_variant_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{self.clinvar_vcf_comparator.target_vcf.release}{self.output_ext_prefix}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    def test_write_gene_comparison(self):
        self.clinvar_vcf_comparator.compare_references()
        compared_genes = self.clinvar_vcf_comparator.compare_genes()
        self.clinvar_vcf_comparator.write_gene_comparison(
            compared_genes, self.output_directory
        )
        comparison_output = os.path.join(
            self.output_directory,
            f"clinvar_gene_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{self.clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    self.output_directory,
                    f"expected_clinvar_gene_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{self.clinvar_vcf_comparator.target_vcf.release}{self.output_ext_prefix}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    def test_write_clinvarome(self):
        self.clinvar_vcf_comparator.compare_references()
        self.clinvar_vcf_comparator.write_clinvarome(self.output_directory)
        comparison_output = os.path.join(
            self.output_directory,
            f"clinvarome_{self.clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    self.output_directory,
                    f"expected_clinvarome_{self.clinvar_vcf_comparator.target_vcf.release}{self.output_ext_prefix}.tsv",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    def test_write_variant_comparison_as_vcf(self):
        self.clinvar_vcf_comparator.compare_references()
        (
            compared_variants,
            variants_have_been_lost,
        ) = self.clinvar_vcf_comparator.compare_variants()
        self.clinvar_vcf_comparator.write_variant_comparison_as_vcf(
            compared_variants, self.output_directory
        )
        comparison_output = os.path.join(
            self.output_directory,
            f"clinvar_variant_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
            f"_to_{self.clinvar_vcf_comparator.target_vcf.release}.vcf",
        )
        self.assertTrue(os.path.exists(comparison_output))
        with open(comparison_output) as co:
            with open(
                os.path.join(
                    self.output_directory,
                    f"expected_clinvar_variant_diff_from_{self.clinvar_vcf_comparator.source_vcf.release}"
                    f"_to_{self.clinvar_vcf_comparator.target_vcf.release}{self.output_ext_prefix}.vcf",
                )
            ) as gold_standard:
                self.assertCountEqual(gold_standard.readlines(), co.readlines())

    @classmethod
    def tearDownClass(cls) -> None:
        tsv_variant_comparison_output = os.path.join(
            cls.output_directory,
            (
                f"clinvar_variant_diff_from_{cls.clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{cls.clinvar_vcf_comparator.target_vcf.release}.tsv"
            ),
        )
        if os.path.exists(tsv_variant_comparison_output):
            os.remove(tsv_variant_comparison_output)

        vcf_variant_comparison_output = os.path.join(
            cls.output_directory,
            (
                f"clinvar_variant_diff_from_{cls.clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{cls.clinvar_vcf_comparator.target_vcf.release}.vcf"
            ),
        )
        if os.path.exists(vcf_variant_comparison_output):
            os.remove(vcf_variant_comparison_output)

        tsv_gene_comparison_output = os.path.join(
            cls.output_directory,
            (
                f"clinvar_gene_diff_from_{cls.clinvar_vcf_comparator.source_vcf.release}_"
                f"to_{cls.clinvar_vcf_comparator.target_vcf.release}.tsv"
            ),
        )
        if os.path.exists(tsv_gene_comparison_output):
            os.remove(tsv_gene_comparison_output)

        tsv_clinvarome_output = os.path.join(
            cls.output_directory,
            f"clinvarome_{cls.clinvar_vcf_comparator.target_vcf.release}.tsv",
        )
        if os.path.exists(tsv_clinvarome_output):
            os.remove(tsv_clinvarome_output)


class MissingVariantClinvarVCFComparatorTestCase(ClinvarVCFComparatorTestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.target_clinvar_vcf = os.path.join(
            TEST_DIR,
            "parsers/vcf/data/clinvar/clinvar_target_with_lost_variants.vcf.gz",
        )
        cls.source_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/clinvar_source.vcf.gz"
        )
        cls.hg38_clinvar_vcf = os.path.join(
            TEST_DIR, "parsers/vcf/data/clinvar/vcf_new_38.vcf.gz"
        )
        cls.output_directory = os.path.join(TEST_DIR, "parsers/vcf/data/clinvar")
        cls.output_ext_prefix = ""

        cls.source_clinvar_vcf_parser = ClinvarVCFParser(cls.source_clinvar_vcf)
        cls.clinvar_vcf_comparator = ClinvarVCFComparator(
            source_vcf=cls.source_clinvar_vcf, target_vcf=cls.target_clinvar_vcf
        )
        cls.clinvar_vcf_comparator_with_wrong_genome = ClinvarVCFComparator(
            cls.source_clinvar_vcf, cls.hg38_clinvar_vcf
        )

    def test_compare_variants(self):
        with open(
            os.path.join(
                TEST_DIR,
                f"parsers/vcf/data/clinvar/expected_compared_variants_{self.clinvar_vcf_comparator.target_vcf.release}.json",
            )
        ) as gold_standard:
            self.clinvar_vcf_comparator.compare_references()
            (
                compared_variants,
                variants_have_been_lost,
            ) = self.clinvar_vcf_comparator.compare_variants()
            self.assertTrue(variants_have_been_lost)
            self.assertCountEqual(json.load(gold_standard), compared_variants)
