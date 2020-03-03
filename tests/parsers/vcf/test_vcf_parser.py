import os
import unittest
from datetime import date

from parsers.vcf import VCFParser
from tests import TEST_DIR


class VCFParserTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.vcf_with_caller = os.path.join(
            TEST_DIR, "parsers/vcf/data/raw/toy_withcaller.vcf.gz"
        )
        self.vcf_without_caller = os.path.join(
            TEST_DIR, "parsers/vcf/data/raw/toy.vcf.gz"
        )

    def test_vcf_parser(self):
        vcf_parser = VCFParser(self.vcf_without_caller)

        self.assertEqual(date.today(), vcf_parser.file_date)
        self.assertEqual("GRCh37", vcf_parser.reference)
        self.assertEqual("VCFv4.2", vcf_parser.vcf_file.header.version)
        self.assertEqual(
            [
                "e83f121a-ad81-4c2f-acc7-eaf050659cf5",
                "0647cfa8-b105-4755-905a-d60e9f033525",
            ],
            vcf_parser.samples,
        )

        first_variant = next(vcf_parser.vcf_file.fetch())
        self.assertEqual("1", first_variant.chrom)
        self.assertEqual(2340056, first_variant.pos)
        self.assertEqual("C", first_variant.ref)
        self.assertEqual("T", first_variant.alts[0])
        self.assertEqual(
            ("Benign/Likely_benign", "_other"), first_variant.info.get("CLNSIG")
        )

        sample_name, format_fields = first_variant.samples.items()[0]
        self.assertEqual("e83f121a-ad81-4c2f-acc7-eaf050659cf5", sample_name)
        format_field_name, format_field_value = format_fields.items()[6]
        self.assertEqual("PL2", format_field_name)
        self.assertEqual((0, 199, 200), format_field_value)

    def test_vcf_parser_with_caller(self):
        vcf_parser = VCFParser(self.vcf_with_caller)

        header_name, header_value = vcf_parser.vcf_file.header.formats.items()[4]
        self.assertEqual("CALLER", header_name)
        self.assertEqual(
            ["CALLER", "String", "Name of caller"],
            [header_value.name, header_value.type, header_value.description],
        )

        first_variant = next(vcf_parser.vcf_file.fetch())
        sample_name, format_fields = first_variant.samples.items()[0]
        format_field_name, format_field_value = format_fields.items()[-1]
        self.assertEqual("CALLER", format_field_name)
        self.assertEqual((".",), format_field_value)

    def test_set_variant_id(self):
        vcf_parser = VCFParser(self.vcf_with_caller)
        expected_variant_id = "1_123456_A_G"

        variant_id = vcf_parser.set_variant_id("1", 123456, "A", "G")

        self.assertEqual(expected_variant_id, variant_id)

    def test_set_genome_variant_id(self):
        vcf_parser = VCFParser(self.vcf_with_caller)
        expected_genome_variant_id = "GRCh37_1_123456_A_G"

        genome_variant_id = vcf_parser.set_genome_variant_id(
            vcf_parser.reference, "1", 123456, "A", "G"
        )

        self.assertEqual(expected_genome_variant_id, genome_variant_id)

    def test_split_piped_annotation_header(self):
        vcf_parser = VCFParser(self.vcf_with_caller)
        expected_info_subfields_by_rank = {
            "CSQ": {
                1: "Allele",
                2: "Consequence",
                3: "IMPACT",
                4: "SYMBOL",
                5: "Gene",
                6: "Feature_type",
                7: "Feature",
                8: "BIOTYPE",
                9: "EXON",
                10: "INTRON",
                11: "HGVSc",
                12: "HGVSp",
                13: "cDNA_position",
                14: "CDS_position",
                15: "Protein_position",
                16: "Amino_acids",
                17: "Codons",
                18: "Existing_variation",
                19: "DISTANCE",
                20: "STRAND",
                21: "FLAGS",
                22: "VARIANT_CLASS",
                23: "SYMBOL_SOURCE",
                24: "HGNC_ID",
                25: "CCDS",
                26: "SWISSPROT",
                27: "TREMBL",
                28: "UNIPARC",
                29: "REFSEQ_MATCH",
                30: "GIVEN_REF",
                31: "USED_REF",
                32: "BAM_EDIT",
                33: "GENE_PHENO",
                34: "SIFT",
                35: "PolyPhen",
                36: "DOMAINS",
                37: "HGVS_OFFSET",
                38: "MaxEntScan_alt",
                39: "MaxEntScan_diff",
                40: "MaxEntScan_ref",
                41: "GeneSplicer",
                42: "AltaiNeandertal",
                43: "Ancestral_allele",
                44: "CADD_phred",
                45: "CADD_raw",
                46: "CADD_raw_rankscore",
                47: "DANN_rankscore",
                48: "DANN_score",
                49: "Denisova",
                50: "Eigen-PC-phred",
                51: "Eigen-PC-raw",
                52: "Eigen-PC-raw_rankscore",
                53: "Eigen-phred",
                54: "Eigen-raw",
                55: "Eigen_coding_or_noncoding",
                56: "Ensembl_geneid",
                57: "Ensembl_proteinid",
                58: "Ensembl_transcriptid",
                59: "FATHMM_converted_rankscore",
                60: "FATHMM_pred",
                61: "FATHMM_score",
                62: "GERP++_NR",
                63: "GERP++_RS",
                64: "GERP++_RS_rankscore",
                65: "GM12878_confidence_value",
                66: "GM12878_fitCons_score",
                67: "GM12878_fitCons_score_rankscore",
                68: "GTEx_V6p_gene",
                69: "GTEx_V6p_tissue",
                70: "GenoCanyon_score",
                71: "GenoCanyon_score_rankscore",
                72: "H1-hESC_confidence_value",
                73: "H1-hESC_fitCons_score",
                74: "H1-hESC_fitCons_score_rankscore",
                75: "HUVEC_confidence_value",
                76: "HUVEC_fitCons_score",
                77: "HUVEC_fitCons_score_rankscore",
                78: "Interpro_domain",
                79: "LRT_Omega",
                80: "LRT_converted_rankscore",
                81: "LRT_pred",
                82: "LRT_score",
                83: "M-CAP_pred",
                84: "M-CAP_rankscore",
                85: "M-CAP_score",
                86: "MetaLR_pred",
                87: "MetaLR_rankscore",
                88: "MetaLR_score",
                89: "MetaSVM_pred",
                90: "MetaSVM_rankscore",
                91: "MetaSVM_score",
                92: "MutPred_AAchange",
                93: "MutPred_Top5features",
                94: "MutPred_protID",
                95: "MutPred_rankscore",
                96: "MutPred_score",
                97: "MutationAssessor_UniprotID",
                98: "MutationAssessor_pred",
                99: "MutationAssessor_score",
                100: "MutationAssessor_score_rankscore",
                101: "MutationAssessor_variant",
                102: "MutationTaster_AAE",
                103: "MutationTaster_converted_rankscore",
                104: "MutationTaster_model",
                105: "MutationTaster_pred",
                106: "MutationTaster_score",
                107: "PROVEAN_converted_rankscore",
                108: "PROVEAN_pred",
                109: "PROVEAN_score",
                110: "Polyphen2_HDIV_pred",
                111: "Polyphen2_HDIV_rankscore",
                112: "Polyphen2_HDIV_score",
                113: "Polyphen2_HVAR_pred",
                114: "Polyphen2_HVAR_rankscore",
                115: "Polyphen2_HVAR_score",
                116: "REVEL_rankscore",
                117: "REVEL_score",
                118: "Reliability_index",
                119: "SIFT_converted_rankscore",
                120: "SIFT_pred",
                121: "SIFT_score",
                122: "SiPhy_29way_logOdds",
                123: "SiPhy_29way_logOdds_rankscore",
                124: "SiPhy_29way_pi",
                125: "Transcript_id_VEST3",
                126: "Transcript_var_VEST3",
                127: "Uniprot_aapos_Polyphen2",
                128: "Uniprot_acc_Polyphen2",
                129: "Uniprot_id_Polyphen2",
                130: "VEST3_rankscore",
                131: "VEST3_score",
                132: "aapos",
                133: "cds_strand",
                134: "clinvar_clnsig",
                135: "clinvar_golden_stars",
                136: "clinvar_rs",
                137: "codon_degeneracy",
                138: "codonpos",
                139: "fathmm-MKL_coding_group",
                140: "fathmm-MKL_coding_pred",
                141: "fathmm-MKL_coding_rankscore",
                142: "fathmm-MKL_coding_score",
                143: "integrated_confidence_value",
                144: "integrated_fitCons_score",
                145: "integrated_fitCons_score_rankscore",
                146: "phastCons100way_vertebrate",
                147: "phastCons100way_vertebrate_rankscore",
                148: "phastCons20way_mammalian",
                149: "phastCons20way_mammalian_rankscore",
                150: "phyloP100way_vertebrate",
                151: "phyloP100way_vertebrate_rankscore",
                152: "phyloP20way_mammalian",
                153: "phyloP20way_mammalian_rankscore",
                154: "refcodon",
                155: "ada_score",
                156: "rf_score",
            },
            "DTG": {1: "Gene", 2: "medgen_name", 3: "medgen_id", 4: "short_name"},
        }

        info_subfields_by_rank = {}
        for info_key, info_items in vcf_parser.vcf_file.header.info.items():
            if info_key in ["CSQ", "DTG"]:
                info_subfields_by_rank.update(
                    vcf_parser._split_piped_annotation_header(info_key, info_items)
                )
        self.assertCountEqual(expected_info_subfields_by_rank, info_subfields_by_rank)

    def test_split_piped_annotation(self):
        vcf_parser = VCFParser(self.vcf_with_caller)

        annot_field = "CSQ"
        # @formatters:off
        annot_field_value = (
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_000141.4|protein_coding|14/18||NM_000141.4:c.1937G>C|NP_000132.3:p.Gly646Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144913.1|protein_coding|13/17||NM_001144913.1:c.1940G>C|NP_001138385.1:p.Gly647Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144914.1|protein_coding|11/15||NM_001144914.1:c.1601G>C|NP_001138386.1:p.Gly534Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144915.1|protein_coding|13/17||NM_001144915.1:c.1670G>C|NP_001138387.1:p.Gly557Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144917.1|protein_coding|12/16||NM_001144917.1:c.1589G>C|NP_001138389.1:p.Gly530Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144918.1|protein_coding|12/16||NM_001144918.1:c.1586G>C|NP_001138390.1:p.Gly529Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001144919.1|protein_coding|13/17||NM_001144919.1:c.1673G>C|NP_001138391.1:p.Gly558Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001320654.1|protein_coding|13/17||NM_001320654.1:c.1253G>C|NP_001307583.1:p.Gly418Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_001320658.1|protein_coding|14/18||NM_001320658.1:c.1931G>C|NP_001307587.1:p.Gly644Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_022970.3|protein_coding|14/18||NM_022970.3:c.1940G>C|NP_075259.4:p.Gly647Ala",  # noqa E501
            "G|missense_variant|MODERATE|FGFR2|2263|Transcript|NM_023029.2|protein_coding|12/16||NM_023029.2:c.1670G>C|NP_075418.1:p.Gly557Ala|1820|167",  # noqa E501
        )
        expected_annot_subvalue_lists = [
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_000141.4",
                "protein_coding",
                "14/18",
                "",
                "NM_000141.4:c.1937G>C",
                "NP_000132.3:p.Gly646Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144913.1",
                "protein_coding",
                "13/17",
                "",
                "NM_001144913.1:c.1940G>C",
                "NP_001138385.1:p.Gly647Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144914.1",
                "protein_coding",
                "11/15",
                "",
                "NM_001144914.1:c.1601G>C",
                "NP_001138386.1:p.Gly534Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144915.1",
                "protein_coding",
                "13/17",
                "",
                "NM_001144915.1:c.1670G>C",
                "NP_001138387.1:p.Gly557Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144917.1",
                "protein_coding",
                "12/16",
                "",
                "NM_001144917.1:c.1589G>C",
                "NP_001138389.1:p.Gly530Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144918.1",
                "protein_coding",
                "12/16",
                "",
                "NM_001144918.1:c.1586G>C",
                "NP_001138390.1:p.Gly529Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001144919.1",
                "protein_coding",
                "13/17",
                "",
                "NM_001144919.1:c.1673G>C",
                "NP_001138391.1:p.Gly558Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001320654.1",
                "protein_coding",
                "13/17",
                "",
                "NM_001320654.1:c.1253G>C",
                "NP_001307583.1:p.Gly418Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_001320658.1",
                "protein_coding",
                "14/18",
                "",
                "NM_001320658.1:c.1931G>C",
                "NP_001307587.1:p.Gly644Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_022970.3",
                "protein_coding",
                "14/18",
                "",
                "NM_022970.3:c.1940G>C",
                "NP_075259.4:p.Gly647Ala",
            ],
            [
                "G",
                "missense_variant",
                "MODERATE",
                "FGFR2",
                "2263",
                "Transcript",
                "NM_023029.2",
                "protein_coding",
                "12/16",
                "",
                "NM_023029.2:c.1670G>C",
                "NP_075418.1:p.Gly557Ala",
                "1820",
                "167",
            ],
        ]
        # @formatters:on

        annot_subvalue_lists = vcf_parser._split_piped_annotation(
            annot_field, annot_field_value
        )
        self.assertCountEqual(expected_annot_subvalue_lists, annot_subvalue_lists)

    def test_get_variant_buffer(self):
        vcf_parser = VCFParser(self.vcf_with_caller)
        expected_variant_buffer = {
            "GRCh37_10_123247554_C_G": {
                "chrom": "10",
                "pos": 123247554,
                "id": None,
                "ref": "C",
                "alt": "G",
            },
            "GRCh37_1_2340056_C_T": {
                "chrom": "1",
                "pos": 2340056,
                "id": None,
                "ref": "C",
                "alt": "T",
            },
            "GRCh37_1_161599571_T_.": {
                "chrom": "1",
                "pos": 161599571,
                "id": "36924",
                "ref": "T",
                "alt": ".",
            },
            "GRCh37_2_238247758_CACCTAAAGAAAAAAA_.": {
                "chrom": "2",
                "pos": 238247758,
                "id": "504370",
                "ref": "CACCTAAAGAAAAAAA",
                "alt": ".",
            },
        }

        variant_buffer = vcf_parser.get_variant_buffer(
            ["CHROM", "POS", "ID", "REF", "ALT"]
        )

        self.assertCountEqual(expected_variant_buffer, variant_buffer)

        variant_buffer_with_no_input_fields = vcf_parser.get_variant_buffer()
        self.assertCountEqual(
            expected_variant_buffer, variant_buffer_with_no_input_fields
        )
