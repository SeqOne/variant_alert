# coding=utf-8
# =========================================================================
# Modules
# =========================================================================
import os

import click

from logger.logger import get_logger
from parsers.vcf.clinvar import ClinvarVCFComparator
from variant_alert import VariantAlertParams

logger = get_logger(__name__)
# Object to pass group command parameters to child commands
pass_variant_alert_params = click.make_pass_decorator(VariantAlertParams, ensure=True)


@click.group(chain=True)  # So people can launch both subcommands at the same time
@click.option(
    "--output_directory",
    "-o",
    "output_directory",
    help="Output directory of variant comparison file. Default current path + "
    "clinvar_gene_diff_from_$release-source_to _$release-target.txt",
    type=click.Path(
        exists=False, file_okay=True, dir_okay=True, writable=True, readable=True
    ),
    default=os.getcwd()
)
@click.argument("vcf_source_path")
@click.argument("vcf_target_path")
@pass_variant_alert_params
def variant_alert(variant_alert_params, output_directory, vcf_source_path, vcf_target_path):
    """
    CLI for variant_alert package.

    VCF_SOURCE_PATH: Used as a reference
    VCF_TARGET_PATH: Compared to the source vcf
    """
    logger.info(" started variant-alert command")
    clinvar_vcf_comparator = ClinvarVCFComparator(vcf_source_path, vcf_target_path)
    genome_reference = clinvar_vcf_comparator.compare_references()
    logger.info(f"... genome reference of compared file is {genome_reference}")

    variant_alert_params.clinvar_vcf_comparator = clinvar_vcf_comparator
    variant_alert_params.output_directory = output_directory


@variant_alert.command()
@click.option(
    "--output_format",
    "-f",
    "output_format",
    help="Format of output file, either tsv or vcf. Default = tsv",
    type=click.Choice(["tsv", "vcf"]),
    default="tsv"
)
@pass_variant_alert_params
def compare_variant(variant_alert_params, output_format):
    """
    Output variants with change of pathogenicity\n
    """
    logger.info("... launched variant comparison")
    output_directory = variant_alert_params.output_directory
    clinvar_vcf_comparator = variant_alert_params.clinvar_vcf_comparator
    logger.info("...... comparing variants")
    compared_variants, has_lost_variants = clinvar_vcf_comparator.compare_variants()
    if has_lost_variants:
        logger.warning("Your new version of clinvar is missing variants from your source file")
    logger.info(f"...... writing output to directory {output_directory}")
    clinvar_vcf_comparator.write_variant_comparison(compared_variants, output_directory, output_format)
    logger.info("... done")


@variant_alert.command()
@pass_variant_alert_params
def compare_gene(variant_alert_params):
    """
    Output genes with change of pathogenicity and/or review status\n
    """
    logger.info("... launched gene comparison")
    output_directory = variant_alert_params.output_directory
    clinvar_vcf_comparator = variant_alert_params.clinvar_vcf_comparator
    logger.info("...... comparing gene")
    compared_genes = clinvar_vcf_comparator.compare_genes()
    logger.info(f"...... writing output to directory {output_directory}")
    clinvar_vcf_comparator.write_gene_comparison(compared_genes, output_directory)
    logger.info("... done")


@variant_alert.command()
@pass_variant_alert_params
def clinvarome(variant_alert_params):
    """
    Output all pathogenic genes of target clinvar vcf\n
    """
    logger.info("... launched clinvarome creation")
    output_directory = variant_alert_params.output_directory
    clinvar_vcf_comparator = variant_alert_params.clinvar_vcf_comparator
    logger.info(f"...... writing output to directory {output_directory}")
    clinvar_vcf_comparator.write_clinvarome(output_directory)
    logger.info("... done")
