# coding=utf-8
# =========================================================================
# Modules
# =========================================================================
import datetime
import os
import re
from datetime import date
from typing import List

import pysam

from logger.logger import get_logger


class MalformedVCFError(Exception):
    pass


class NotSameReferenceError(Exception):
    pass


class VCFParser(object):
    @property
    def get_field_values(self):
        return self._get_fields_dict_methods

    def __init__(self, vcf_file_path):
        self.logger = get_logger(__name__)

        # vcf file infos
        self.vcf_file_path = vcf_file_path
        self.vcf_file_name = os.path.basename(vcf_file_path)
        self.vcf_type_name = "VCF"
        self.reference = None
        self.file_date = None
        self.samples = []
        self._get_fields_dict_methods = {
            "CHROM": self._get_variant_record_chrom,
            "POS": self._get_variant_record_pos,
            "ID": self._get_variant_record_id,
            "REF": self._get_variant_record_ref,
            "ALT": self._get_variant_record_alt,
        }

        self.logger.info(f" VCF Parser initialized with input {self.__str__()}")
        # Reading VCF file with pysam
        self._get_records_from_file()
        self._set_date()
        self._set_reference()
        self._set_samples()

    @staticmethod
    def _get_variant_record_chrom(variant_record):
        return {"chrom": variant_record.chrom}

    @staticmethod
    def _get_variant_record_pos(variant_record):
        return {"pos": variant_record.pos}

    @staticmethod
    def _get_variant_record_id(variant_record):
        return {"id": variant_record.id}

    @staticmethod
    def _get_variant_record_ref(variant_record):
        return {"ref": variant_record.ref}

    @staticmethod
    def _get_variant_record_alt(variant_record):
        if variant_record.alts is None:
            alt = "."
        else:
            alt = variant_record.alts[0]
        return {"alt": alt}

    @staticmethod
    def set_genome_variant_id(
        genome_version: str, chrom: str, pos: int, ref: str, alt: str
    ) -> str:
        return "_".join([genome_version, chrom, str(pos), ref, alt])

    @staticmethod
    def set_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
        return "_".join([chrom, str(pos), ref, alt])

    @staticmethod
    def _split_piped_annotation_header(info_key: str, info_items):
        """
        Function to split Header used to annotate variants so when parsing the variant records we always know
        to which annotation a column relates. The format of these annotations is ANNOT=a|b|c|d,a|b|e|d,f|b|e|g
        Values and rank are mapped in a dict attribute so we can store multiple fields of this type at the same time.
        """
        if info_items.description is not None:
            try:
                part_thrown_away, piped_annotations = info_items.description.split(
                    ":", 1
                )
            except ValueError:
                raise MalformedVCFError(
                    f"Description of Annotator info field {info_key} is not has expected, "
                    f"can't split sub-annotations: {info_items.record}"
                )
            else:
                # Remove space and single quotes around SnpEff annots
                annotations = piped_annotations.replace("'", "")
                annotation_list = [
                    annotation.strip() for annotation in annotations.split("|")
                ]
                info_subfields_by_rank = {}
                for rank, annotation in enumerate(annotation_list, 1):
                    info_subfields_by_rank.setdefault(info_key, {}).setdefault(
                        rank, annotation
                    )
                return info_subfields_by_rank
        else:
            raise MalformedVCFError(
                f"Annotator info field {info_key} has no Description, "
                f"can't split annotation: {info_items.header}"
            )

    def _split_piped_annotation(self, annot_field, annot_field_value):
        # Fields are of type FIELD_NAME=(a|b|c|d|e,a|h|k|d|u)
        try:
            annot_subvalue_lists = []
            for annot_subvalue_list in annot_field_value:  # Pysam set it as a tuple
                if annot_subvalue_list is None:
                    self.logger.warning(
                        "Found allele variant annotations with a comma seperating nothing "
                        f"for info field {annot_field} of variant {annot_field_value}"
                    )
                elif annot_subvalue_list == ".":
                    continue
                else:
                    annot_subvalue_lists.append(annot_subvalue_list.split("|"))
            return annot_subvalue_lists
        except AttributeError as e:
            raise MalformedVCFError(e)

    def __str__(self):
        return f"file: {self.vcf_file_path}"

    def _get_records_from_file(self):
        """
        method to instantiate pysam VariantFile class from input VCF file,
        check if the VCF format is supported
        """
        self.logger.info(f"... VCFParser: Reading input vcf file {self.vcf_file_path}")
        self.vcf_file = pysam.VariantFile(self.vcf_file_path, "r")

        if self.vcf_file.header.version in ["VCFv4.1", "VCFv4.2"]:
            pass
        else:
            raise IOError(
                f"VCF format {self.vcf_file.header.version} is not supported. Supported ones are VCFv4.2 and VCFv4.1"
            )

    def _set_date(self):
        """
        Look for a file date in the VCF header, if not available,
        send a meany warning to the user
        """
        self.logger.info("... Setting date")

        for header in self.vcf_file.header.records:
            if header.key == "fileDate":
                self.logger.info("... Found a defined date")
                if self.file_date:
                    raise MalformedVCFError("VCF file date has already been set")
                else:
                    try:
                        self.file_date = datetime.datetime.strptime(
                            header.value, "%Y%m%d"
                        ).strftime("%Y-%m-%d")
                    except ValueError:
                        try:
                            self.file_date = datetime.datetime.strptime(
                                header.value, "%Y-%m-%d"
                            )
                        except ValueError:
                            self.logger.error(
                                f"Failed to get VCF file date from {header.value}. Setting it to today's date"
                            )
                            self.file_date = date.today()
                    else:
                        self.logger.info(f"...... VCF file date is: {self.file_date}")
        else:
            if self.file_date is None:
                self.logger.warning(
                    "... No file date set in your VCF header to help us tracking your analysis,"
                    " Shame on you! ... Shame ! ... Shame ! ..."
                )
                self.file_date = date.today()

    def _set_reference(self):
        """
        Look for a genome reference file path in the VCF header, if not available,
        send a meany warning to the user
        """
        self.logger.info("... Setting reference")

        for header in self.vcf_file.header.records:
            if header.key == "reference":
                self.logger.info("... Found a defined reference")
                if self.reference:
                    raise MalformedVCFError("VCF reference file has already been set")
                else:
                    self.reference = re.search(
                        r"(?P<genome>(GRCh|hg)\d+)", header.value
                    ).group("genome")
                    self.logger.info(f"...... VCF reference file is: {self.reference}")

    def _set_samples(self):
        self.logger.info("... Setting VCF samples")
        for sample in self.vcf_file.header.samples:
            self.samples.append(sample)

    # WARNING : does not currently handle multi-allelic alts
    def get_variant_buffer(self, fields: List[str] = None):
        variants = {}
        for variant_record in self.vcf_file.fetch():
            if variant_record.alts is None:
                alt = "."
            elif len(variant_record.alts) > 1:
                raise ValueError(
                    f"Variant buffer generation does not currently support {self.vcf_type_name} file"
                    " with multi-allelic alts on a same line."
                )
            else:
                alt = variant_record.alts[0]

            variant_annotations = {}
            # by default, use all available fields
            if fields is None:
                fields = self.get_field_values.keys()
            for field in fields:
                try:
                    variant_annotations.update(
                        self.get_field_values[field](variant_record)
                    )
                # If one of the given fields does not exists log an explicit message and re-raise
                except KeyError:
                    self.logger.exception(
                        "given field {} is not available. Either choose one of the available ones "
                        "or implement your field. Available fields: \n {}".format(
                            field, "\n".join(self.get_field_values.keys())
                        )
                    )
                    raise
            variants[
                self.set_genome_variant_id(
                    self.reference,
                    variant_record.chrom,
                    variant_record.pos,
                    variant_record.ref,
                    alt,
                )
            ] = variant_annotations

        return variants
