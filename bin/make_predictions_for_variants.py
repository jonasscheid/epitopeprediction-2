#!/usr/bin/env python

import sys
import logging
import vcf

from epytope.Core.Variant import Variant, VariationType, MutationSyntax

__author__ = "Christopher Mohr"
VERSION = "1.1"

# instantiate global logger object
logger = logging.getLogger(__name__)
# turn off passing of messages to root logger
# logger.propagate = False
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)

def read_vcf(filename: str):
    """
    reads vcf files
    returns a list of epytope variants
    :param filename: /path/to/file
    :param boolean pass_only: only consider variants that passed the filter (default: True)
    :return: list of epytope variants
    """

    #parse variant of vcf file into variants list
    vcf_reader = vcf.Reader(open(filename, 'rt'))
    variants = [record for record in vcf_reader]
    for record in variants:
        #print(record.CHROM, record.POS, record.ID, record.REF, record.ALT, record.QUAL, record.FILTER, record.INFO)
        pass

    #determine filetype
    FILE_TYPE = None
    VEP_KEY = "CSQ"
    SNPEFF_KEY = "ANN"

    if VEP_KEY in vcf_reader.infos:
        FILE_TYPE = "VEP"
    elif SNPEFF_KEY in vcf_reader.infos:
        FILE_TYPE = "SNPEFF"
    else:
        logging.error("Cannot determine the filetype vep or snpeff.")
        sys.exit(2)

    if FILE_TYPE == "VEP":
        parse_vep_record_info(variants, VEP_KEY)
    elif FILE_TYPE == "SNPEFF":
        parse_snpeff_record_info(variants, SNPEFF_KEY)

    return variants


def parse_vep_record_info(variants, key):
    '''description = ''
    description_list = description.split(" | ")
    for record in variants:
        relevant_information = record.INFO[key][0]
        relevant_information_list = relevant_information.split("|")
        #put all information with the description into a dict that
        record.INFO = {k: v for k, v in zip(description_list, relevant_information_list)}'''
    pass

def parse_snpeff_record_info(variants, key):
    description = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
    description_list = description.split(" | ")
    for record in variants:
        relevant_information = record.INFO[key][0]
        relevant_information_list = relevant_information.split("|")
        #put all information with the description into a dict that
        record.INFO = {k: v for k, v in zip(description_list, relevant_information_list)}

def __main__():
    #read_vcf(filename = 'variants.vcf')
    #gerade implementiert: Lesen von VCF files und parsing von dem record.INFO für snpeff vcfs aber noch nicht für vep
    pass

if __name__ == "__main__":
    __main__()
