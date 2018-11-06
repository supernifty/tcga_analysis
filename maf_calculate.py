#!/usr/bin/env python
'''
  stats on tcga maf files
'''

import argparse
import collections
import csv
import gzip
import logging
import sys

CAPTURE_SIZE=36000000

# Hugo_Symbol     Entrez_Gene_Id  Center  NCBI_Build      Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status        Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode     Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1        Tumor_Validation_Allele2        Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status     Validation_Status       Mutation_Status Sequencing_Phase        Sequence_Source Validation_Method       Score   BAM_File        Sequencer       Tumor_Sample_UUID       Matched_Norm_Sample_UUID        HGVSc   HGVSp   HGVSp_Short     Transcript_ID   Exon_Number     t_depth t_ref_count     t_alt_count     n_depth n_ref_count     n_alt_count     all_effects     Allele  Gene    Feature Feature_type    One_Consequence Consequence     cDNA_position   CDS_position    Protein_position        Amino_acids Codons  Existing_variation      ALLELE_NUM      DISTANCE        TRANSCRIPT_STRAND       SYMBOL  SYMBOL_SOURCE   HGNC_ID BIOTYPE CANONICAL       CCDS    ENSP    SWISSPROT       TREMBL  UNIPARC RefSeq  SIFT    PolyPhen        EXON    INTRON  DOMAINS GMAF    AFR_MAF AMR_MAF ASN_MAF EAS_MAF EUR_MAF SAS_MAF AA_MAF  EA_MAF  CLIN_SIG        SOMATIC PUBMED  MOTIF_NAME      MOTIF_POS       HIGH_INF_POS    MOTIF_SCORE_CHANGE      IMPACT  PICK    VARIANT_CLASS   TSL     HGVS_OFFSET     PHENO   MINIMISED       ExAC_AF ExAC_AF_Adj     ExAC_AF_AFR     ExAC_AF_AMR     ExAC_AF_EAS     ExAC_AF_FIN     ExAC_AF_NFE     ExAC_AF_OTH     ExAC_AF_SAS     GENE_PHENO      FILTER  CONTEXT src_vcf_id      tumor_bam_uuid  normal_bam_uuid case_id GDC_FILTER      COSMIC  MC3_Overlap     GDC_Validation_Status
# C1orf94 84970   BCM     GRCh38  chr1    34218878        34218878        +       3'UTR   SNP     C       C       T       novel           TCGA-AG-3727-01A-01W-0899-10    TCGA-AG-3727-10A-01W-0901-10 Somatic                                         Illumina HiSeq 2000     7cbdbb30-2abc-4440-9db3-f219bc358aad    44fbb860-5ffb-458b-a1a6-a52b749d9f9c    c.*117C>T                       ENST00000488417 7/7     263     249     14      175 C1orf94,3_prime_UTR_variant,,ENST00000488417,NM_001134734.1,c.*117C>T,MODIFIER,YES,,,1;C1orf94,3_prime_UTR_variant,,ENST00000373374,NM_032884.4,c.*117C>T,MODIFIER,,,,1 T       ENSG00000142698 ENST00000488417 Transcript      3_prime_UTR_variant     3_prime_UTR_variant     2034/2287       -/1797  -/598                           1               1       C1orf94 HGNC    HGNC:28250      protein_coding  YES     CCDS44108.1     ENSP00000435634 Q6P1W5          UPI0000D4BFB0   NM_001134734.1 7/7                                                                                                                                                     MODIFIER        1       SNV     1                       1 PASS    GGTCACAGACA     a6945d7d-472c-4da7-a122-de6d014564e9    12f0b883-900e-41b4-bbdf-a66b96e825ca    8d36e3fc-bd2c-441d-a4e1-ed7c94d2cb3e    e6827400-0d95-46b0-8874-6ce9e9d5011b    wga_pair                True    Unknown

SILENT = set([
  "3'Flank",
  "3'UTR",
  "5'Flank",
  "5'UTR",
  "Intron",
  "Silent",
  "Splice_Region",
  "Splice_Site",
  "Translation_Start_Site"
])

NON_SILENT=set([
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "IGR", # ???
  "In_Frame_Del",
  "In_Frame_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "RNA" # ???
])

def get_value(header, col, row):
  return row[header.index(col)]

def main(mafs, prefix):
  logging.info('starting...')

  # stuff to compute:
  # - mutations per sample
  sample_mutations_silent = collections.defaultdict(int)
  sample_mutations_non_silent = collections.defaultdict(int)

  # - number of samples mutated per gene
  mutated_genes = collections.defaultdict(set)

  samples_seen = set()
  genes_seen = set()

  for maf in mafs:
    logging.info('reading %s...', maf)
    header = None
    idx = 0
    for idx, row in enumerate(csv.reader(gzip.open(maf, 'r'), delimiter='\t')):
      if row[0].startswith('#'):
        continue
      if header is None:
        header = row
        continue
      
      gene = get_value(header, 'Hugo_Symbol', row)
      sample = get_value(header, 'Tumor_Sample_Barcode', row)
      mutation_type = get_value(header, 'Variant_Classification', row)

      genes_seen.add(gene)
      samples_seen.add(sample)

      if mutation_type in SILENT:
        sample_mutations_silent[sample] += 1
      elif mutation_type in NON_SILENT:
        sample_mutations_non_silent[sample] += 1
      else:
        logging.warn('unrecognised variant type %s', mutation_type)
      mutated_genes[gene].add(sample)

      if idx % 100000 == 0:
        logging.debug('%i lines', idx)

    logging.info('reading %s: %i lines done', maf, idx + 1)

  logging.info('%i genes. %i samples', len(genes_seen), len(samples_seen))

  if prefix is not None:
    logging.info('reading mafs: done. writing...')
    with open('{}.sample_mutation_rate.tsv'.format(prefix), 'w') as fh:
      fh.write('Sample\tSilent\tNon-silent\n')
      for sample in samples_seen:
        #fh.write('{}\t{:.2f}\t{:.2f}\n'.format(sample, sample_mutations_silent[sample], sample_mutations_non_silent[sample]))
        fh.write('{}\t{:.2f}\t{:.2f}\n'.format(sample, 1000000.0 * sample_mutations_silent[sample] / CAPTURE_SIZE, 1000000.0 * sample_mutations_non_silent[sample] / CAPTURE_SIZE))
  
    with open('{}.gene_mutation_frequency.tsv'.format(prefix), 'w') as fh:
      for gene in genes_seen:
        fh.write('{}\t{:.2f}\n'.format(gene, 1.0 * len(mutated_genes[gene]) / len(samples_seen)))
        #fh.write('{}\t{:.2f}\n'.format(gene, len(mutated_genes[gene])))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generate stats on tcga maf files')
  parser.add_argument('--mafs', required=True, nargs='+', help='maf files to analyse')
  parser.add_argument('--prefix', help='filename prefix')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.mafs, args.prefix)
