'''
RiboSeq Pipeline
Author: Thomas Goodman

Pipeline to process RiboSeq fastqs for data analysis Alignment + QC

Pipeline created using Chunky pipes
'''

## TODO: Clean up code a little
## TODO: Make usable for paired end samples
## Get versions of software used 
## Fix picard

import matplotlib
matplotlib.use('Agg')
import os
import multiprocessing
import subprocess
import re
import time
import pysam
from chunkypipes.components import Software, Parameter, Redirect, Pipe, BasePipeline
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

defaultThreads = multiprocessing.cpu_count()


def split_bid(fastq):
  return fastq.split('/')[-1].split('.')[0]
def new_dir(outputDir, folder):
  return os.path.join(outputDir, folder)


class Pipeline(BasePipeline):
  def description(self):
    return ['Pipeline to process RiboSeq Data and do QC']

  # Required paths to run commands
  def configure(self):
    return{
      'cutadapt':{
        'path': 'Full path to cutadapt executable'
      },
      'STAR':{
        'path': 'Full path to STAR executable',
        'genomeDir': 'Full path to STAR indexed genome directory',
        'GTF_ref': 'Path to gtf annotation for genome'
      },
      'samtools':{
        'path': 'Full path to samtools'
      },
      'bowtie2':{
        'path': 'Full path to bowtie2',
        'genome_ref': 'Path to index of genome for bowtie2'
      },
      'bedtools':{
        'path': 'Full path to bedtools',
        'hg19_start100': 'Path to sorted start100 annotation bed file',
        'grch37_bed': 'Path to sorted grch37 bed file'
      },
      'read_distribution':{ # RSeq QC package
        'path': 'Full path to read_distribution.py', 
        'hg19_bed': 'Path to hg19_refseq bed12 file'
      },
      'FastQC':{
        'path' : 'Full path to FastQC software'
      },
      'featureCounts':{
        'path' : 'Full path to featureCounts software'
      },
      'picard':{
        'path'      : 'Full path and command to run picard ex: java -jar /mnt/picard/dist/picard.jar',
        'genomeFasta' : 'Full path to the genome fasta reference file'
      }
    }

  # Required input for pipeline software
  # Next time put all these paths to annotations in configure
  def add_pipeline_args(self, parser):
    parser.add_argument('--bam:lib', required=True, nargs='*', help='Fastq input for pipeline:library name(prefix for files)')
    parser.add_argument('--output', required=True, help='Where pipeline output should go')
    parser.add_argument('--adapter', default='AGATCGGAAGAGCACACGTCT', help='Adapter sequence for trimming')
    parser.add_argument('--threads', default=defaultThreads, help='Threads to be used for multi-threaded programs. Default is 8')
    

    # chunky run RiboSeq_pipe.py --fastqs 
    #  /mnt/cinder/thomas/RiboSeq/Lane5/AWS-3_S3_L005_R1_001.fastq.gz
    #  --output /mnt/cinder/thomas/RiboSeq/test --threads

  def run_pipeline(self, pipeline_args, pipeline_config):
    # create variables from parser if wanted
    bamFiles = pipeline_args['bam:lib']
    outputDir = pipeline_args['output']
    adapter = pipeline_args['adapter']
    numThreads = pipeline_args['threads']

    # Create output directory
    subprocess.call(['mkdir', outputDir])

    # Software
    cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
    star = Software('STAR', pipeline_config['STAR']['path'])
    bedtools = Software('bedtools', pipeline_config['bedtools']['path'])
    bowtie2 = Software('bowtie2', pipeline_config['bowtie2']['path'])
    samtools = Software('samtools', pipeline_config['samtools']['path'])
    samtools_sort = Software('samtools sort', pipeline_config['samtools']['path'])
    read_distribution = Software('read_distribution.py', 
      pipeline_config['read_distribution']['path'])
    featureCounts = Software('featureCounts', pipeline_config['featureCounts']['path'])
    fastQC = Software('FastQC', pipeline_config['FastQC']['path'])
    picard = Software('picard', pipeline_config['picard']['path'])

    # Change these to just be done in python script?

    # Common software tools 
    awk = Software('awk', 'awk')
    sort = Software('sort', 'sort')
    uniq = Software('uniq', 'uniq')
    paste = Software('paste', 'paste')
    cat = Software('cat', 'cat')
    grep = Software('grep', 'grep')

    # Directories and Files
    pathToGenomeDir = pipeline_config['STAR']['genomeDir']
    pathToGenome = pipeline_config['bowtie2']['genome_ref']
    pathToGtf = pipeline_config['STAR']['GTF_ref']
    pathTo_hg19_bed = pipeline_config['read_distribution']['hg19_bed']
    pathTo_hg19_bed_start100 = pipeline_config['bedtools']['hg19_start100']
    pathTo_grch37_bed = pipeline_config['bedtools']['grch37_bed']
    pathTo_genomeFasta = pipeline_config['picard']['genomeFasta']


    '''

      remove adaptor and trim
      adaptor sequence: AGATCGGAAGAGCACACGTCT
      -m 25 discard any reads shorter than 25 nucleotides
      keep only reads that had the adaptor sequence --discard-untrimmed

      cutadapt -a AGATCGGAAGAGCACACGTCT -m 25 --discard-untrimmed {filename}.fastq.gz
       > {filename}_trimmed.fastq.gz 2> {filename}_report.txt
      
      Remove adapters
      Only keep reads with adapters, otherwise artifact
      Discard reads shorter than 25 bp
      
    '''

    # Keep track of Bids in pipeline

    bid_list = []
    for bamLib in bamFiles:
      bid_list.append(bamLib.split(':')[-1])

    '''
      Sort and extract uniquely mapped reads for QC and further analyses
        samtools view -H $file > header.sam
        samtools view $file | grep -w NH:i:1 | cat header.sam - | samtools view -bS - | samtools sort - ${filename}_uniq_sorted
        rm header.sam

      Using this file for the rest of the analysis
    '''

    for bamLib in bamFiles:
      bam, bid = bamLib.split(':')
      newDir = new_dir(outputDir, bid)
      samtools.run(
        Parameter('view'),
        Parameter('-H'),
        Parameter(bam), # star outfile name
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.header.sam'.format(bid)))
      )
      samtools.run(
        Parameter('view'),
        Parameter(bam), # star outfile name
        Pipe(
          grep.pipe(
            Parameter('-w'),
            Parameter('NH:i:1'),
            Pipe(
              cat.pipe(
                Parameter(os.path.join(newDir, '{}.header.sam'.format(bid)), '-'),
                Pipe(
                  samtools.pipe(
                    Parameter('view'),
                    Parameter('-bS', '-'),
                    Pipe(
                      samtools.pipe(
                        Parameter('sort'),
                        Parameter('-', '-o', '{}/{}.uniq_sorted.bam'.format(newDir, bid))
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
      # subprocess.call(['rm', '{}/{}.header.sam'.format(newDir, bid)])

    '''
      SeQC to evaluate percent reads mapped to each genomic features
        read_distribution.py -r hg19_RefSeq.bed12 -i $file
    '''

    for bid in bid_list:
      newDir = new_dir(outputDir, bid)
      read_distribution.run(
        Parameter('-r'),
        Parameter(pathTo_hg19_bed),
        Parameter('-i'),
        Parameter('{}/{}.uniq_sorted.bam'.format(newDir, bid)),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.read_distribution.log'.format(bid))),
        shell=True
      )
    
    '''
      codon periodicity
        annotation=/glusterfs/users/ashieh/annotations/hg19_ccds_exons_plus_start100.bed

        bedtools intersect -a {annotation} -b {uniq.bam} -s -bed -wa -wb > intersect_start100
        awk -v OFS='\t' '{print ($2-($14+100))}' ${filename}_intersect_start100.bed
         | sort | uniq -c > ${filename}_relative_pos_aggregate.table
    ''' 

    # bedtools intersect -a {annotation} -b {uniq.bam} -s -bed -wa -wb > intersect_start100

    for bid in bid_list:
      newDir = new_dir(outputDir, bid)
      bedtools.run(
        Parameter('intersect'),
        Parameter('-a {}'.format(pathTo_hg19_bed_start100)),
        Parameter('-b {}/{}.uniq_sorted.bam'.format(newDir, bid)),
        Parameter('-s'),
        Parameter('-bed'),
        Parameter('-wa'),
        Parameter('-wb'),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.intersect_start100.bed'.format(bid))),
        shell=True
      )
      awk.run(
        Parameter('-v'),
        Parameter("OFS='\\t'"),
        Parameter('{print ($8-($2+100))}'),
        Parameter('{}/{}.intersect_start100.bed'.format(newDir, bid)),
        Pipe(
          sort.pipe(
            Pipe(
              uniq.pipe(
                Parameter('-c'),
                Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
                  '{}_relative_pos_aggregate.table'.format(bid)))
              )
            )  
          )
        )
      )

    for bid in bid_list:
      newDir = new_dir(outputDir, bid)
      rpaFile = open('{dir}/{bid}_relative_pos_aggregate.table'.format(dir=newDir, bid=bid), 'rb')
      myDict = {}

      for i in range(-30,31):
        myDict[i] = 0

      for line in rpaFile:
        Frequency, start = line.strip().split(' ')
        if int(start) >= -30 and int(start) <= 30:
          print start
          myDict[int(start)] = Frequency

      # print times

      freqs = []
      starts = []
      for i in range(-30, 31):
        starts.append(i)
        freqs.append(myDict[i])

      # print freqs


      fig, ax = plt.subplots()
      # plt.set_title('{} codon periodicity'.format(bid))
      plt.xlabel("-30 to 30 relative position")
      plt.ylabel("Frequency")
      plt.bar(starts, freqs)
      fig.savefig('{dir}/{bid}_codon_periodicity_plot.png'.format(dir=newDir, bid=bid))
    

    '''
    Picard tools

    java -jar picard.jar CollectMultipleMetrics 
    I=2017-221.uniq_sorted.bam 
    O= multiple_metrics 
    R=GRCh37.p13.genome.fa

    java -jar picard.jar CollectGcBiasMetrics
    I= .uniq
    O=gc_bias_metrics.txt 
    CHART=gc_bias_metrics.pdf 
    S=summary_metrics.txt 
    R=reference_sequence.fasta

    java -jar picard.jar CollectRnaSeqMetrics
    I=input.bam 
    O=output.RNA_Metrics 
    REF_FLAT=ref_flat.txt 
    STRAND=FIRST_READ_TRANSCRIPTION_STRAND

    java -jar picard.jar MarkDuplicates
    I=input.bam 
    O=marked_duplicates.bam 
    M=marked_dup_metrics.txt
    ASSUME_SORTED=true
    '''

    for bid in bid_list:
      newDir = new_dir(outputDir, bid)

      picard.run(
        Parameter('CollectMultipleMetrics'),
        Parameter('I={}/{}.uniq_sorted.bam'.format(newDir, bid)),     # input
        Parameter('O={}/{}.multiple_metrics'.format(newDir, bid)),    # output
        Parameter('R={}'.format(pathTo_genomeFasta))                  # genomeReference
      )

      picard.run(
        Parameter('CollectGcBiasMetrics'),
        Parameter('I={}/{}.uniq_sorted.bam'.format(newDir, bid)),          # input
        Parameter('O={}/{}.gc_bias_metrics'.format(newDir, bid)),           # output
        Parameter('CHART={}/{}.gc_bias_metrics.pdf'.format(newDir, bid)),   # chart
        Parameter('S={}/{}.summary_metrics'.format(newDir, bid)),           # summary metrics
        Parameter('R={}'.format(pathTo_genomeFasta))                        # genome reference
      )

      picard.run(
        Parameter('CollectRnaSeqMetrics'),
        Parameter('I={}/{}.uniq_sorted.bam'.format(newDir, bid)),     # input
        Parameter('O={}/{}.RNA_Metrics'.format(newDir, bid)),         # output
        Parameter('REF_FLAT={}/{}'.format(newDir, bid)),              # ref_flat
        Parameter('STRAND=FIRST_READ_TRANSCRIPTION_STRAND')           # strandedness
      )

      picard.run(
        Parameter('MarkDuplicates'),
        Parameter('I={}/{}.uniq_sorted.bam'.format(newDir, bid)),       # input
        Parameter('O={}/{}.marked_duplicates.bam'.format(newDir, bid)), # output
        Parameter('M={}/{}.marked_dup_metrics.txt'),                    # marked dup metrics
        Parameter('ASSUME_SORTED=true')                                # sorted
      )


    '''
    subread: featureCounts

      featureCounts -a /path_to_gtf/gencode.v19.annotation.gtf -o <bid>.featureCounts <bid>.uniq_sorted.bam
    '''

    for bid in bid_list:
      newDir = new_dir(outputDir, bid)
      featureCounts.run(
        Parameter('-a', '{}'.format(pathToGtf)),                    # gtf
        Parameter('-s', '1'),                                       # strand-specific read counting 
        Parameter('-o', '{}/{}.featureCounts'.format(newDir, bid)), # output
        Parameter('{}/{}.uniq_sorted.bam'.format(newDir, bid))      # input
      )

    # '''
    # FastQC

    #   fastqc --outdir=/path_to/<bid>/ /path_to_fastq/<bid>.fastq.gz
    # '''

    # for fastqlib in fastqFiles:
    #   fastq, bid = fastqlib.split(':')
    #   newDir = new_dir(outputDir, bid)
    #   fastQC.run(
    #     Parameter('--outdir={}'.format(newDir)),   # output
    #     Parameter('--t', numThreads),           
    #     Parameter(fastq)                           # input
    #   )
