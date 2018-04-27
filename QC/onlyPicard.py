'''
RiboSeq Picard QC Pipeline
Author: Thomas Goodman
Pipeline created using Chunky pipes

Pipeline to process Picard QC for Rna/Ribo seq bams 


'''

import matplotlib
matplotlib.use('Agg')
import os
import multiprocessing
import subprocess
import re
import time
from chunkypipes.components import Software, Parameter, Redirect, Pipe, BasePipeline

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
      'picard':{
        'path'        : 'Full path and command to run picard ex: java -Xms5g -Xmx6g -jar /mnt/picard/dist/picard.jar \n \tNote: -Xms5g -Xmx6g may need to be changed depending on your memory needs ',
        'genomeFasta' : 'Full path to the genome fasta reference file'
      }
    }

  # Required input for pipeline software
  # Next time put all these paths to annotations in configure
  def add_pipeline_args(self, parser):
    parser.add_argument('--bam:lib', required=True, nargs='*', help='Fastq input for pipeline:library name(prefix for files)')
    parser.add_argument('--output', required=True, help='Where pipeline output should go')
    

    # chunky run RiboSeq_pipe.py --fastqs 
    #  /mnt/cinder/thomas/RiboSeq/Lane5/AWS-3_S3_L005_R1_001.fastq.gz
    #  --output /mnt/cinder/thomas/RiboSeq/test --threads

  def run_pipeline(self, pipeline_args, pipeline_config):
    # create variables from parser if wanted
    bamFiles = pipeline_args['bam:lib']
    outputDir = pipeline_args['output']

    # Create output directory
    subprocess.call(['mkdir', outputDir])

    # Software
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
    pathTo_genomeFasta = pipeline_config['picard']['genomeFasta']


    # Keep track of Bids in pipeline

    # bid_list = []
    # bam_list = []
    # for bamLib in bamFiles:
    #   bid_list.append(bamLib.split(':')[1])
    #   bam_list.append(bamLib.split(':')[0])


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

    for bamLib in bamFiles:
      bam, bid = bamLib.split(':')
      newDir = new_dir(outputDir, bid)
      subprocess.call(['mkdir', newDir])

      # consider multithreading?

      picard.run(
        Parameter('CollectMultipleMetrics'),
        Parameter('I={}'.format(bam)),     # input
        Parameter('O={}/{}.multiple_metrics'.format(newDir, bid)),    # output
        Parameter('R={}'.format(pathTo_genomeFasta))                  # genomeReference
      )

      picard.run(
        Parameter('CollectGcBiasMetrics'),
        Parameter('I={}'.format(bam)),          # input
        Parameter('O={}/{}.gc_bias_metrics'.format(newDir, bid)),           # output
        Parameter('CHART={}/{}.gc_bias_metrics.pdf'.format(newDir, bid)),   # chart
        Parameter('S={}/{}.summary_metrics'.format(newDir, bid)),           # summary metrics
        Parameter('R={}'.format(pathTo_genomeFasta))                        # genome reference
      )

      picard.run(
        Parameter('CollectRnaSeqMetrics'),
        Parameter('I={}'.format(bam)),     # input
        Parameter('O={}/{}.RNA_Metrics'.format(newDir, bid)),         # output
        Parameter('REF_FLAT={}/{}'.format(newDir, bid)),              # ref_flat
        Parameter('STRAND=FIRST_READ_TRANSCRIPTION_STRAND')           # strandedness
      )

      picard.run(
        Parameter('MarkDuplicates'),
        Parameter('I={}'.format(bam)),       # input
        Parameter('O={}/{}.marked_duplicates.bam'.format(newDir, bid)), # output
        Parameter('M={}/{}.marked_dup_metrics.txt'.format(new, bid)),                    # marked dup metrics
        Parameter('TMP_DIR={}'.format(newDir)),
        Parameter('ASSUME_SORTED=true')                                # sorted
        Parameter('VALIDATION_STRINGENCY=LENIENT'),
        Redirect(stream=Redirect.BOTH, dest=os.path.join(newDir, 'mark_dup.log'))
      )
