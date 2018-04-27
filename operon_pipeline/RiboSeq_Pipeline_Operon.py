# RiboSeq_pipeline_Operon

'''
RiboSeq Pipeline using operon https://github.com/djf604/operon
Author: Thomas Goodman

Pipeline to process RiboSeq fastqs for data analysis Alignment + QC

Pipeline created using Operon

using python3
'''

## Get versions of software used 
## Fix picard


import time, subprocess, os

from operon.components import ParslPipeline, Software, Parameter, Redirect, Data, CodeBlock, CondaPackage, Pipe

defaultThreads = 1

def relative_pos_table(intersect_bedtools_filepath, relativePos_filepath):
  from collections import Counter
  start100_file = open(intersect_bedtools_filepath, 'r')
  relativePos_file = open(relativePos_filepath, 'w')
  distanceList = []
  for line in start100_file:
    splitLine = line.split('\t')
    # Really is relative start
    if len(splitLine) >= 7:
      distance = int(splitLine[7]) - (int(splitLine[1]) + 100)
      distanceList.append(distance)
  distanceList.sort()
  distanceCounting = Counter(distanceList)
  for key, value in distanceCounting.items():
    relativePos_file.write("{}\t{}\n".format(value,key))
  start100_file.close()

# Create plots for showing codon periodicity
def create_codon_periodicity(relativePos_filepath, codon_periodicity_filepath):
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.mlab as mlab
  import matplotlib.pyplot as plt
  rpaFile = open(relativePos_filepath, 'r')
  myDict = {}

  for i in range(-30,31):
    myDict[i] = 0

  for line in rpaFile:
    Frequency, start = line.strip().split('\t')
    if int(start) >= -30 and int(start) <= 30:
      # print start
      myDict[int(start)] = Frequency

  # Change to log scaling?

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
  fig.savefig(codon_periodicity_filepath)

class Pipeline(ParslPipeline):
  def description(self):
    return 'Pipeline to process RiboSeq fastqs for data analysis Alignment + QC'

  def dependencies(self):
    depend = []
    # depend.append('cutadapt')
    # depend.append('RSeQC') # for read_distribution.py
    return depend

  def conda(self):
    return {
      'packages': [
        CondaPackage(tag='samtools=1.6', config_key='samtools'),
        # CondaPackage(tag='picard=2.9.2', config_key='picard'),
        CondaPackage(tag='bedtools=2.25.0', config_key='bedtools'),
        CondaPackage(tag='star=2.5.3a', config_key='star', executable_path='bin/STAR'),
        CondaPackage(tag='bowtie2=2.3.3.1', config_key='bowtie2'),
        CondaPackage(tag='fastqc', config_key='FastQC'),
        CondaPackage(tag='subread', config_key='featureCounts', executable_path='bin/featureCounts'), # feature_counts
        CondaPackage(tag='cutadapt', config_key='cutadapt'),
        CondaPackage(tag='rseqc', config_key='read_distribution', executable_path='bin/read_distribution.py')
      ]
  }

  def arguments(self, parser):
    parser.add_argument('--fastq:lib', required=True, nargs='*', help='Fastq input for pipeline:library name(prefix for files)')
    parser.add_argument('--output', required=True, help='Where pipeline output should go')
    parser.add_argument('--adapter', default='AGATCGGAAGAGCACACGTCT', help='Adapter sequence for trimming')
    parser.add_argument('--threads', default=defaultThreads, help='Threads to be used for multi-threaded programs. Default is 8')

  def configuration(self):
    return{
      'cutadapt':{
        'path': 'Full path to cutadapt executable'
      },
      'star':{
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
        'path' : 'Full path to featureCounts'
      },
      'picard':{
        'path'      : 'Full path to picard jar ex: /mnt/picard/dist/picard.jar',
        'genomeFasta' : 'Full path to the genome fasta reference file',
        'refFlat'   : 'Full path to the refFlat.txt reference file'
      }
    }

  def pipeline(self, pipeline_args, pipeline_config):
    # chunky run RiboSeq_pipe.py --fastqs 
    #  /mnt/cinder/thomas/RiboSeq/Lane5/AWS-3_S3_L005_R1_001.fastq.gz
    #  --output /mnt/cinder/thomas/RiboSeq/test --threads

    # create variables from parser if wanted
    fastqFiles = pipeline_args['fastq:lib']
    outputDir = pipeline_args['output']
    adapter = pipeline_args['adapter']
    numThreads = pipeline_args['threads']

    # Create output directory
    subprocess.call(['mkdir', outputDir])

    # Software
    cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
    star = Software('star', pipeline_config['star']['path'])
    bedtools = Software('bedtools', pipeline_config['bedtools']['path'])
    bowtie2 = Software('bowtie2', pipeline_config['bowtie2']['path'])
    samtools = Software('samtools', pipeline_config['samtools']['path'])
    samtools_header = Software('samtools', pipeline_config['samtools']['path'])
    samtools_uniq = Software('samtools', pipeline_config['samtools']['path'])
    samtools_sort = Software('samtools sort', pipeline_config['samtools']['path'])
    read_distribution = Software('read_distribution', 
      pipeline_config['read_distribution']['path'])
    featureCounts = Software('featureCounts', pipeline_config['featureCounts']['path'])
    fastQC = Software('FastQC', pipeline_config['FastQC']['path'])
    picard = Software('picard', "java -Xms8g -Xmx9g -jar {}".format(pipeline_config['picard']['path']))

    # Change these to just be done in python script?

    # Common software tools 
    awk = Software('awk', 'awk')
    sort = Software('sort', 'sort')
    uniq = Software('uniq', 'uniq')
    paste = Software('paste', 'paste')
    cat = Software('cat', 'cat')
    grep = Software('grep', 'grep')

    # Directories and Files
    pathToGenomeDir = pipeline_config['star']['genomeDir']
    pathToGenome = pipeline_config['bowtie2']['genome_ref']
    pathToGtf = pipeline_config['star']['GTF_ref']
    pathTo_hg19_bed = pipeline_config['read_distribution']['hg19_bed']
    pathTo_hg19_bed_start100 = pipeline_config['bedtools']['hg19_start100']
    pathTo_grch37_bed = pipeline_config['bedtools']['grch37_bed']
    pathTo_genomeFasta = pipeline_config['picard']['genomeFasta']
    pathTo_ref_flat = pipeline_config['picard']['refFlat'] 


    # bid_list = []
    # for fastqlib in fastqFiles:
    #   bid_list.append(fastqlib.split(':')[-1])

    for fastqlib in fastqFiles:
      fastq, bid = fastqlib.split(':')

      # Make new directories to store data
      newDir = os.path.join(outputDir, bid)
      subprocess.call(['mkdir', newDir])

      trimmed_read_filename = '{}/{}.trimmed.fastq.gz'.format(newDir, bid)

      cutadapt.register(
        Parameter('--quality-base=33'),
        Parameter('--minimum-length=25'),
        Parameter('--discard-untrimmed'),
        Parameter('--output={}'.format(trimmed_read_filename)),
        # Parameter('--cores', numThreads),
        Parameter('-a', adapter),
        Parameter(Data(fastq).as_input()),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.cutadapt.summary.log'.format(bid))),
        extra_outputs=[Data(trimmed_read_filename).as_output(tmp=True)]
      )
      
      bowtie2.register(
        Parameter('--seedlen=23'),
        Parameter('--threads', numThreads),
        Parameter('--un-gz', Data('{}/{}.filtered.fastq.gz'.format(newDir, bid)).as_output()),
        Parameter('-x', Data(pathToGenome).as_input()), # Path to rtsRNA_seqs files
        Parameter('-U', Data(trimmed_read_filename).as_input()),
        Parameter('-S'),
        Parameter(Data('{}/{}.rts.sam'.format(newDir, bid)).as_output(tmp=True)),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.bowtie2.log'.format(bid))),
        Redirect(stream=Redirect.STDERR, dest=os.path.join(newDir, 
          '{}.bowtie2.log2'.format(bid)))  
      )

      samtools.register(
        Parameter('view'),
        Parameter('-Sb'),
        Parameter(Data('{}/{}.rts.sam'.format(newDir, bid)).as_input()),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.rts.bam'.format(bid))),
      )

      star.register(
        Parameter('--runThreadN', numThreads), # Change to command line parameter --threads
        Parameter('--sjdbGTFfile', pathToGtf),
        Parameter('--outSAMtype', 'BAM', 'Unsorted'),
        Parameter('--outFileNamePrefix', '{}/{}_'.format(newDir, bid)),
        Parameter('--genomeDir', pathToGenomeDir),
        # Parameter('--genomeLoad', genomeLoad), broken
        Parameter('--readFilesIn', Data('{}/{}.filtered.fastq.gz'.format(newDir, bid)).as_input()),
        Parameter('--readFilesCommand zcat'), # reads gzipped files
        extra_outputs=[Data('{}/{}.Aligned.bam'.format(newDir, bid)).as_output()]
      )

      samtools_header.register(
        Parameter('view'),
        Parameter('-H'),
        Parameter(Data('{}/{}.Aligned.bam'.format(newDir, bid)).as_input()), # star outfile name
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.header.sam'.format(bid))),
        extra_outputs=[Data('{}/{}.header.sam'.format(newDir, bid)).as_output()]
      )

      uniq_bam = '{}/{}.uniq_sorted.bam'.format(newDir, bid)

      samtools_uniq.register(
        Parameter('view'),
        Parameter(Data('{}/{}.Aligned.bam'.format(newDir, bid)).as_input()), # star outfile name
        Pipe(
          grep.prep(
            Parameter('-w'),
            Parameter('NH:i:1'),
            Pipe(
              cat.prep(
                Parameter(os.path.join(newDir, '{}.header.sam'.format(bid)), '-'),
                Pipe(
                  samtools.prep(
                    Parameter('view'),
                    Parameter('-bS', '-'),
                    Pipe(
                      samtools.prep(
                        Parameter('sort'),
                        Parameter('-', '-o', Data(uniq_bam).as_output())
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )

# QC

# Codon_periodicity
      
      intersectBed_filepath = '{}/{}.intersect_start100.bed'.format(newDir, bid)
      relative_pos_table_filepath = '{}/{}_relative_pos_aggregate.table'.format(newDir, bid)

      bedtools.register(
        Parameter('intersect'),
        Parameter('-a {}'.format(pathTo_hg19_bed_start100)),
        Parameter('-b', Data(uniq_bam).as_input()),
        Parameter('-s'),
        Parameter('-bed'),
        Parameter('-wa'),
        Parameter('-wb'),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.intersect_start100.bed'.format(bid))),
        extra_outputs=[Data(intersectBed_filepath).as_output()]
      )
    
      CodeBlock.register(
        func=relative_pos_table,
        args=[],
        kwargs={'intersect_bedtools_filepath' : intersectBed_filepath,
                'relativePos_filepath': relative_pos_table_filepath},
        inputs=[Data(intersectBed_filepath).as_input()],
        outputs=[Data(relative_pos_table_filepath).as_output()]
      )

      codon_periodicity_filepath = '{}/{}_codon_periodicity_plot.png'.format(newDir, bid)

      CodeBlock.register(
        func=create_codon_periodicity,
        args=[],
        kwargs={'relativePos_filepath': relative_pos_table_filepath,
                'codon_periodicity_filepath' : codon_periodicity_filepath},
        inputs=[Data(relative_pos_table_filepath).as_input()],
        outputs=[Data(codon_periodicity_filepath).as_output()]
      )

# Picard

      picard.register(
        Parameter('CollectMultipleMetrics'),
        Parameter('I={}'.format(uniq_bam)),                           # input
        Parameter('O={}/{}.multiple_metrics'.format(newDir, bid)),    # output
        Parameter('R={}'.format(pathTo_genomeFasta)),                 # genomeReference
        extra_inputs=[Data(uniq_bam)]
      )

      picard.register(
        Parameter('CollectGcBiasMetrics'),
        Parameter('I={}'.format(uniq_bam)),       
        Parameter('O={}/{}.gc_bias_metrics'.format(newDir, bid)),           # output
        Parameter('CHART={}/{}.gc_bias_metrics.pdf'.format(newDir, bid)),   # chart
        Parameter('S={}/{}.summary_metrics'.format(newDir, bid)),           # summary metrics
        Parameter('R={}'.format(pathTo_genomeFasta)),                       # genome reference
        extra_inputs=[Data(uniq_bam)]
      )

      picard.register(
        Parameter('CollectRnaSeqMetrics'),
        Parameter('I={}'.format(uniq_bam)),    
        Parameter('O={}/{}.RNA_Metrics'.format(newDir, bid)),          # output
        Parameter('REF_FLAT={}'.format('{}'.format(pathTo_ref_flat))), # ref_flat
        Parameter('STRAND=FIRST_READ_TRANSCRIPTION_STRAND'),           # strandedness
        extra_inputs=[Data(uniq_bam)]
      )

      picard.register(
        Parameter('MarkDuplicates'),
        Parameter('I={}'.format(uniq_bam)),       
        Parameter('O={}/{}.marked_duplicates.bam'.format(newDir, bid)),  # output
        Parameter('M={}/{}.marked_dup_metrics.txt'.format(newDir, bid)), # marked dup metrics
        Parameter('ASSUME_SORTED=true'),                                 # It is sorted
        extra_inputs=[Data(uniq_bam)]
      )

# featureCounts
      
      featureCounts.register(
        Parameter('-a', Data('{}'.format(pathToGtf)).as_input()),   # gtf
        Parameter('-s', '1'),                                       # strand-specific read counting 
        Parameter('-o', '{}/{}.featureCounts'.format(newDir, bid)), # output
        Parameter(Data(uniq_bam).as_input())                        # input
      )

# fastQC

      fastQC.register(
        Parameter('--outdir={}'.format(newDir)),   # output
        Parameter('--t', numThreads),           
        Parameter(Data(fastq).as_input())
      )

# read_distribution
      
      read_distribution.register(
        Parameter('-r'),
        Parameter(Data(pathTo_hg19_bed).as_input()),
        Parameter('-i'),
        Parameter(Data(uniq_bam).as_input()),
        Redirect(stream=Redirect.STDOUT, dest=os.path.join(newDir, 
          '{}.read_distribution.log'.format(bid)))
      )
    