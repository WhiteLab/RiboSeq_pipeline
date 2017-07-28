# Author: Thomas Goodman 28 July 2017
#
# This program should be placed and ran in a directory where a the RiboSeq pipeline output all the librarys
# Issues?
# could possibly be changed to just take in each file, but that's a lot more work than I want to do
# this also assumes the naming conventions have not changed. 
# Would be easy to add what the extensions for each file type will be after the library name

import argparse
import sys
import subprocess
import csv

parser = argparse.ArgumentParser()

parser.add_argument('--samtoolsPath')
parser.add_argument('--libs', nargs='*')
# parser.add_argument('--uniqReads', nargs='*')
# parser.add_argument('--caLogs', nargs='*')
# parser.add_argument('--bt2L', nargs='*')
# parser.add_argument('--starL', nargs='*')

args = vars(parser.parse_args())

samtools = args['samtoolsPath']
libs = args['libs']
# caLogs = args['caLogs']
# uniqFiles = args['uniqReads']
# bt2Logs = args['bt2L']
# starLogs = args['starL']

# Dictionary to keep track of everything
lib_dict = {}

for lib in libs:
  lib_dict[lib] = {}

# Average Read Length

for lib in libs:

  uniqFile = '{l}/{l}.uniq_sorted.bam'.format(l=lib)

  subprocess.call([samtools, 'view', uniqFile, '-o', '{}.txt'.format(uniqFile)])

  uniq_sorted = open('{}.txt'.format(uniqFile), 'rb')

  numLines = 0
  totalLength = 0
  for line in uniq_sorted:
    read = line.split('\t')[9]
    totalLength += len(read)
    numLines += 1

  averageReadLength = float(totalLength) / numLines

  lib_dict[lib]['arl'] = averageReadLength

  subprocess.call(['rm', '{}.txt'.format(uniqFile)])

# Cutadapt Logs

  # firstLine = ['totalReads', 'readsWritten', 'totalBasesWritten']

  cutadaptLog = '{l}/{l}.cutadapt.summary.log'.format(l=lib)

  with open(cutadaptLog, 'rb') as aLog:
    totalReads = ''
    readsWritten = ''
    totalBasesWritten = ''
    i = 0
    for line in aLog:
      if i == 7:
        totalReads = line.strip().split(':')[-1]
      if i == 10:
        readsWritten = line.strip().split(':')[-1]
      if i == 13:
        totalBasesWritten = line.strip().split(':')[-1]
        break
      i += 1
      
  lib_dict[lib]['tr'] = totalReads.replace(',', '')
  lib_dict[lib]['rw'] = readsWritten.replace(',', '')
  lib_dict[lib]['tbw'] = totalBasesWritten.replace(',', '')


  # Bowtie  STAR

  '''
  totalTrimmedReads  = total reads of trimmed fq file (bowtie2 input)
  rtsReads           = number of reads mapped to rts sequence and filtered (can be found in bowtie2 mapping statistics, = # reads mapped 1 time + # reads mapped > 1 time)
  totalKeptReads     = number of reads aligned 0 times from bowtie2 filtering step (also = totalReads - rtsReads, also = tophat Input)
  primaryMappedReads = tophat alignment result: can be found in align_summary.txt file of tophat output (= number of 'Mapped' reads in align_summary.txt) 
  multiMappedReads   = number of reads with multiple alignments in align_summary.txt
  uniqMappedReads    = primaryMappedReads - multiMappedReads (also = number of reads remaining after extracting 'NH:i:1' reads)
  '''

  nextSectionLine = ['arl','tr','rw','tbw','ttr','rtsr','tkr','pmr','mmr','umr']

  # Bowtie

  with open('{l}/{l}.bowtie2.log2'.format(l=lib), 'rb') as btlog:
    btlogReader = csv.reader(btlog, dialect='excel', delimiter='\t')
    i = 0
    for row in btlogReader:
      theRow = row[0]
      splitRow = theRow.split()
      if i == 0:
        totalTrimmedReads = splitRow[0].replace(',', '')
      elif i == 2:
        totalKeptReads = splitRow[0].replace(',', '')
      elif i == 3:
        rtsReads = int(splitRow[0])
      elif i == 4:
        rtsReads += int(splitRow[0])
        break
      i += 1

  # STAR

  with open('{l}/{l}_Log.final.out'.format(l=lib), 'rb') as starLog:
    starLogReader = csv.reader(starLog, dialect='excel', delimiter='\t')
    i = 0
    for row in starLogReader:
      if i == 8:
        uniqMappedReads = int(row[1])
      if i == 23:
        multiMappedReads = int(row[1])
        break
      i += 1
    primaryMappedReads = int(multiMappedReads) + int(uniqMappedReads)

  nextRow = [totalTrimmedReads, rtsReads, totalKeptReads, primaryMappedReads, multiMappedReads, uniqMappedReads]

  for i, abv in enumerate(nextSectionLine[4:]):
    lib_dict[lib][abv] = nextRow[i]

firstLine = ['Bid', 'averageReadLength', 'totalReads', 'readsWritten', 'totalBasesWritten','totalTrimmedReads','rtsReads', 'perc_rts', 
             'totalKeptReads', 'perc_kept','primaryMappedReads','multiMappedReads', 'perc_multi', 'uniqueMappedReads', 'perc_uniq']

i = 0

for lib in lib_dict:
  if i == 0:
    i+=1
    sys.stdout.write(','.join(firstLine) + '\n')
    continue

  sys.stdout.write(lib + ',')
  for key in nextSectionLine:
    sys.stdout.write(str(lib_dict[lib][key]) + ',')
    if key == 'rtsr':
      perc_rts = float(lib_dict[lib][key]) / float(lib_dict[lib]['ttr'])
      sys.stdout.write(str(perc_rts) + ',')
    if key == 'tkr':
      perc_tkr = float(lib_dict[lib][key]) / float(lib_dict[lib]['ttr'])
      sys.stdout.write(str(perc_tkr) + ',')
    if key == 'mmr':
      perc_mmr = float(lib_dict[lib][key]) / float(lib_dict[lib]['ttr'])
      sys.stdout.write(str(perc_mmr) + ',')
    if key == 'umr':
      perc_umr = float(lib_dict[lib][key]) / float(lib_dict[lib]['tr'])      
      sys.stdout.write(str(perc_umr) + '\n')