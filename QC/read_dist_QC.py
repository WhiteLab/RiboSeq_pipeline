# Author: Thomas Goodman 

import argparse
import sys
import subprocess
import csv

parser = argparse.ArgumentParser()

parser.add_argument('--samtoolsPath')
parser.add_argument('--libs', nargs='*')

args = vars(parser.parse_args())

samtools = args['samtoolsPath']
libs = args['libs']

# Dictionary to keep track of everything
lib_dict = {}

for lib in libs:
  lib_dict[lib] = {}

# Read_distribution

# Fill out data
# read distribution
read_distList = []
for lib in libs:
  read_dist_file = open('{l}/{l}.read_distribution.log'.format(l=lib), 'rb')
  i = 0
  #sys.stdout.write(lib + '\t')
  nextLine = []
  nextLine.append(lib)
  # Enumerate might have been good here. LUL
  for line in read_dist_file:
    splitLine = line.split()
    if i < 2:
      #sys.stdout.write(splitLine[2] + '\t')
      nextLine.append(splitLine[2])
    elif i == 2:
      #sys.stdout.write(splitLine[3] + '\t')
      nextLine.append(splitLine[3])
    elif i > 4:
      if i < 9 or i == 11 or i == 14:
        #sys.stdout.write(('\t').join(splitLine[1:]) + '\t')
        for x in splitLine[1:]:
          nextLine.append(x)
    i += 1
  #sys.stdout.write('\n')
  read_distList.append(nextLine)

# Add up all tags, and find percent of tags in each region
percentageRegions = []

# Tag spot in array
CDS = 6        #G
UTR5 = 9       #J
UTR3 = 12      #M
INTRONS = 15   #P
TSS = 18       #S
TES = 21       #V

featureList = [[CDS,0,"CDS"],[UTR5,0,"5UTR"],[UTR3,0,"3UTR"],[INTRONS,0,"INTRONS"],[TSS,0,"TSS"],[TES,0,"TES"]]

read_distList.sort()

sys.stdout.write("bid,totalTags,CDS,5UTR,3UTR,INTRONS,TSS,TES")
for sample in read_distList:
  totalTags = 0
  for feature in featureList:
    totalTags += float(sample[feature[0]])
  sys.stdout.write("\n{},".format(sample[0]))
  sys.stdout.write('{},'.format(float(sample[feature[0]])/totalTags))
  for feature in featureList:
    sys.stdout.write('{},'.format(float(sample[feature[0]])/totalTags))
    feature[1] += float(float(sample[feature[0]])/totalTags)

sys.stdout.write("\nfeature,Tags_percentage\n")
is100 = 0
for feature in featureList:
  is100 += feature[1]/len(read_distList)
  sys.stdout.write("{},{}\n".format(feature[2],str(feature[1]/len(read_distList))))

sys.stdout.write("{}\n".format(str(is100)))

for sample in read_distList:
  sys.stdout.write("{}\n".format(','.join(sample)))

