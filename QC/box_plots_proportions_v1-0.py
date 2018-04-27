'''
Author: Thomas Goodman

This program creates boxplots out of data in batches

TSV input format should be:

Bid averageReadLength totalReads  readsWritten  totalBasesWritten rtsReads  perc_rts  totalKeptReads  perc_kept primaryMappedReads  multiMappedReads  perc_multi  uniqueMappedReads perc_uniq
BN  221-222 g                  
2017-221  28.58787339 261413021 154599597 4487140371  98618295  0.6378949034  55981302  0.3621050966  51145848  7408282 0.0479191547  43737566  0.1673121172
2017-222  28.13942448 268064617 142959424 4101511402  103809826 0.7261488826  39149598  0.2738511174  33032328  7966543 0.055725903 25065785  0.0935065033
BN  543-544 #1da7d1                    
2017-543  28.55728576 251538680 147866979 4291229021  109869511 0.7430293886  37997468  0.2569706114  32343870  5585964 0.0377769536  26757906  0.1063769039
2017-544  28.25695512 251857442 130464343 3788336616  99833718  0.7652184168  30630625  0.2347815832  25358343  4788864 0.0367063053  20569479  0.0816711185

available colors

b: blue
g: green
r: red
c: cyan
m: magenta
y: yellow
k: black
w: white

or use html hex #1da7d1
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import sys
import subprocess
import os
import random

parser = argparse.ArgumentParser()
parser.add_argument('--tsv', help='tsv in form of column Names\nBN x-x color\nbid data\nBN y-y\nbid data\nand so on...')
parser.add_argument('--out', help='Name of a new dir to save plots')

sampleData = []

samplesRanges = []
sampleNames = []
sampleColors = []
bids = []
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
color = random.choice(colors)

args = vars(parser.parse_args())
sampleSheet = args['tsv']
output = args['out']

if not os.path.exists(output):
  subprocess.call(['mkdir', output])

first = True
with open(sampleSheet, 'rb') as ss:
  for line in sampleSheet:
    bid, proportion

# print len(sampleData)
# print samplesRanges

flierprops = dict(marker='o', markerfacecolor='m', markersize=6.5,
                  linestyle='none', markeredgecolor='m')
whiskerprops = dict(linestyle='solid',color='k')

size = 50

jj = []
for i in range(1, len(sampleNames)+1):
  jj.append(i * 10)

for i, header in enumerate(headers[1:]):
  plotData = []

  for y, x in enumerate(sampleData):
    plotData.append(x[i])

  fig, ax = plt.subplots(figsize=(12, 12))
  
  # bplot1 = plt.boxplot([plotData[0:6], plotData[7:15], plotData[16:25], plotData[26:35], plotData[36:47], plotData[48:]], meanline=True, widths=4,
  #             positions=jj,labels=sampleNames, flierprops=flierprops, whiskerprops=whiskerprops, patch_artist=True)

  # batches = [plotData[0:6], plotData[7:15], plotData[16:25], plotData[26:35], plotData[36:47], plotData[48:]]
  # sampleNames= ['221 - 232', '543 - 554', '787 - 797', '798 - 809', '900 - 911', '981 - 992']

  # Set range of batches
  batches = []
  for batch in samplesRanges:
    batches.append(plotData[batch[0]:batch[1]])

  bplot1 = plt.boxplot(batches, meanline=True, widths=4, positions=jj,labels=sampleNames, 
                       flierprops=flierprops, whiskerprops=whiskerprops, patch_artist=True)


  ax.set_xlim([-1,len(plotData)+5])
  for patch, color in zip(bplot1['boxes'], sampleColors):
    patch.set(color=color)
    patch.set(facecolor='none')

  # myLegend = []
  # for label, color in samples.iteritems():
  #   nextLegend = mpatches.Patch(color=color, label=label)
  #   myLegend.append(nextLegend)
  
  plt.title(header)
  plt.xlabel("Batch")

  plt.tight_layout()
  fig.savefig('{o}/{h}'.format(o=output, h=header))

