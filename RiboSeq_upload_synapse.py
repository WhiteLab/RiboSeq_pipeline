from copy import copy
import multiprocessing as mp
import synapseclient
from synapseclient import File
import argparse, sys

# Ribo-seq FASTQ_FOLDER_SYNAPSE = 'syn00000000'

# RNA-seq format
#   PEC_BrainGVEX_UIC-UChicago_FC_mRNA_HiSeq2000_2014-2194_AC5K6UACXX_6_1.fastq.gz

# Ribo-seq format
#   PEC_BrainGVEX_UIC-UChicago_FC_mRNA_HiSeq2000_2017-1_AC5K6UACXX_6_1.fastq.gz

parser = argparse.ArgumentParser()
parser.add_argument('--fastqs', nargs='*')
parser.add_argument('--synID')
parser.add_argument('--threads')
args = vars(parser.parse_args())

fastqs = args['fastqs']
FASTQ_FOLDER_SYNAPSE = args['synID']
NUM_PROCESSORS = int(args['threads'])

# fastqFilepath = '/filestore/RSYNC_SHARE/PEC/UC-UIC_ATACseq_Stanley/RAW'

platformDict = {'K00242':'HiSeq4000',
                'SN1070':'HiSeq2500',
                '700819F':'HiSeq2500'}

def upload_fastq(fastq, bid, annot, theFilepath):
    print 'Uploading {}'.format(fastq)
    barcode = '_'.join(fastq.split('_')[4:7])
    platform = annot['platform']
    filename = synapse_filename_template.format(bid=bid,platform=platform,barcode=barcode)
    fastq_synapse = File(theFilepath, name=filename, parent=FASTQ_FOLDER_SYNAPSE)

    # Deprecated. Will error out trying to upload_fastq
    #   fastq_synapse.properties.fileNameOverride = filename
    
    fastq_synapse.annotations = annot
    print "Uploading"
    syn.store(fastq_synapse)

default_annotations = {
    'fileType': 'fastq',
    'assay': '',
    'grant': '', # grant ID
    'runType': '',
    'study': '', 
    'organism': '', 
    'consortium': '',
    'organ': '',
    'PI': '',
    'tissueType': '',
    'tissueTypeAbrv': '',
    'readLength': '', 
    'group': 'UIC-UChicago'
}

synapse_filename_template = 'PEC_BrainGVEX_UIC-UChicago_FC_RiboSeq_{platform}_{bid}_{barcode}.fastq.gz'
    
# Start upload to Synapse
# TODO Multiprocess this
print 'Logging into Synapse'
syn = synapseclient.Synapse()
syn.login('email', 'pass')

upload_pool = mp.Pool(processes=NUM_PROCESSORS)

print 'Iterating through samples'
for myFile in fastqs:
    fastq = myFile.split('/')[-1]
    annot = copy(default_annotations)
    bid = fastq.split('_')[0]
    platform = platformDict[fastq.split('_')[2]]
    annot.update({
        'platform': platform
    })
    upload_pool.apply_async(upload_fastq, args=(fastq, bid, annot, myFile))
    # upload_fastq(fastq, bid, annot, myFile)

upload_pool.close()
upload_pool.join()


