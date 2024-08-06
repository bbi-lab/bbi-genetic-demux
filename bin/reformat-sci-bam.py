#!/usr/bin/env python3

import argparse
import glob
import pysam

parser = argparse.ArgumentParser(description="convert sci data to a format souporcell can use")
parser.add_argument("-i", "--input_dir", required=True, help="directory with sci bams")
parser.add_argument("-r", "--run", required= True, help = "Run number to append")
parser.add_argument("-o", "--output_prefix", required=True, help="outputs bam and barcode files with this prefix")
args = parser.parse_args()

bamfiles = glob.glob(args.input_dir + "/*.bam")
print(bamfiles)
assert len(bamfiles) > 0, "I don't see bam files in directory" + args.input_dir
template = pysam.AlignmentFile(bamfiles[0])
outbam = pysam.AlignmentFile(args.output_prefix + "_" + args.run + ".bam", 'wb', template=template)
cell_barcode = 0

for (cell_barcode, bamfile) in enumerate(bamfiles):
    inbam = pysam.AlignmentFile(bamfile)
    for read in inbam:
        barcode = '_'.join(read.query_name.split('|')[2:5]) + "_" + args.run
        umi = read.query_name.split('|')[-1]
        read.set_tag("UB", str(umi))
        read.set_tag("CB", str(barcode))
        outbam.write(read)

