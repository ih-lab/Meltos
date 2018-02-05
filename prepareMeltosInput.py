import sys
import os
import numpy as np
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="file path to SV regions in bed file")
parser.add_argument("-o", "--output", help="output file name")
parser.add_argument("-b", "--bams", help="list of bam files corresponding to tumor regions")
parser.add_argument("-f", "--frag", help="comma separated list of average fragment lengths for each bam")
parser.add_argument("-vcf", help="comma separated list of vcf tardis files with deletions in each bam")
parser.add_argument("-vcfbed", help="comma separated list of bed files associated with tardis vcfs")

args = parser.parse_args()

input_bed = args.input
outputfile = args.output
bamlist = args.bams.split(',')
frag_length = args.frag.split(',')
vcflist = args.vcf.split(',')
bedlist = args.vcfbed.split(',')

samples = []

for each in bamlist:
	temp1 = os.path.basename(each)
	temp2 = temp1.replace(".bam", "")
	samples.append(temp2)


df1 = pd.read_csv(input_bed, sep="\t", header=None, names=['Chr', 'Start', 'End'])

#samples = ['PD9770A','PD9770C','PD9770D']

df1.Start = df1.Start + 1

#import VCF files
s0_VCF = pd.read_csv(vcflist[0], sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])
s1_VCF = pd.read_csv(vcflist[1], sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])
s2_VCF = pd.read_csv(vcflist[2], sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])

#s0_VCF = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/a/PD9770A_deletions.vcf", sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])
#s1_VCF = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/c/PD9770C_deletions.vcf", sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])
#s2_VCF = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/d/PD9770D_deletions.vcf", sep="\t", index_col=False, skiprows = 40, usecols = [0, 1, 7, 9], names = ['#CHROM', 'POS', 'INFO', 'Sample'])


#s0SVs_df = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/a/PD9770A_deletions.bed", sep="\t", header=None, names=['Chr', 'Start', 'End'])
#s1SVs_df = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/c/PD9770C_deletions.bed", sep="\t", header=None, names=['Chr', 'Start', 'End'])
#s2SVs_df = pd.read_csv("/pbtech_mounts/ihlab_store_athena/cnr38/BreastCancer/PD9770/d/PD9770D_deletions.bed", sep="\t", header=None, names=['Chr', 'Start', 'End'])

s0SVs_df = pd.read_csv(bedlist[0], sep="\t", header=None, names=['Chr', 'Start', 'End'])
s1SVs_df = pd.read_csv(bedlist[1], sep="\t", header=None, names=['Chr', 'Start', 'End'])
s2SVs_df = pd.read_csv(bedlist[2], sep="\t", header=None, names=['Chr', 'Start', 'End'])

s0SVs_df.Start = s0SVs_df.Start + 1
s1SVs_df.Start = s1SVs_df.Start + 1
s2SVs_df.Start = s2SVs_df.Start + 1
#s3SVs_df.Start = s3SVs_df.Start + 1
#s6SVs_df.Start = s6SVs_df.Start + 1
#s9SVs_df.Start = s9SVs_df.Start + 1
#s10SVs_df.Start = s10SVs_df.Start + 1


#for each in samples:
allSVs = df1.iloc[:,0].map(str) + ':' + df1.iloc[:,1].map(str) + '-' + df1.iloc[:,2].map(str) 
allSVs = np.array(allSVs)

s0SVs = s0SVs_df.iloc[:,0] + ':' + s0SVs_df.iloc[:,1].map(str) + '-' + s0SVs_df.iloc[:,2].map(str) 
s0SVs = np.array(s0SVs)
s1SVs = s1SVs_df.iloc[:,0] + ':' + s1SVs_df.iloc[:,1].map(str) + '-' + s1SVs_df.iloc[:,2].map(str) 
s1SVs = np.array(s1SVs)
s2SVs = s2SVs_df.iloc[:,0] + ':' + s2SVs_df.iloc[:,1].map(str) + '-' + s2SVs_df.iloc[:,2].map(str) 
s2SVs = np.array(s2SVs)
#s3SVs = s3SVs_df.iloc[:,0] + ':' + s3SVs_df.iloc[:,1].map(str) + '-' + s3SVs_df.iloc[:,2].map(str) 
#s3SVs = np.array(s3SVs)
#s6SVs = s6SVs_df.iloc[:,0] + ':' + s6SVs_df.iloc[:,1].map(str) + '-' + s6SVs_df.iloc[:,2].map(str) 
#s6SVs = np.array(s6SVs)
#s9SVs = s9SVs_df.iloc[:,0] + ':' + s9SVs_df.iloc[:,1].map(str) + '-' + s9SVs_df.iloc[:,2].map(str) 
#s9SVs = np.array(s9SVs)
#s10SVs = s10SVs_df.iloc[:,0] + ':' + s10SVs_df.iloc[:,1].map(str) + '-' + s10SVs_df.iloc[:,2].map(str) 
#s10SVs = np.array(s10SVs)

s0_VCF.drop(s0_VCF.columns[0:2], axis=1, inplace=True)
s0_VCF.insert(0, 'SV', s0SVs)
s1_VCF.drop(s1_VCF.columns[0:2], axis=1, inplace=True)
s1_VCF.insert(0, 'SV', s1SVs)
s2_VCF.drop(s2_VCF.columns[0:2], axis=1, inplace=True)
s2_VCF.insert(0, 'SV', s2SVs)
#s3_VCF.drop(s3_VCF.columns[0:2], axis=1, inplace=True)
#s3_VCF.insert(0, 'SV', s3SVs)
#s6_VCF.drop(s6_VCF.columns[0:2], axis=1, inplace=True)
#s6_VCF.insert(0, 'SV', s6SVs)
#s9_VCF.drop(s9_VCF.columns[0:2], axis=1, inplace=True)
#s9_VCF.insert(0, 'SV', s9SVs)
#s10_VCF.drop(s10_VCF.columns[0:2], axis=1, inplace=True)
#s10_VCF.insert(0, 'SV', s10SVs)

sampleSVs = [s0SVs, s1SVs, s2SVs]
sampleVCFs = [s0_VCF, s1_VCF, s2_VCF]

for s in range(0,3):

	os.system("bedtools intersect -wa -wb -a pd9770_tardis_validated.bed -b %s -r -f 0.70 > tmp.overlap.s%s" %(bedlist[s], samples[s]))
	print samples[s]
	tmp_df = pd.read_csv("tmp.overlap.s%s" %(samples[s]), sep="\t", header=None)
	os.system("rm tmp.overlap.s%s" %(samples[s]))

	sampledata = np.zeros(shape=(len(allSVs),2), dtype=object)
	#sampledata[:] = np.NAN

	seen_loc = []

	for i in range(0, len(tmp_df.iloc[:,0])):
		chr = tmp_df.iloc[i,0]
		start = tmp_df.iloc[i,1] + 1
		end = tmp_df.iloc[i,2]

		#sv as represented in allSamples file
		querysv = str(chr) + ':' + str(start) + '-' + str(end)

		chr_overlap = tmp_df.iloc[i,3]
		start_overlap = tmp_df.iloc[i,4] + 1
		end_overlap = tmp_df.iloc[i,5]

		path = ''

		if samples[s] == 'PD9770A':
			path = 'PD9770/a/PD9770A'
		elif samples[s] == 'PD9770C':
			path = 'PD9770/c/PD9770C'
		elif samples[s] == 'PD9770D':
			path = 'PD9770/d/PD9770D'

		print "running samtools depth"
		#calculate n_a and n_b read counts
		os.system("~/tools/samtools/bin/samtools depth -r %s:%s-%s %s > depth.tmp" %(chr_overlap, start_overlap, start_overlap, bamlist[s]))
		
		file = open("depth.tmp", "r")
		n_a = ''
		n_b = ''

		for line in file:
			n_a = line.strip().split('\t')[2]
		os.system("rm depth.tmp")

		#check if n_a is empty
		if n_a == "":
			n_a = 0

		os.system("~/tools/samtools/bin/samtools depth -a -r %s:%s-%s %s > depth.tmp" %(chr_overlap, end_overlap, end_overlap, bamlist[s]))
		
		file = open("depth.tmp", "r")
		
		for line in file:
			n_b = line.strip().split('\t')[2]
		os.system("rm depth.tmp")

		#check if n_b is empty
		if n_b == "":
			n_b = 0

		#Psuedo-extension of breakpoints
		start_overlap_ext = start_overlap - 151
		end_overlap_ext = end_overlap + 151

		#overlap SV as represented in the sample specific bed file
		querysv_overlap = str(chr_overlap) + ':' + str(start_overlap) + '-' + str(end_overlap)

		crp = ''
		print "running samtools view"
		os.system("~/tools/samtools/bin/samtools view -f 0x2 %s %s | awk '!seen[$1]++' | wc -l > counts.tmp" %(bamlist[s], querysv_overlap))

		file = open("counts.tmp", "r")
		
		#all_reads = 0
		for line in file:
			all_reads = int(line.strip())
			#all_reads + 1
		os.system("rm counts.tmp")

		#reset to correct endpoints
		querysv_overlap = str(chr_overlap) + ':' + str(start_overlap) + '-' + str(end_overlap)

		#pull information from the vcf files
		info = sampleVCFs[s].iloc[np.where(sampleSVs[s] == querysv_overlap)[0][0],1:3]

		res = re.search('(?<=SVLEN=)\w*', info.INFO)
		length_sv = res.group(0)
		info.INFO  = 'SIZE;N_A;N_B;SR;DRP;CRP;FRAG_LEN'

		#separate counts
		drp = info.Sample.strip().split(':')[4]
		splitreads = info.Sample.strip().split(':')[5]

		#if all_reads/2 > int(drp):
		crp = all_reads
		#else:
			#crp = 0

		all_counts = length_sv + ';' + str(n_a) + ';' + str(n_b) + ';' + str(splitreads) + ';' + str(drp) + ';' + str(crp) + ';' + frag_length[s]
		info.Sample = all_counts

		#find location of that SV in allSamplesSv file
		loc = np.where(allSVs == querysv)[0][0]

		if loc in seen_loc:	
			if querysv == querysv_overlap:

				#assign info to the sample data array to keep track of it
				sampledata[loc,:] = info
		else:
			seen_loc.append(loc)
		
			#assign info to the sample data array to keep track of it
			sampledata[loc,:] = info

	s_df = pd.DataFrame(data=sampledata, columns=['S%s_FORMAT' %(samples[s]), 'S%s_GenomeCounts' %(samples[s])])

	df1 = pd.concat([df1, s_df], axis=1)

df1.drop(df1.columns[0:3], axis=1, inplace=True)
df1.insert(0, 'SV', allSVs)
df1.to_csv(outputfile, sep='\t', index=False, header=True)

		

	




