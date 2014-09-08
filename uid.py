#! /usr/bin/env python

import argparse
import re
import sys
from itertools import izip
from Bio import SeqIO

__author__ = "Gene Blanchard"
__email__ = "me@geneblanchard.com"

'''
DOC
'''
def re_compiler(primer):
	# IUPAC lookup table
	IUPAC_dict = {'R':"[A,G]"\
				,'Y':"[C,T]"\
				,'S':"[G,C]"\
				,'W':"[A,T]"\
				,'K':"[G,T]"\
				,'M':"[A,C]"\
				,'B':"[C,G,T]"\
				,'D':"[A,G,T]"\
				,'H':"[A,C,T]"\
				,'V':"[A,C,G]"\
				,'N':"[A,C,G,T]"\
				}
	for key in IUPAC_dict.iterkeys():
		if key in primer:
			primer = primer.replace(key, IUPAC_dict[key])
	return primer

def main():
	parser = argparse.ArgumentParser(description='Calculate stats on the input fastqs')

	#Files
	# Forward file
	parser.add_argument('-f','--forward',dest='forward', help='The forward fastq', required=True)
	# Reverse file
	parser.add_argument('-r','--reverse',dest='reverse', help='The reverse fastq', required=True)

	#Primer
	# Forward primer
	parser.add_argument('-fp','--fprimer',dest='fprimer', help='The forward primer', required=True)
	# Reverse primer
	parser.add_argument('-rp','--rprimer',dest='rprimer', help='The reverse primer', required=True)

	args = parser.parse_args()

	forward_file = args.forward
	reverse_file = args.reverse
	forwad_primer = args.fprimer
	reverse_primer = args.rprimer

	# Hardcoded primer values (See http://www.bioinformatics.org/sms/iupac.html)
	v1_primer = "AGAGTTTGATYMTGGCTCAG"
	v3_primer = "ATTACCGCGGCTGCTGGC"

	v1_re = re_compiler(v1_primer)
	v3_re = re_compiler(v3_primer)

	# Open both the files for reading at the same time
	forward_file_handle = open(forward_file, 'rU')
	reverse_file_handle = open(reverse_file, 'rU')
	output_r1_handle = open(forward_file+'.UID', 'w+')
	output_r2_handle = open(reverse_file+'.UID', 'w+')


	# Counters
	ndians = 0 # How many reads with N bases
	forward_primer_mismatches = 0 # How many reads with forward primer mismatches
	reverse_primer_mismatches = 0 # How many reads with reverse primer mismatches
	complete_primer_mismatches = 0 # How many reads with complete primer mismatches
	uid_q_score_low = 0 # How many reads had too low of q scores
	uid_length_mismatches = 0
	uid_length_matches = 0
	total_reads = 0
	mismatch_offset = {}
	forward_uid_length = {}
	reverse_uid_length = {}
	total_uid_length = {}
	uid_length_10 = 0
	uid_length_12 = 0
	forward_UID_quality_total = 0
	reverse_UID_quality_total = 0
	forward_min_quality = {}
	reverse_min_quality = {}
	uid_passing_filter = 0

	# Zip both files together as a SeqIO object so we can easily process both files at once
	for forward_record, reverse_record in izip(SeqIO.parse(forward_file_handle,"fastq"), SeqIO.parse(reverse_file_handle,"fastq")):
		# Make sure the records match, otherwise fail and check files
		total_reads = total_reads + 1
		if forward_record.id == reverse_record.id:
			# Check for N bases
			if not 'N' in forward_record.seq and not 'N' in reverse_record.seq:
				# Check for primer match using regular expressions (see magic)
				# Seach forward record
				forward_record_primer = re.findall(v1_re, str(forward_record.seq))
				reverse_record_primer = re.findall(v3_re, str(reverse_record.seq))
				# If there is a match continue
				if len(forward_record_primer) > 0 and len(reverse_record_primer) > 0:
					# Identify UID as the substring before the primer
					forward_UID = str(forward_record.seq).split(forward_record_primer[0])[0]
					reverse_UID = str(reverse_record.seq).split(reverse_record_primer[0])[0]
					#if (len(forward_UID) + len(reverse_UID)) == 10 or  (len(forward_UID) + len(reverse_UID)) == 12:
					if (len(forward_UID) == 0 and (len(reverse_UID) == 10 or len(reverse_UID) ==12)) or (len(reverse_UID) == 0 and (len(forward_UID) == 10 or len(forward_UID) ==12)):
						uid_length_matches = uid_length_matches + 1
						if len(forward_UID) == 10 or len(reverse_UID) == 10:
							uid_length_10 += 1
						if len(forward_UID) == 12 or len(reverse_UID) == 12:
							uid_length_12 += 1

						# Get UID quality
						if len(forward_UID) == 0:
							forward_UID_quality = forward_record.letter_annotations["phred_quality"][:len(forward_UID)]
						else:
							forward_UID_quality = forward_record.letter_annotations["phred_quality"][:(len(forward_UID)-1)]

						if len(reverse_UID) == 0:

							reverse_UID_quality = reverse_record.letter_annotations["phred_quality"][:len(reverse_UID)]
						else:
							reverse_UID_quality = reverse_record.letter_annotations["phred_quality"][:(len(reverse_UID)-1)]
						# If the Q-score in the UID is 2 or less toss it (Probably not needed)

						if not False in [q>=11 for q in reverse_UID_quality]:
							# Write our output, QA/QC done
							# Forward file
							# Add UID to end of description
							forward_record.description = forward_record.description +' '+forward_UID
							#output_r1_handle.write(forward_record.format('fastq'))
							# Reverse file
							# Add UID to end of description
							reverse_record.description = reverse_record.description +' '+reverse_UID
							#output_r2_handle.write(reverse_record.format('fastq'))
							uid_passing_filter += 1
						# UID is low quality, toss it
						else:
							if False in [q<11 for q in forward_UID_quality] and not forward_UID_quality:
								forward_UID_quality_total = forward_UID_quality_total + 1
								try:
									forward_min_quality[min(forward_UID_quality)] = forward_min_quality[min(forward_UID_quality)] +1
								except KeyError:
									forward_min_quality[min(forward_UID_quality)] = 1
							if False in [q<11 for q in reverse_UID_quality] and not reverse_UID_quality:
								reverse_UID_quality_total = reverse_UID_quality_total + 1
								try:
									reverse_min_quality[min(reverse_UID_quality)] = reverse_min_quality[min(reverse_UID_quality)] +1
								except KeyError:
									reverse_min_quality[min(reverse_UID_quality)] = 1
							pass
					# UID not 10 or 12
					else:
						# UID mismatch total
						uid_length_mismatches = uid_length_mismatches + 1
						# Forward UID length distribution
						try:
							forward_uid_length[len(forward_UID)] = forward_uid_length[len(forward_UID)] + 1
						except KeyError:
							forward_uid_length[len(forward_UID)] = 1
						# Reverse UID length distribution
						try:
							reverse_uid_length[len(reverse_UID)] = reverse_uid_length[len(reverse_UID)] + 1
						except KeyError:
							reverse_uid_length[len(reverse_UID)] = 1
						# Total UID length distribution
						total_UID = len(forward_UID) + len(reverse_UID)
						try:
							total_uid_length[total_UID] = total_uid_length[total_UID] + 1
						except KeyError:
							total_uid_length[total_UID] = 1
						pass
				# Primer mismatch
				else:
					if len(forward_record_primer) == 0 and len(reverse_record_primer) > 0:
						forward_primer_mismatches = forward_primer_mismatches + 1
					if len(forward_record_primer) > 0 and len(reverse_record_primer) == 0:
						reverse_primer_mismatches = reverse_primer_mismatches + 1
					if len(forward_record_primer) == 0 and len(reverse_record_primer) == 0:
						complete_primer_mismatches = complete_primer_mismatches + 1
					pass
			# N base found move on
			else:
				ndians = ndians +1
				pass

		else:
			print "ID mismatch error! Make sure you put in a correct pair of files"
			forward_file_handle.close()
			reverse_file_handle.close()
			output_r1_handle.close()
			output_r2_handle.close()
			sys.exit()


	# Stats
	print "Total reads: %s" % total_reads
	print "Reads with N bases: %s" % ndians
	print "\n***Primer Filtering***"
	print "Forward primer mismatches: %s" % forward_primer_mismatches
	print "Reverse primer mismatches: %s" % reverse_primer_mismatches
	print "Complete primer mismatches: %s" % complete_primer_mismatches
	print "Total primer mismatches: %s" % (forward_primer_mismatches+reverse_primer_mismatches+complete_primer_mismatches)
	print "\n***UID Matching***"
	print "UID length matches: %s" % uid_length_matches
	print "UID length mismatches: %s " % uid_length_mismatches
	print "Forward UID length distribution: %s" % forward_uid_length
	print "Reverse UID length distribution: %s" % reverse_uid_length
	print "Total UID length distribution: %s" % total_uid_length
	print "UID's of length 10: %s" % uid_length_10

	print "UID's of length 12: %s" % uid_length_12
	print "\n***UID Quality Filtering**"
	print "Forward UID bad quality: %s" % forward_UID_quality_total
	print "Reverse UID bad quality: %s" % reverse_UID_quality_total
	print "Forward_UID_quality distribution %s" % forward_min_quality
	print "Reverse_UID_quality distribution %s" % reverse_min_quality
	print "UID Passing filter: %s" % uid_passing_filter

	# Close down the files
	forward_file_handle.close()
	reverse_file_handle.close()
	output_r1_handle.close()
	output_r2_handle.close()

if __name__ == '__main__':
	main()
