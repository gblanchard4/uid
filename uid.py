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
	# IUPAC lookup table, see http://www.bioinformatics.org/sms/iupac.html
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

# Check that file are a pair by comparing sequence headers
def matching_files(forward_id, reverse_id):
	return forward_id == reverse_id

# Check if their are any N bases in a sequence
def any_N_bases(sequence):
	return 'N' in sequence

# Regular expression match primer
def match_primer(regular_expression,sequence):
	re_search_result = re.findall(regular_expression, str(sequence))
	return (re_search_result)

# Get uid from the sequence, split on the first finding
def get_uid(primer,sequence):
	return str(sequence.split(primer,1)[0])


def main():

	'''
	Get command line arguments
	'''
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

	# Quality filter
	parser.add_argument('-q','--qual',dest='qual', help='The quality score to filter. Keep reads with bases <= Q', required=True)

	args = parser.parse_args()

	# Input files
	forward_file = args.forward
	reverse_file = args.reverse
	# Primer files
	forward_primer = args.fprimer
	reverse_primer = args.rprimer
	# Qual
	qual = int(args.qual)


	# Compile the regular expression for the primers
	forward_re = re_compiler(forward_primer)
	reverse_re = re_compiler(reverse_primer)

	# Counters
	ndians = 0 # How many reads with N bases
	forward_ndians = 0
	reverse_ndians = 0
	forward_primer_mismatches = 0 # How many reads with forward primer mismatches
	reverse_primer_mismatches = 0 # How many reads with reverse primer mismatches
	multi_forward_primer_match = 0
	multi_reverse_primer_match = 0
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
	uid_passing_filter = 0
	forward_uid_quality_total = 0
	reverse_min_quality = {}
	quality_distribution = {}
	reverse_uid_quality_total = 0

	#Open our files for IO
	with open(forward_file, 'rU') as forward_file_handle, open(reverse_file, 'rU') as reverse_file_handle, open(forward_file+'.UID', 'w+') as output_r1_handle, open(reverse_file+'.UID', 'w+') as output_r2_handle:
		# Zip both files together as a SeqIO object so we can easily process both files at once
		for forward_record, reverse_record in izip(SeqIO.parse(forward_file_handle,"fastq"), SeqIO.parse(reverse_file_handle,"fastq")):
			# Increment read counter
			total_reads = total_reads + 1

			# Make sure the files are a matching pair
			if not matching_files(forward_record.id, reverse_record.id):
				print "ID mismatch error! Make sure you put in a correct pair of files"
				sys.exit()

			# Check for a N bases, increment counter to keep track of stats
			if any_N_bases(forward_record.seq) or any_N_bases(reverse_record.seq):
				if any_N_bases(forward_record.seq):
					forward_ndians += 1
				if any_N_bases(reverse_record.seq):
					reverse_ndians += 1
				continue

			# Regular Expression search for a primer match
			#>> Forward
			forward_primer_list = match_primer(forward_re, forward_record.seq)
			#>> Reverse
			reverse_primer_list = match_primer(reverse_re, reverse_record.seq)

			# Handle more than one match
			if  len(forward_primer_list) > 1 or len(reverse_primer_list) > 1:
				if len(forward_primer_list) > 1 :
					print "MULTIPLE PRIMER MATCH WARINING IN FORWARD FILE FOR:\n\t%s" % forward_record.id
					#multiple_forward_primer_match += 1
				if len(reverse_primer_list) > 1:
					print "MULTIPLE PRIMER MATCH WARINING IN REVERSE FILE FOR:\n\t%s" % reverse_record.id
					#multi_reverse_primer_match += 1
				continue
			# Handle No Match
			if len(forward_primer_list) == 0 or len(reverse_primer_list)== 0:
				if len(forward_primer_list) == 0 and len(reverse_primer_list)== 0:
					complete_primer_mismatches += 1
				if len(forward_primer_list) == 0:
					forward_primer_mismatches += 1
				if len(reverse_primer_list) == 0:
					reverse_primer_mismatches += 1
				continue


			# Make sure UID is of the right length
			forward_uid = get_uid(forward_primer_list[0], forward_record.seq)
			reverse_uid = get_uid(reverse_primer_list[0], reverse_record.seq)
			uid_tuple = (forward_uid,reverse_uid)

			# Length not 10 or 12
			if not ((len(forward_uid) == 0 and (len(reverse_uid) == 10 or len(reverse_uid) ==12)) or (len(reverse_uid) == 0 and (len(forward_uid) == 10 or len(forward_uid) ==12))):
				# UID mismatch total
				uid_length_mismatches = uid_length_mismatches + 1
				# Forward UID length distribution
				try:
					forward_uid_length[len(forward_uid)] = forward_uid_length[len(forward_uid)] + 1
				except KeyError:
					forward_uid_length[len(forward_uid)] = 1
				# Reverse UID length distribution
				try:
					reverse_uid_length[len(reverse_uid)] = reverse_uid_length[len(reverse_uid)] + 1
				except KeyError:
					reverse_uid_length[len(reverse_uid)] = 1
				# Total UID length distribution
				total_UID = len(forward_uid) + len(reverse_uid)
				try:
					total_uid_length[total_UID] = total_uid_length[total_UID] + 1
				except KeyError:
					total_uid_length[total_UID] = 1
				continue
			else:
				uid_length_matches += 1
				if (len(forward_uid + reverse_uid)) == 10:
					uid_length_10 += 1
				if (len(forward_uid + reverse_uid)) == 12:
					uid_length_12 += 1

			# Get UID quality
			forward_uid_quality = forward_record.letter_annotations["phred_quality"][:len(forward_uid)]
			reverse_uid_quality = reverse_record.letter_annotations["phred_quality"][:len(reverse_uid)]

			try:
				quality_distribution[min(reverse_uid_quality)] = quality_distribution[min(reverse_uid_quality)] +1
			except KeyError:
				quality_distribution[min(reverse_uid_quality)] = 1
			"""
			Need to generalize
			"""
			if min(reverse_uid_quality) >= qual:
				# Write our output, QA/QC done
				uid = forward_uid + reverse_uid
				# Forward file
				fheader,fseq,fplus,fqual =  forward_record.format("fastq").rstrip().split('\n')
				fheader = fheader +' UID='+uid
				fseq = fseq[len(forward_uid):]
				fqual = fqual[len(forward_uid):]
				fstring = '%s\n%s\n%s\n%s\n' % (fheader, fseq, fplus, fqual)
				output_r1_handle.write(fstring)
				# Reverse file
				rheader,rseq,rplus,rqual =  reverse_record.format("fastq").rstrip().split('\n')
				rheader = rheader +' UID='+uid
				rseq = rseq[len(reverse_uid):]
				rqual = rqual[len(reverse_uid):]
				rstring = '%s\n%s\n%s\n%s\n' % (rheader, rseq, rplus, rqual)
				output_r2_handle.write(rstring)
				uid_passing_filter += 1
			# UID is low quality, toss it
			else:
				if min(reverse_uid_quality) < qual:
					reverse_uid_quality_total = reverse_uid_quality_total + 1
					try:
						reverse_min_quality[min(reverse_uid_quality)] = reverse_min_quality[min(reverse_uid_quality)] +1
					except KeyError:
						reverse_min_quality[min(reverse_uid_quality)] = 1

	# Stats
	print "Files: %s\t%s" % (forward_file, reverse_file)

if __name__ == '__main__':
	main()
