#!/usr/bin/env python
from fluConsensus import run_settings, command_runner, log_writer
from Bio import SeqIO
import sys, subprocess, argparse, os, logging

def parse_arguments():
    parser= argparse.ArgumentParser(description="""This script repeatedly aligns either a
                                    single or paired fastq file(s) to a
                                    reference sequence(s), correcting all
                                    variable positions in the reference at each
                                    iteration.
                                    """)
    parser.add_argument('--sample_name', help="""Name of sample. This will be used
                        for file and temporary directory naming purposes.
                        defaults to tmp/ under work_dir""")
    parser.add_argument('--work_dir', default = os.curdir, help="""The
                        working directory where the final edited reference
                        sequences will be saved. Defaults to the current
                        directory""")
    parser.add_argument('--threshold', type = float, default = 0.5, help =
                        """Minimum variant frequency threshold for editing the
                        reference sequence""")
    parser.add_argument('--min_correction_coverage', type = float, default =
                        100, help = """Minimum coverage of a variable position
                        in order to induce a reference sequence change""")
    parser.add_argument('--min_coverage', type = float, default =
                        10, help = """Minimum coverage of position in order
                        to be considered covered. Covered positions are uppercase,
                        wherease uncovered positions are lowercase in the final
                        consensus sequence.""")
    parser.add_argument('--max_iterations', type = int, default = 10, help =
                        """Maximum number of alignment-consensus iterations to perform""")
    parser.add_argument('reference_filename', type = str, default = """""")
    parser.add_argument('fastq_filenames', nargs = '+', type = str, default =
                        """Either one or two fastq filenames that will be
                        aligned against the reference_filename and used for
                        interative correction.""")
    args = parser.parse_args()
    run_settings.init(args.sample_name, args.work_dir, args.reference_filename,
                      args.fastq_filenames, args.threshold, args.min_correction_coverage, args.max_iterations, args.min_coverage)

parse_arguments()

#create logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

#create temp directory and working directory if they don't already exist
if not os.path.exists(run_settings.global_args['work_dir']):
    os.mkdir(run_settings.global_args['work_dir'])
if not os.path.exists(run_settings.global_args['temp_dir']):
    os.mkdir(run_settings.global_args['temp_dir'])

#Allow first round of alignment
iteration_counter = 0
another_iteration = True
aligned_bam_filename = None
target_ref_filename = run_settings.global_args['reference_filename']

#Iteratively align and correct until no further improvement
while(another_iteration and (iteration_counter < run_settings.global_args['max_iterations'])):
    logging.info("Performing a round of alignment and variant calling")
    aligned_bam_filename = command_runner.align_sequences(run_settings.global_args['fastq_filenames_list'],
                                   target_ref_filename)
    variants_df = command_runner.get_variants_from_bam(aligned_bam_filename,
                                                       target_ref_filename,
                                         run_settings.global_args['threshold'])
    iteration_counter += 1
    if not variants_df.empty:
        #log_writer.update_variants(variants_df)
        logging.info("Identified {} variants".format(variants_df.shape[0]))
        target_ref_filename = command_runner.edit_reference_sequences(target_ref_filename, variants_df)
    else:
        logging.info("No More variants identified")
        another_iteration = False

[cov_df, stats_df] = command_runner.get_bam_stats(aligned_bam_filename, target_ref_filename, run_settings.global_args['min_coverage'])
cov_df.to_csv(os.path.join(run_settings.global_args['temp_dir'], run_settings.global_args['sample_name'] + "_cov.tsv"), sep = '\t')
stats_df.to_csv(os.path.join(run_settings.global_args['temp_dir'], run_settings.global_args['sample_name'] + "_stats.tsv"), sep = '\t')

#convert low-coverage regions to lowercase and save to fasta file
try:
    print "Just a BAR!"
    consensus_seqrecs = command_runner.get_lowercase_seq_coverage(aligned_bam_filename,
                                                           target_ref_filename,
                                             run_settings.global_args['min_coverage'])
    print consensus_seqrecs
    print "FOO and a BAR!"
    SeqIO.write(consensus_seqrecs, os.path.join(run_settings.global_args['work_dir'],
                                               "".join([run_settings.global_args['sample_name'],
                                               '_consensus_seqs.fa'])),"fasta" )

except Exception as e:
    sys.exit("Problem converting consensus sequences to lower case or problem saving file")
