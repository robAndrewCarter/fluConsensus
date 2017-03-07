import pandas as pd
from fluConsensus import command_runner
from fluConsensus import run_settings
import pytest, os, sys, tempfile
from Bio import SeqIO, SeqRecord, Seq



ref_filename = "test_na_ref.fa"
fwd_fastq_filename = "test_na_variant_seq_1.fastq"
rev_fastq_filename = "test_na_variant_seq_2.fastq"
module_dir = os.path.dirname(os.path.abspath(sys.modules[__name__].__file__))

#Create global variables for running all the functions
run_settings.init(sample_name = 'test', work_dir = module_dir,
                  reference_filename = os.path.join(module_dir, ref_filename), fastq_filenames_list =
                  [os.path.join(module_dir, fwd_fastq_filename),
                   os.path.join(module_dir, rev_fastq_filename)], threshold =
                  0.5, min_correction_coverage = 2 )
print run_settings.global_args
print "Here are the global settings:"

def test_alignment_and_variant_identification():
    run_settings.init(sample_name = 'test', work_dir = module_dir,
                  reference_filename = os.path.join(module_dir, ref_filename), fastq_filenames_list =
                  [os.path.join(module_dir, fwd_fastq_filename),
                   os.path.join(module_dir, rev_fastq_filename)], threshold =
                  0.5, min_correction_coverage = 2 )
    print run_settings.global_args
    test_bam_filename = command_runner.align_sequences(run_settings.global_args['fastq_filenames_list'],
                                   run_settings.global_args['reference_filename'])
    print("Testing command_runner.align_sequences function to determine if bam file is generated")
    assert(isinstance(test_bam_filename, basestring))
    print("Testing command_runner.get_variants_from_bam function to determine if the variants can be identified and returned as a pandas DataFrame")
    print(test_bam_filename)
    var_df = command_runner.get_variants_from_bam(test_bam_filename,
                                         run_settings.global_args['reference_filename'])
    target_var_df = pd.DataFrame({'reference': ['test_na_variant_seq', 'test_na_variant_seq',
                                'test_na_variant_seq'], 'position': [280, 560,
                                                                     840], 'ref_base':['G', 'T', 'A'], 'var_base':['T', 'C', 'T']})
    var_df = var_df.sort_index(axis = 1)
    target_var_df = target_var_df.sort_index(axis = 1)
    assert(var_df.equals(target_var_df))
    print """Successfully aligned sequences and detected the variants"""

def test_edit_sequence():
    seq1 = SeqRecord.SeqRecord(id = 'temp_seq1', seq = Seq.Seq('ATCGATCG'))
    seq2 = SeqRecord.SeqRecord(id = 'temp_seq2', seq = Seq.Seq('GCTAGCTA'))
    tmp_file = tempfile.NamedTemporaryFile()
    SeqIO.write([seq1, seq2], tmp_file.name, 'fasta')
    assert os.path.getsize(tmp_file.name) > 0
    edit_df = pd.DataFrame({'reference' : ['temp_seq1','temp_seq1', 'temp_seq2'],
                            'position' : [2, 4, 6], 'ref_base' : ['C', 'A',
                                                                 'T'],
                            'var_base' : ['G', 'G', 'G']})
    edited_filename  = command_runner.edit_reference_sequences(tmp_file.name,  edit_df)
    assert(isinstance(edited_filename, basestring))
    seqrec_id_to_seqrec_dict = {sr.id : str(sr.seq) for sr in SeqIO.parse(edited_filename, 'fasta')}
    print seqrec_id_to_seqrec_dict
    seq1_str = seqrec_id_to_seqrec_dict['temp_seq1']
    seq2_str = seqrec_id_to_seqrec_dict['temp_seq2']
    assert seq1_str == 'ATGGGTCG'
    assert seq2_str == 'GCTAGCGA'

def run_full_pipeline_and_compare_sequences():
   ref_seqrec = SeqIO.read(ref_filename)
