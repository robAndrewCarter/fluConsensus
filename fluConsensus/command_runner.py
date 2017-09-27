import os, tempfile, sys, subprocess, re    
from fluConsensus import run_settings
from fluConsensus.utils import binner
import pandas as pd
import pysam
import numpy as np
from Bio import SeqIO, SeqRecord, Seq
import logging


#Required functions for python script:

def align_sequences(fastq_filename_list, target_reference_filename):
    '''
    This function takes a single or paired fastq filename as a list, a reference sequence
    filename, performs the alignment, and returns the name of the sorted,
    aligned bam. If the reference sequence does not have an indexed genome, one
    will be created and used throughout the function call.

    INPUT:

    fastq_filename_list: a list of length one or two that contains the
    filennames of the sample to align.

    reference_filename: a filename containing the reference sequences that the
    fastq filenames will be aligned against

    RETURNS:

    the filename of the resulting bam.
    '''

    #check that files exist
    for _filename in fastq_filename_list:
        if not os.path.exists(_filename):
            sys.exit("{} does not exist.".format(_filename))
    if not os.path.exists(target_reference_filename):
            sys.exit("{} does not exist.".format(target_reference_filename))

    #Check that the reference sequence is indexed.
    ensure_indexed_reference_file(target_reference_filename)

    #align sequences
    aln_filenames = []
    for _filename in fastq_filename_list:
        tmpfilename = tempfile.NamedTemporaryFile(dir=run_settings.global_args['temp_dir'], delete = False)
        tmpfilename.close()
        try:
            aln_command = [binner('bwa'), 'aln', '-f', tmpfilename.name,
                                   target_reference_filename, _filename]
            subprocess.check_call(aln_command)
            aln_filenames.append(tmpfilename.name)
        except subprocess.CalledProcessError as e:
            sys.exit("FAILED: " + " ".join(aln_command + ["\nEXCEPTION:".format(e.output)]))
    
    #Run SAMPE
    print run_settings.global_args['temp_dir']
    tmpsamfile = tempfile.NamedTemporaryFile(dir=run_settings.global_args['temp_dir'], delete = False)
    tmpsamfile.close()
    if len(aln_filenames) == 2:
        try:
            sampe_command = [binner('bwa'), 'sampe', '-f', tmpsamfile.name,
                                   target_reference_filename, aln_filenames[0],
                                   aln_filenames[1], fastq_filename_list[0],
                                   fastq_filename_list[1]]
            subprocess.check_call(sampe_command)
        except subprocess.CalledProcessError as e:
            sys.exit("FAILED: " + " ".join(sampe_command + "\nEXCEPTION: ".format(e.output)))
    elif len(aln_filenames) == 1:
        try:
            samse_command = [binner('bwa'), 'samse', '-f',
                                   tmpsamfile.name,target_reference_filename,
                                   aln_filenames[0],fastq_filename_list[0]]
            subprocess.check_call(samse_command)
        except subprocess.CalledProcessError as e:
            sys.exit("FAILED: " + " ".join(samse_command + "\nEXCEPTION: ".format(e.output)))

    #Convert to bam
    tmpbamfilename = tempfile.NamedTemporaryFile(dir=run_settings.global_args['temp_dir'], delete = False)
    tmpbamfilename.close()
    try:
        samtools_view_command = [binner('samtools'), 'view', "-o", tmpbamfilename.name,
                               '-Sb', tmpsamfile.name]
        subprocess.check_call(samtools_view_command)
    except subprocess.CalledProcessError as e:
        sys.exit("FAILED: " + " ".join(samtools_view_command + "\nEXCEPTION: ".format(e.output)))

    #Sort file
    tmpsortedbamfilename = tempfile.NamedTemporaryFile(dir=run_settings.global_args['temp_dir'], delete = False)
    tmpsortedbamfilename.close()
    try:
        samtools_sort_command =[binner('samtools'), 'sort', '-o',
                               tmpsortedbamfilename.name, tmpbamfilename.name] 
        subprocess.check_call(samtools_sort_command)
    except subprocess.CalledProcessError as e:
        sys.exit("FAILED: " + " ".join(samtools_sort_command + "\nEXCEPTION: ".format(e.output)))

#index file
    try:
        samtools_index_command = [binner('samtools'), 'index', tmpsortedbamfilename.name]
        subprocess.check_call(samtools_index_command)
    except subprocess.CalledProcessError as e:
        sys.exit("FAILED: " + " ".join(samtools_index_command + "\nEXCEPTION: ".format(e.output)))

    return tmpsortedbamfilename.name

def create_reference_file_from_accessions(accessions_list, sql_database):
    '''
    This function extracts the fasta files corresponding to the accessions_list
    from the sql_database, writes them to a randomly named file, then returns the filename

    INPUT:

    accessions_list: a list of accessions that map to sequences in the
    sql_database

    sql database: a database that contains a mapping of accessions to fasta
    sequences (as strings)

    RETURNS:

    a filename that contains the all reference sequences from accessions_list
    '''
    pass


def get_lowercase_seq_coverage(bam_filename, fasta_filename, min_coverage = None):
    ''' 
    This function reads a bam_filename bam and the corresponding fasta_filename fasta file.
    It converts the sequences from the FASTA file to a combination of uppercase and lowercase characters.
    Lowercase characters are placed wherever the coverage of the corresponding base is less than min_coverage.
    Uppercase characters are placed at all other character positions.

    INPUT:

    bam_filename: a bam filename that will be used to extract all variable
    positions

    fasta_filename: a fasta filename containing the reference sequences from
    bam_filename. This has to be included because bam files do not store the
    reference sequences.

    min_coverage: the minimum coverage of a position  in order for it to be considered covered.
    if None, then all positions are considered covered
    
    RETURNS:

    a dict mapping each sequence name to a seqrecord of the corresponding sequence.  
    The Seq object in the SeqRecord is converted to the upper/lowercase.
    '''

    if not os.path.exists(bam_filename):
        sys.exit("Bam file {} does not exist".format(bam_filename))
    if not os.path.exists(fasta_filename):
        sys.exit("Fasta file {} does not exist".format(fasta_filename))
    if min_coverage == None:
        min_coverage = 0
    try:
        bam_obj = pysam.Samfile(bam_filename, 'rb')
    except Exception:
        sys.exit("Could not load bam file {}".format(bam_filename))
    
    ref_name_to_length_dict = dict(zip(bam_obj.references, bam_obj.lengths))

    seg_to_seq_string_dict = {}
    try:
        for _seqrec in SeqIO.parse(fasta_filename, "fasta"):
            seg_to_seq_string_dict[_seqrec.id] = str(_seqrec.seq)
    except Exception:
        sys.exit("Could not parse " + fasta_filename)

    seqrec_list = []
    for _refname in seg_to_seq_string_dict.keys():
        tolower_inds_list = []
        counts = bam_obj.count_coverage(reference = _refname, start = 0, end = ref_name_to_length_dict[_refname])
        cov_array = np.asarray(counts)
        seq_base_array = list(seg_to_seq_string_dict[_refname])
        for _i in range(0,cov_array.shape[1]):
            if(sum(cov_array[:,_i]) < min_coverage):
                tolower_inds_list.append(_i)
        for _ind in tolower_inds_list:
            seq_base_array[_ind] = seq_base_array[_ind].lower()
        seqrec_list.append(SeqRecord.SeqRecord(seq = Seq.Seq("".join(seq_base_array)), id = _refname))
    return seqrec_list
    

def get_variants_from_bam(bam_filename, fasta_filename, threshold = None, min_coverage = None):
    '''
    This function reads a bam and gets the variants at each site that are
    greater than threshold in frequency. It only considers positions with
    minimum coverage of at least min_coverage. If min_coverage == None, the
    global run_settings.global_args['min_coverage'] is used

    INPUT:

    bam_filename: a bam filename that will be used to extract all variable
    positions

    fasta_filename: a fasta filename containing the reference sequences from
    bam_filename. This has to be included because bam files do not store the
    reference sequences.

    threshold: a float between 0 and 1.0, denoting the minimum frequency of a
    non-reference variant in order to be reported. If threshold == None, then the global settings are
    used, which is present in run_settings.global_args['threshold']


    min_coverage: the minimum coverage of a position in order to consider
    editing the sequence. If min_coverage == None, then the global settings are
    used, which is present in run_settings.global_args['min_coverage']

    RETURNS:

    a dataFrame with a row for each site in the bam_filename with a variant of
    higher frequency than threshold. The following columns are present:
        * reference - the reference name as it is contained in the bam
        * position - the position (in zero-based coordinates) of the variable
        position from the reference in the bam
        * ref_base - the reference base, as present in the bam
        * var_base - the variant base, with frququency above threshold
    '''
    #Make sure bam file exists
    if not os.path.exists(bam_filename):
        sys.exit("Bam file {} does not exist".format(bam_filename))
    if not os.path.exists(fasta_filename):
        sys.exit("Fasta file {} does not exist".format(fasta_filename))
    try:
        bam_obj = pysam.Samfile(bam_filename, 'rb')
    except Exception:
        sys.exit("Could not load bam file {}".format(bam_filename))

    ref_name_to_length_dict = dict(zip(bam_obj.references, bam_obj.lengths))
    #create mapping between base and integer
    #'N' is used for ALL non-ATCG bases
    base_to_ind_dict = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    ind_to_base_array = ['A','C','G','T','N']
    #load reference_sequences
    seg_to_seq_string_dict = {}
    for _seqrec in SeqIO.parse(fasta_filename, "fasta"):
        seg_to_seq_string_dict[_seqrec.id] = str(_seqrec.seq)

	#create empty array to fill with variable positions
	array_of_pos_to_change = []

    #Identify variable positions for editing
    for _ref_name, _ref_length in ref_name_to_length_dict.iteritems():
        base_freqs_array = np.array([list(_i) for _i in
                                  bam_obj.count_coverage(reference = _ref_name,
                                                         start = 0, end =
                                                         _ref_length)])
        ref_base_inds = []
        for _base in seg_to_seq_string_dict[_ref_name]:
          try:
            ref_base_inds.append(base_to_ind_dict[_base])
          except Exception:
            ref_base_inds.append(base_to_ind_dict['N'])
             
        base_counts_and_ref_ind_array = np.append(base_freqs_array,
                                                np.reshape(ref_base_inds,
                                                           [1,_ref_length]), axis = 0)
        for _pos in range(0, base_counts_and_ref_ind_array.shape[1]):
            temp_array = base_counts_and_ref_ind_array[0:4,_pos].transpose()
            if sum(temp_array) < min_coverage:
                continue
            max_args_inds = np.where(temp_array == np.max(temp_array))[0]
            if len(max_args_inds) > 1:
                continue
            max_args_inds = np.asscalar(max_args_inds)
            if (max_args_inds == base_counts_and_ref_ind_array[4,_pos]):
                continue
            else:
                print _pos
                base_freq = float(temp_array[max_args_inds])/sum(temp_array)
                if base_freq >= threshold:
                    array_of_pos_to_change.append([_ref_name, _pos, ind_to_base_array[base_counts_and_ref_ind_array[4,_pos]], ind_to_base_array[max_args_inds]])
    var_df = pd.DataFrame(array_of_pos_to_change, columns = ['reference', 'position', 'ref_base', 'var_base'])
    var_df.to_csv(os.path.join(run_settings.global_args['temp_dir'], re.sub("\.[^.]+$", "", bam_filename) + "_variants.tsv"), sep = "\t")
    return(var_df)

def get_bam_stats(bam_filename, fasta_filename, min_coverage = None):
    '''
    This function compiles a list of statistics for a bam file.
    Currently, it reports the fraction of each reference covered at or above min_coverage,
    the mean read coverage across each segment, the mean read coverage across the covered 
    portion of each reference (bases with at least min_coverage coverage)

    INPUT

    bam_filename: the name of the bam file to generate stats from

    fasta_filename: a fasta filename containing the reference sequences from
    bam_filename. This has to be included because bam files do not store the
    reference sequences.

    min_coverage: the minimum coverage required for a base to eb considered covered. 
    If None, then the global min_coverage is used

    RETURNS

    A dataframe containing the following rows:
        *ref, the reference name form the bam
        *pos, The nucleotide position of the ref
        *A, the A count at pos
        *C, the C count at pos
        *G, the G count at pos
        *T, the T count at pos
        *cov, the total coverage at pos

    '''

    #Make sure bam file exists
    if not os.path.exists(bam_filename):
        sys.exit("Bam file {} does not exist".format(bam_filename))
    if not os.path.exists(fasta_filename):
        sys.exit("Fasta file {} does not exist".format(fasta_filename))
    try:
        bam_obj = pysam.Samfile(bam_filename, 'rb')
    except Exception:
        sys.exit("Could not load bam file {}".format(bam_filename))

    ref_name_to_length_dict = dict(zip(bam_obj.references, bam_obj.lengths))
    
    seg_to_seq_string_dict = {}
    for _seqrec in SeqIO.parse(fasta_filename, "fasta"):
        seg_to_seq_string_dict[_seqrec.id] = str(_seqrec.seq)

    ref_cov_df_list = []
    ref_stats_list = []
    for _ref_name, _ref_length in ref_name_to_length_dict.iteritems():
        base_freqs_array = np.array([list(_i) for _i in
                                  bam_obj.count_coverage(reference = _ref_name,
                                                         start = 0, end =
                                                         _ref_length)])
        cov_df = pd.DataFrame(base_freqs_array.T, columns = ['A', 'C', 'G', 'T'])
        cov_df['cov'] = cov_df.apply(lambda x: sum(x.iloc[[0,1,2,3]]), axis = 1)
        cov_df['ref'] = _ref_name
        cov_df['pos'] = range(0, ref_name_to_length_dict[_ref_name])
        covered_mean_cov = cov_df['cov'][cov_df['cov'] >= min_coverage].mean()
        mean_cov = cov_df['cov'].mean()
        cov_pct = (cov_df['cov'] >= min_coverage).mean()
        ref_cov_df_list.append(cov_df)
        ref_stats_list.append([_ref_name, cov_pct, mean_cov, covered_mean_cov])
    cov_df = pd.concat(ref_cov_df_list, axis = 0)
    stats_df = pd.DataFrame(ref_stats_list, columns = ['ref', 'cov_pct', 'mean_cov', 'covered_mean_cov'])
    return(cov_df, stats_df) 

def edit_reference_sequences(reference_fasta_filename, edits_df):
    '''
    This function reads reference_fasta_filename, edits the corresponding
    positions from the edits_df dataframe, saves the resulting fasta file,
    then returns its filename.

    INPUT:

    reference_fasta_filename: a filename of a reference fasta file
    fcontaining the references from the bam file from which  the edit_df was
    derived.

    edits_df: a dataFrame with a row for each site to be edited in the sequences
    contained in the reference_fasta_filename. The following columns are present:
        * reference - the reference name as it is contained in the bam
        * position - the position (in zero-based coordinates) of the variable
        position from the reference in the bam
        * ref_base - the reference base, as present in the bam
        * var_base - the variant base, with frququency above threshold

    RETURN:

    the fasta filename of the edited sequences
    
    '''
    if not os.path.exists(reference_fasta_filename):
        sys.exit("File {} does not exist".format(reference_fasta_filename))
    if os.path.getsize(reference_fasta_filename) == 0:
        sys.exit("File {} is empty".format(reference_fasta_filename))

    #Load reference sequences and convert them into ndarrays of letters
    seg_to_seq_ndarray_dict = {}

    for _seqrec in SeqIO.parse(reference_fasta_filename, "fasta"):
        seg_to_seq_ndarray_dict[_seqrec.id] = np.array(list(str(_seqrec.seq)))
        print 'foo'
    print seg_to_seq_ndarray_dict
    for _ref_name, group in edits_df.groupby(['reference']):
        #assertion removed because there will be mismatches between N ands some other non-ACGT bases
        #assert(np.array_equal(np.array(list(seg_to_seq_ndarray_dict[_ref_name][group.loc[:,'position']])),
        #               np.array(list(group.loc[:,"ref_base"].values))))
        seg_to_seq_ndarray_dict[_ref_name][np.array(group.loc[:,'position'])] = np.array(group.loc[:,'var_base'])
    edited_fasta_tempfile = tempfile.NamedTemporaryFile(dir =
                                                        run_settings.global_args['temp_dir'],
                                                        delete=False)
    seqrec_list = []
    for _ref_name, _seq_ndarray in seg_to_seq_ndarray_dict.iteritems():
        new_seqrec = SeqRecord.SeqRecord(id = _ref_name, seq=
                                             Seq.Seq(''.join(seg_to_seq_ndarray_dict[_ref_name])),
                                             description = '')
        seqrec_list.append(new_seqrec)
        print(new_seqrec)
    try:
        SeqIO.write(seqrec_list, handle = edited_fasta_tempfile, format = 'fasta')
    except Exception:
        sys.exit("Could not write sequences to file {}".format(edited_fasta_tempfile.name))
    edited_fasta_tempfile.close
    return(edited_fasta_tempfile.name)
        

def ensure_indexed_reference_file(fasta_filename):
    '''
    This functions ensures that fasta_filename is indexed by bwa. Specifically,
    it looks for files with the same name as the fasta_filename, but with the
    following suffixes: [bwt, pac, ann, amb, sa]

    INPUT:
        fasta_filename: path to reference fasta sequence
    
    RETURN:
        True if fasta_filename is indexed by bwa by the end of the function
        call
    '''
    if all([os.path.exists(".".join([fasta_filename, _bwa_index_suffix])) for _bwa_index_suffix in
            ['bwt','pac','ann','amb','sa']]):
        return True
    else:
        try:
            index_command = [binner('bwa'), 'index', fasta_filename]
            subprocess.check_call(index_command)
        except Exception as e:
            sys.exit('Cannot create index for file {}'.format(fasta_filename))

