import run_settings
#Required functions for python script:

def align_sequences(fastq_filename_list, reference_filename):
    '''
    This function takes a single or paired fastq filename, a reference sequence
    filename, performs the alignment, and returns the name of the sorted, aligned bam
    
    INPUT:

    fastq_filename_list: a list of length one or two that contains the
    filennames of the sample to align.

    reference_filename: a filename containing the reference sequences that the
    fastq filenames will be aligned against

    RETURNS:
        
    the filename of the resulting bam.
    '''
    pass

def create_reference_file_from_accessions(accessions_list, sql_database):
    '''
    This function extracts the fasta files corresponding to the accessions_list
    from the sql_database, writes them to a file, then returns the filename

    INPUT:

    accessions_list: a list of accessions that map to sequences in the
    sql_database

    sql database: a database that contains a mapping of accessions to fasta
    sequences (as strings)

    RETURNS:

    a filename that contains the all reference sequences from accessions_list
    '''
    pass

def get_variants_from_bam(bam_filename, threshold = 0.5):
    '''
    This function reads a bam and gets the variants at each site that are
    greater than threshold in frequency

    INPUT:

    bam_filename: a bam filename that will be used to extract all variable
    positions

    threshold: a float between 0 and 1.0, denoting the minimum frequency of a
    non-reference variant in order to be reported

    RETURNS:

    a dataFrame with a row for each site in the bam_filename with a variant of
    higher frequency than threshold. The following columns are present:
        * reference - the reference name as it is contained in the bam
        * position - the position (in zero-based coordinates) of the variable
        position from the reference in the bam
        * ref_base - the reference base, as present in the bam
        * var_base - the variant base, with frququency above threshold
        * ref_base - the reference base 
    '''
    pass

def edit_reference_sequences(reference_fasta_filename, edits_df):
    '''
    This function reads reference_fasta_filename, edits the corresponding
    positions from the edits_df dataframe, saves the resulting fasta file,
    thenn returns its filename.

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
        * ref_base - the reference base 

    RETURN:

    the fasta filename of the edited sequences
    
    '''
