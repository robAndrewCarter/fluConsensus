import os, sys

def init(sample_name, work_dir, reference_filename, fastq_filenames_list, threshold, min_correction_coverage, max_iterations, min_coverage):
    global global_args
    global_args = {}
    global_args['sample_name'] = sample_name
    global_args['work_dir'] = work_dir
    global_args['reference_filename'] = reference_filename
    global_args['fastq_filenames_list'] = fastq_filenames_list
    global_args['threshold'] = threshold
    global_args['min_correction_coverage'] = min_correction_coverage
    global_args['max_iterations'] = max_iterations
    global_args['min_coverage'] = min_coverage

    global_args['bin_dir'] = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(sys.modules[__name__].__file__))),
        "external")

    #Define global variables that are composed if input variables
    global_args['temp_dir'] = os.path.join(work_dir, sample_name)
