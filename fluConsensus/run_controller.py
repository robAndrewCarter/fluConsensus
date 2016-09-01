import run_settings
import command_runner
import log_writer
import sys, subprocess
run_settings.init(args.work_dir, args.fastq_filenames, args.min_correction_coverage )
#Basic steps that need to be performed


#Allow first round of alignment
another_iteration = True;
target_ref_filename = args.reference_filename
while(another_iteration):
    aligned_bam_filename = command_runner.align_sequences(run_settings.fastq_list,
                                   target_ref_filename)
    variants_df = command_runner.get_variants_from_bam(aligned_bam_filename,
                                         run_settings.threshold)
    if not variants_df.empty:
        log_writer.update_variants(variants_df)
        target_ref_filename = command_runner.edit_reference_sequences(target_ref_filename, variants_df)
    else:
        another_iteration = False
if subprocess.check_call(["mv", target_ref_filename,
                              args.output_filename]):
    sys.exit("Problem moving file")
    return(0)
else:
    log_writer.log_result(args.output_filename)
