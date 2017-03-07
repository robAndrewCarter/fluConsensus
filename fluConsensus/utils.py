from fluConsensus import run_settings
import os

def binner(command):
    '''
    This command merely adds the bin_dir folder path to the command and
    returns the full path. This makes the program more self-contained since it
    does not rely on any of the programs, such as BWA or Samtools, being on the
    search path.
    INPUT:
        command: a string of the command to be executed
    RETURN:
        the full path to the command within the external directory
    '''
    return os.path.join(run_settings.global_args['bin_dir'], command)
