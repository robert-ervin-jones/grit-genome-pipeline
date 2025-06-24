import subprocess
import argparse
import os

def make_config_1():
    '''
    Make config list to be passed to

    :return: snakemake config
    :rtype: list
    '''
    ret = [
        f'qscores={args.q_filter}',
        f'lengths={args.lengths}',
        f'wrapper_dir={os.path.dirname(os.path.abspath(__file__))}',    
    ]

    return ' '.join(ret)


def run_snakemake(length_filter=False):
    '''
    Combine snakemake cmd frags and runs snakemake
    '''
    smk_frags = [
        f'snakemake',
        f'--config {make_config_1()}',
        f'--cores {args.threads}',
        f'--rerun-incomplete',
        f'--use-conda',
        f'--use-apptainer',
        f'--nolock'
    ]

    if args.dry_run:
        smk_frags.append('-n')

    smk_cmd = ' '.join(smk_frags)
    subprocess.run(smk_cmd, shell=True, executable='/bin/bash', check=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='Run metagenome snakemake pipeline',
                    epilog='')
    # Optional Arguments
    parser.add_argument('-t', '--threads', metavar='N', type=int, default=1, help='Desired number of threads to be utilized. Default: 1')
    parser.add_argument('--q_filter', metavar='N', type=str, default=20, help='Q score filter value. Default: 20')
    parser.add_argument('--lengths', metavar='N', type=str, default=5000, help='Length filter value. Default: 5000')
    # TODO: Add option set flye to use nano-raw
    parser.add_argument('-n', '--dry_run', action='store_true', default=False)

    # Changes help descriptions from the default input and output help descriptions
    args = parser.parse_args()
    
    
    run_snakemake()
    run_mfannot()
