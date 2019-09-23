import os
import subprocess as sp
from typing import Iterable

import pandas as pd


def read_kraken_report(path):
    fields = 'perc n_reads n_reads_specific rank taxid sciname'.split()
    return pd.read_csv(path, sep='\t', header=None, names=fields)


def read_kraken2_results(path: str) -> pd.DataFrame:
    kraken2_fields = [('is_classified', 'category'),
                      ('readID', str),
                      ('taxID', 'uint32'),
                      ('queryLength', 'uint16'),
                      ('LCA_mapping', str)]
    return pd.read_csv(path, sep='\t', header=None,
                       names=[k for k, v in kraken2_fields],
                       dtype={k: v for k, v in kraken2_fields}) \
        .set_index('readID')


def read_centrifuge_results(path: str) -> pd.DataFrame:
    centrifuge_results_dtypes = {
        'readID': str,
        'seqID': 'category',
        'taxID': 'uint32',
        'score': 'uint32',
        '2ndBestScore': 'uint32',
        'hitLength': 'uint16',
        'queryLength': 'uint16',
        'numMatches': 'uint8', }
    return pd.read_csv(path,
                       sep='\t',
                       dtype=centrifuge_results_dtypes) \
        .set_index('readID')


def write_reads_seqtk(reads_path: str, names: Iterable[str],
                      output_path: str) -> None:
    """Write reads with specified read names to an output file with seqtk and compress with pbgzip

    Using the following command-line
    `seqtk subseq reads.fq - | pbgzip -c > output.fq.gz`
    use seqtk to pull out a set of reads by name into a parallel block gzipped
    FASTQ file. Read names are provided via stdin.

    Args:
        reads_path: FASTQ file path
        names: read names
        output_path: Output block Gzipped FASTQ file
    Raises:
        subprocess.CalledProcessError: if the subprocess returns a non-zero exit code
        FileNotFoundError: if the output file doesn't exist
    """
    # TODO: make pbgzip optional? (pkruczkiewicz|2019-09-23 21:59:43.015)
    cmd = f'seqtk subseq {reads_path} - | pbgzip -c > {output_path}'
    p = sp.Popen(cmd,
                 stderr=sp.PIPE,
                 stdout=sp.PIPE,
                 stdin=sp.PIPE,
                 shell=True,
                 encoding='utf8')
    for name in names:
        print(name, file=p.stdin)
    stdout, stderr = p.communicate()
    if p.returncode > 0:
        raise sp.CalledProcessError(f'seqtk returned non-zero exit code '
                                    f'({p.returncode}) with command\n{cmd}\n'
                                    f'STDERR: {stderr}\n'
                                    f'STDOUT: {stdout}')
    if not os.path.exists(output_path):
        raise FileNotFoundError(f'seqtk subseq did not output file at '
                                f'"{output_path}" with command "{cmd}"')
