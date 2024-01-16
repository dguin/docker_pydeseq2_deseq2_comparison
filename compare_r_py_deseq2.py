"""Compares current version of pyDESeq2 with R stable DESeq2

Usage:
  python3 -m analysis.sequencing.deseq2.compare_r_py_deseq2 \
    --data_directory_path=${DATA_DIR_PATH} \
    --output_directory_path=${OUPUT_DIR_PATH}
"""

import argparse
import os
import subprocess
from typing import Sequence

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

_MEAN_COUNT_COLNAME = 'baseMean'
_FC_COLNAME = 'log2FoldChange'
_Q_VALUE_COLNAME = 'padj'


def run_r_deseq2(data_directory_path: str, output_directory_path: str,
                 deseq2_script_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
  """Runs DESeq2 using R stable version.

  Args:
    data_directory_path: Path to folder containing input data.
    output_directory_path: Path to folder for output data.
    deseq2_script_path: R DESeq2 script path.

  Returns:
    Tuple of unshrunk and shrunk fold changes.
  """
  # Get outputs from R
  subprocess.run([
      'Rscript', deseq2_script_path, data_directory_path, output_directory_path
  ],
                 check=True)
  # Load R results
  r_unshrunk_fc = pd.read_csv(os.path.join(output_directory_path,
                                           'r_results_pre_shrunk.csv'),
                              index_col=0)
  r_shrunk_fc = pd.read_csv(os.path.join(output_directory_path,
                                         'r_results_post_shrunk.csv'),
                            index_col=0)

  return r_unshrunk_fc, r_shrunk_fc


def run_pydeseq2(
    data_directory_path: str,
    output_directory_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
  """Runs pyDESeq2 on a set of counts.

  Args:
    data_directory_path: Path to folder containing input data.
    output_directory_path: Path to folder for output data.

  Returns:
    Tuple of unshrunk and shrunk fold changes.
  """
  # Load counts and metadata
  counts = pd.read_csv(os.path.join(data_directory_path, 'counts.csv'),
                       index_col=0)
  metadata = pd.read_csv(os.path.join(data_directory_path, 'metadata.csv'),
                         index_col=0)
  counts = counts.T

  dds = DeseqDataSet(counts=counts,
                     metadata=metadata,
                     design_factors=list(metadata.columns),
                     refit_cooks=True,
                     ref_level=['condition', 'control'],
                     quiet=False)
  dds.deseq2()

  # Get fold changes from pydeseq2
  deseq_stats = DeseqStats(dds, quiet=True)
  deseq_stats.summary()
  py_unshrunk_fc = deseq_stats.results_df.copy(deep=True)

  # perform shrinkage
  deseq_stats.lfc_shrink()
  py_shrunk_fc = deseq_stats.results_df.copy(deep=True)

  py_unshrunk_fc.to_csv(
      os.path.join(output_directory_path, 'py_results_pre_shrunk.csv'))
  py_shrunk_fc.to_csv(
      os.path.join(output_directory_path, 'py_results_post_shrunk.csv'))
  pd.DataFrame(dds.obsm['size_factors'],
               index=counts.index,
               columns=['size_factor']).to_csv(
                   os.path.join(output_directory_path, 'py_size_factors.csv'))

  return py_unshrunk_fc, py_shrunk_fc


def main(data_directory_path: str, output_directory_path: str,
         deseq2_rscript_path: str) -> None:
  """Performs comparative runs using R DESeq2 and pyDESeq2 on count data.

  Args:
    data_directory_path: Folder path containing input data.
    output_directory_path: Folder to write output files in.
    deseq2_rscript_path: Path to R DESeq2 script.
  """

  if not os.path.exists(output_directory_path):
    os.makedirs(output_directory_path)

  # Run pydeseq2
  py_unshrunk_fc, py_shrunk_fc = run_pydeseq2(data_directory_path,
                                              output_directory_path)
  r_unshrunk_fc, r_shrunk_fc = run_r_deseq2(data_directory_path,
                                            output_directory_path,
                                            deseq2_rscript_path)

  # Merge datasets and format
  shrunk = py_shrunk_fc.join(r_shrunk_fc, lsuffix='_py', rsuffix='_r')
  unshrunk = py_unshrunk_fc.join(r_unshrunk_fc, lsuffix='_py', rsuffix='_r')


if __name__ == '__main__':

  parser = argparse.ArgumentParser(
      prog='CompareDESeq2',
      description='Fold change comparison generated from R DESeq2'
      ' and its pythonic counterpart pyDESeq2',
  )
  parser.add_argument(
      '--data_directory_path',
      metavar='i',
      type=str,
      required=True,
      dest='data_path',
      help=
      'Path for the input directory, should contain counts.csv and metadata.csv '
      'from analysis/sequencing/deseq2/set_up_deseq2_counts.py')
  parser.add_argument(
      '--output_directory_path',
      metavar='o',
      type=str,
      required=True,
      dest='output_path',
      help='Path for the output directory, directory will be made if it does not '
      'exist.')
  parser.add_argument('--deseq2_rscript_path',
                      metavar='r',
                      type=str,
                      default='deseq2.R',
                      dest='deseq2_rscript_path',
                      help='Path to DESeq2 Rscript.'
                      'exist.')

  args = parser.parse_args()
  main(args.data_path, args.output_path, args.deseq2_rscript_path)
