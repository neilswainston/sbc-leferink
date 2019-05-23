'''
SYNBIOCHEM (c) University of Manchester 2019

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import itertools
import os.path
import sys

from scipy.spatial.distance import pdist, squareform

import numpy as np
import pandas as pd


def analyse_comp_profiles(df):
    '''Analyse compound profiles.'''
    # Scale values to between 0 and 1:
    df = df / 100.00

    # Calculate Euclidean distances between all compound profiles:
    metric = 'euclidean'
    distances = pdist(df.values, metric=metric)
    dist_matrix = squareform(distances)

    # Update DataFrame with Euclidean distance between VAR3 and all compound
    # profiles:
    df[metric] = dist_matrix[1]  # 1 is index of VAR3 in data

    return df.sort_values(metric)


def analyse_plasticity(df):
    '''Analyse positional plasticity.'''
    # Scale values between 0 and 1:
    df = df / df.sum()

    # Consolidate all residue frequencies, ignoring codons:
    df.index = df.index.droplevel(level=0)
    df = df.groupby(level=0).sum().reset_index()
    df.index = df[('Residue', '', '')]
    df.index.name = 'Residue'
    df = df.drop(('Residue', '', ''), axis=1)

    # Remove silent mutations:
    # df = _remove_silent_mutations(df)

    # Calculate plasticity, as overlap between distribution and
    # uniform distribution:
    uniform = [1 / len(df)] * len(df)

    name = 'plasticity'
    return df.apply(lambda col: _get_intersection(col, uniform)
                    ).to_frame(name=name)


def _get_intersection(dist1, dist2):
    '''Get intersection between distributions.'''
    minima = np.minimum(dist1, dist2)
    return np.true_divide(np.sum(minima), np.sum(dist2))


def _remove_silent_mutations(df):
    '''Remove silent mitations.'''
    var3 = df.columns.get_level_values('VAR3')

    sim_df = pd.DataFrame(
        np.array(
            [val[0] != val[1] for val in itertools.product(df.index, var3)]
        ).reshape(len(df.index), len(var3)))

    return pd.DataFrame(df.values * sim_df.values,
                        columns=df.columns, index=df.index)


def main(args):
    '''main method.'''
    out_dir = args[0]
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(curr_dir, 'data')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Calculate distances between native and variant compound profiles:
    prof_in_df = pd.read_excel(os.path.join(data_dir,
                                            'Variant product profiles.xlsx'),
                               sheet_name='Relative_incl ger',
                               index_col='Variant')

    prof_out_df = analyse_comp_profiles(prof_in_df)
    prof_out_df.to_csv(os.path.join(out_dir, 'Variant product profiles.csv'))

    # Analyse position plasticity:
    aa_in_df = pd.read_excel(os.path.join(data_dir,
                                          'Amino acid occurence.xlsx'),
                             sheet_name='Sheet1',
                             header=[0, 1, 2],
                             index_col=[0, 1])

    aa_out_df = analyse_plasticity(aa_in_df)
    aa_out_df.to_csv(os.path.join(out_dir, 'plasticity.csv'))


if __name__ == '__main__':
    main(sys.argv[1:])
