'''
SYNBIOCHEM (c) University of Manchester 2019

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import os.path
import sys

from scipy.spatial.distance import pdist, squareform

import pandas as pd


def analyse_profiles(df):
    '''Analyse compound profiles.'''
    # Scale values to between 0 and 1:
    df = df / 100.00

    # Calculate Euclidean distances between all compound profiles:
    metric = 'euclidean'
    distances = pdist(df.values, metric=metric)
    dist_matrix = squareform(distances)

    # Update DataFrame with Euclidean distance between native and all compound
    # profiles:
    df[metric] = dist_matrix[0]

    return df.sort_values(metric)


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

    prof_out_df = analyse_profiles(prof_in_df)
    prof_out_df.to_csv(os.path.join(out_dir, 'Variant product profiles.csv'))


if __name__ == '__main__':
    main(sys.argv[1:])
