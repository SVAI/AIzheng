import numpy as np
import pandas as pd

def load_data(file_name):
    data = pd.read_csv(file_name, sep='\t')
    data.fillna(0.0, inplace=True)
    data_t = np.transpose(data.as_matrix()[:,1:])
    return data_t

filenames = ['stage_i_kich_data.txt', 'stage_i_kirc_data.txt', 'stage_i_kirp_data.txt',
             'stage_ii_kich_data.txt', 'stage_ii_kirc_data.txt', 'stage_ii_kirp_data.txt',
             'stage_iii_kich_data.txt', 'stage_iii_kirc_data.txt', 'stage_iii_kirp_data.txt',
             'stage_iv_kich_data.txt', 'stage_iv_kirc_data.txt', 'stage_iv_kirp_data.txt']

data = load_data('healthy_kidney_data.txt')

for fname in filenames[1:]:
    new_data = load_data(fname)
    data = np.concatenate((data, new_data), axis=0)

np.save('kidney-data.npy',data)
