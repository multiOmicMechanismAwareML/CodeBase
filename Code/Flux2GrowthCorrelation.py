import numpy as np
import pandas as pa
import scipy
from scipy.stats import pearsonr
from functools import partial

def check_file_read_ok(data, name):
    if data is None:
        print("error with " + name + " data not read")
    else:
        print("Success in loading " + name)
        print(data.shape)

def remove_zero_entry_columns(data):
    return data.loc[:, (data != 0).any(axis = 0)]

#Name of the row we are trying to predict (growth rate)
TARGET_NAME = 'log2relT'

#Load the data removing any
full_data = pa.read_csv('completeDataset.csv')
check_file_read_ok(full_data, "full data")
full_data = remove_zero_entry_columns(full_data)


expression_data = pa.read_csv('expressionOnly.csv')
check_file_read_ok(expression_data, "expression data")
#Extract the target and drop target column from main data
target_data = full_data[TARGET_NAME]

full_data = full_data.drop(columns=TARGET_NAME)
expression_data = expression_data.drop(columns=TARGET_NAME)

genes = full_data['Row']
full_data = full_data.drop(columns = 'Row')
reaction_data = full_data.drop(columns = expression_data.columns.values)
print(reaction_data.columns)
pearson = partial(pearsonr, target_data)
flux_correlation = reaction_data.apply(pearson, axis = 0)
flux_correlation = [l[0] for l in flux_correlation]
print(flux_correlation[0])
df = pa.DataFrame(flux_correlation, columns=['cor_with_growth'], index= list(reaction_data.columns))
df.transpose
df.to_csv('flux_correlation.csv', index=True, header=True, sep=' ')
