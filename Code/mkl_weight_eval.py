import pandas as pa
import matplotlib.pyplot as plts
import seaborn as sns
import os.path
mm_ge_flu_weights = pa.read_csv('MM-GE(1:6)-Flu(7:12)_Weights.csv')
mm_gem_flu_weights = pa.read_csv('MM-GEM(1:6)-Flu(7:12)_Weights.csv')
sns.set(color_codes = True)
sns.set_style('white')

print("ge expression : " , mm_ge_flu_weights[1:6].sum() / mm_ge_flu_weights.sum())
print("ge fluxes : " , mm_ge_flu_weights[7:12].sum() / mm_ge_flu_weights.sum())
print("gem expression : " , mm_gem_flu_weights[1:6].sum() / mm_gem_flu_weights.sum())
print("gem fluxes : " , mm_ge_flu_weights[7:12].sum() / mm_gem_flu_weights.sum())
