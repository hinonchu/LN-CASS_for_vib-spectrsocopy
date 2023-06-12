import os, glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

directory = r'Path\to\your\folder'
os.chdir(directory)
files = glob.glob('*.txt')



df  = pd.concat([pd.read_csv(file, sep = '\t', header = None, skiprows=[0], names = ['Wavenumber', file + ' Intensity']) for file in files], axis = 1)

df.to_csv('all_dat.csv')