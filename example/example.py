#!/usr/bin/env python
# coding: utf-8

"""
Example Karstolution run, using input and output files
"""


from __future__ import print_function, division, unicode_literals
import yaml
import pandas as pd
import matplotlib.pyplot as plt

# this lets us use karstolution without installing it
import sys
import os
sys.path.append(os.path.abspath('..'))

from Karstolution import karstolution


#
# ... specify file names
#
config_filename = 'config.yaml'
input_filename = 'input.csv'
output_filename = 'output__.csv'


#
# ... load inputs
#
config=yaml.safe_load(open(config_filename, 'rt').read())
df_input = pd.read_csv(input_filename)

#
# ... run model
#
model_output = karstolution(config, df_input, calculate_drip=True)


#
# ... save output to disk
#
model_output.to_csv(output_filename, index=False)
