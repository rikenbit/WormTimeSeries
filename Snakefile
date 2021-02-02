import pandas as pd

configfile: "config.yaml"

# WTS1 sample_sheet
SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')

#

# print(SAMPLE_SHEET)