import pandas as pd
df = pd.read_csv('datasets/all_stocks_2017-2018.csv')

## check the columns names
df.columns

## make a Series into a dataframe
df.Name.value_counts().to_frame().reset_index()

## filter based on a column value, then reset the index
df_tech = df[df.Name.isin(['AAPL','GOOGL','MSFT','AMZN'])].reset_index()

#Pivot dataframe, so that each stock is a column, more like tidyr::spread in R
#only showing the info for Close, not Open, High Low
pivoted_df = df_tech.pivot(index='Date', columns='Name', values='Close')

# Repeat the shifting and calculation for the 14, 21, 28-day offset
# calculate returns over each window (offset) and store them in a dictionary with a loop
# a smart way, and create the dict name on the go
delta_dict = {}
for offset in [7, 14, 21, 28]:
    delta_dict['delta_{}'.format(offset)] = pivoted_df / pivoted_df.shift(offset) - 1.0
    
    
# df.columns 
# brand speed
# df.brand.name will give you "brand"


## combine newly created columns with previously exist columns in dataframe
cols = df.columns
df = df[list(cols) + ['well', 'row_sort', 'column_sort', '% of the plate']]
