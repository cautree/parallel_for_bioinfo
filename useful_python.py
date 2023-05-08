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

## longer index of a string will not give an error
project = "long_sequence"
project[:31] # will not give an error


## save a dataframe into different sheet of excels based on a categorical variable, such as project here

writer = pd.ExcelWriter(output_file)
for project in sorted([p for p in set(raw.Project) if p != 'default']):
    df = raw[raw.Project == project].reset_index(drop = True)
    df.to_excel(writer, sheet_name = project[:31], index = False)
writer.close()

## get summary stats by group, and get the total to merge [concat] into the same data frame
output_stats = 'seq_stats.txt'

stats = raw.groupby('Project')['Yield'].sum().reset_index()
stats['Project'] = stats.Project.apply(lambda x: 'Undetermined' if x=='default' else x)
sum_stats = pd.DataFrame([['Total', stats['Yield'].sum()]], columns = ['Project', 'Yield'])
## concat, vertically
stats = pd.concat([stats, sum_stats], ignore_index=True)
stats['Yield'] = stats['Yield']/1000
stats['Yield'] = stats.Yield.apply(lambda x: '%.3f' % x) + 'Gb'
stats[['Project', 'Yield']].to_csv(output_stats, index = False, header = False, sep = '\t')



## create an empty python data frame, with column names defined
temp = pd.DataFrame( columns = ["a","b","c","d"])

## get the unique values for a specific columns, return a list
df.Sample.unique()

## get the string length groups for a specific column, like there are two groups, one group of primer 18bp, and another 20bp
df.primer.str.len().unique()

## check if the value of a column start with something
primers=df.primer
primers.startswith("A")
            

## check if there are missing value in the groups, if there is missing values in any column, the count value is smaller
## this returns a data frame
df.groupby('Sample_Plate').count()


## add an addtioanl row at the end of a dataframe, the row index is "Empty"
df.loc['Empty'] = ''


## check for a combined duplicates in two columns
if df.duplicated(subset=['a', 'b']).any():
    raise ValueError('a b pair duplicated')
