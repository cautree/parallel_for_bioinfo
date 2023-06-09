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
output_stats = 'stats.txt'

stats = raw.groupby('group')['count'].sum().reset_index()
stats['group'] = stats.Project.apply(lambda x: 'Undefined' if x=='.' else x)

## create a smaller data frame with the same columns to concat with the stats df above
sum_stats = pd.DataFrame([['Total', stats['count'].sum()]], columns = ['group', 'count'])
#  group	count
#0	Total	5072

## concat, vertically
stats = pd.concat([stats, sum_stats], ignore_index=True)
stats['count'] = stats['count']/1000
stats['count'] = stats.count.apply(lambda x: '%.3f' % x) + 'oz'
stats[['group', 'count']].to_csv(output_stats, index = False, header = False, sep = '\t')



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
    
    
## check how many groups are there after groupby
df.groupby(['a', 'b']).ngroups


## save each column of data into a csv file, transforming long data into wide data using pivot
x = '20220504_sequencerA'
date, sequencer = x.split('_')
df = pd.read_excel('{}_{}.xlsx'.format(date, sequencer), sheet_name='a')

df['forward'] = df.Sample.str.split('_').str[1]
df['reverse'] = df.Sample.str.split('_').str[2]

for i in df.Project.unique():
    if i == 'default':
        continue
    else:
        temp_df=df[df.Project==i]
        names = zip(['% A', '% B', '$ C'],
                    ['A', 'B', 'C'])
    for val, name in names:
        out_path = '{:}_{:}_{:}_{:}.sorted.csv'.format(date, sequencer, i, name)
        tmp = temp_df.groupby(['forward', 'reverse'])[val].sum().reset_index().pivot(index='forward', columns='reverse', values=val)
        tmp=tmp[tmp.index.notnull()]
        tmp=tmp.loc[:, tmp.columns.notnull()]
        tmp.to_csv(out_path)
        
## python way to read merge csv file into one big file
out_file =  'merged.csv'

paths = [path for path in sorted(os.listdir('.')) if path.endswith('.csv')]

df = pd.DataFrame(range(1,1001), columns = ['count'])

for path in paths:
    samp = path.replace('.csv','')

    try:
        tmp = pd.read_csv(path)
        tmp.columns = ['count', samp]
        df = df.merge(tmp, how = 'left').fillna(0)
    except:
        tmp = pd.DataFrame(range(1,1001), columns = ['count'])
        tmp[samp] = 0

df.to_csv(out_file, float_format='%.0f', index = False)



## one metric a excel sheet, save data into one excel file
import os
import pandas as pd

paths = [ "metrics/" + path for path in sorted(os.listdir('metrics')) if path.endswith('.txt')]
output_path = "metrics.xlsx"

# define the output file name
writer = pd.ExcelWriter(output_path)

# create list of file extensions and dict of sheetnames for excel file
file_groups = ['.metricA.txt',
           '.metricB.txt']

sheet_ID = ['CollectMetricA',
              'CollectMetricB']

sheet_dict = dict(zip(file_groups, sheet_ID))

# iterate through file_groups 
for group in file_groups:
    df = pd.DataFrame([])
    for path in [path for path in paths if path.endswith(ending)]:
        # for each path within the file group, create sample name
        sample = path.replace(group,'').split('/')[-1

        # read in data from that path
        try:
            tmp = pd.read_csv(path, delimiter='\t', skiprows=6, nrows=1)
            tmp = tmp.T
            tmp.columns = [sample]
        except:
            tmp = pd.DataFrame([], columns = [sample])

        # append data to dataframe, notice it is full join, both left index and right index kept
        # for the first round, df is empty, then tmp file's index becomes the merged file index                                            
        df = df.merge(tmp, left_index=True, right_index=True, how='outer')

    if len(df) > 0:
        df.T.sort_index().to_excel(writer, sheet_name = sheet_dict[group][:30])

writer.close()


## fillna, use a dictionary to fill na
#    A    B   C    D
#0  NaN  2.0 NaN  0.0
#1  3.0  4.0 NaN  1.0
#2  NaN  NaN NaN  NaN
#3  NaN  3.0 NaN  4.0

#In [136]:df.index
#Out[136]: RangeIndex(start=0, stop=4, step=1)
df.C.fillna(dict(zip([0,1,2,3],[0,10,20,30])  
                 
                 
## parse excel files
path = "fileA.xlsx"
excel = pd.ExcelFile(path)
for sheet in excel.sheet_names:
    df = excel.parse(sheet)
                 
                 
## if there is ? in the report for numerical variable, how to replace it as 0, then calculate the mean
df.A.astype(str).str.replace('?','0').astype(float).mean())
            
            
## str connected
 tmp.metric_name.str.lstrip().str.replace('|','')
            
            
##### python for loop and while loop work together            
from scripts.FASTA import ReadFASTA

dna, sub_seq = [fasta[1] for fasta in ReadFASTA('data/rosalind_sseq_v2.txt')]
print(dna)
print(sub_seq)

#sseq_indicies, i = [], 0
sseq_indicies = []
i = 0

for nucleotide in sub_seq:
    # In practice: Use exception handling/additional constraints as such a subsequence does not necessarily exist.
    while dna[i] != nucleotide:
        print("this is in while")
        print ("this is nucleotide in sub_seq {} ".format(nucleotide))
        print(dna[i])
        print(i)
        
        i += 1
        

    # Use i+1 as the indicies because Rosalind starts at i=1 instead of i=0.
    sseq_indicies.append(str(i+1))
    i += 1
    print("this is the position (1 based) of nucletodie that matched to the subseq: {}". format(i))

print (' '.join(sseq_indicies))
with open('output/030_SSEQ.txt', 'w') as output_data:
    output_data.write(' '.join(sseq_indicies))
            
            
##Check if the nucleotides are in the same purine/pyrimidine group.
[['A', 'G'],['C', 'T']][dna2[i] in ['C', 'T']]

[['A', 'G'],['C', 'T']][True]   ## both are in CT
[['A', 'G'],['C', 'T']][False]  ## both are in AG

            
## use filter with lambda function, and reduce with lammda function
dna_groups_new = []
a = list(filter(lambda x: len(x)>1, dna_groups))
b = list(filter(lambda x: len(x)==1, dna_groups))
from functools import reduce
b = reduce(lambda x,y:x+y, b)
dna_groups_new.append(a)
dna_groups_new.append(b)



# Merge two dictionaries with a comprehension
team1 = {"Jones": 24, "Jameson": 18, "Smith": 58, "Burns": 7}
team2 = {"White": 12, "Macke": 88, "Perce": 4}
newTeam = {k: v for team in (team1, team2) for k, v in team.items()}
print(newTeam)



newTeam_1 = {k: v  for k, v in team1.items()}
newTeam_2 = {k: v  for k, v in team2.items()}
new = newTeam_1 | newTeam_2
print(new)




# define a list of temperature data points
ctemps = [5, 10, 12, 14, 10, 23, 41, 30, 12, 24, 12, 18, 29]

# build a set of unique Fahrenheit temperatures
ftemps1 = [(t * 9/5) + 32 for t in ctemps]
print(ftemps1)
ftemps2 = {(t * 9/5) + 32 for t in ctemps}
print(ftemps2)

# build a set from an input source
s_temp = "The quick brown fox jumped over the lazy dog"
chars = {c.upper() for c in s_temp if not c.isspace()}
print(chars)



import string
import pprint


test_str = "2 apples, 9 oranges?, 4 pears, Mike's 1 egg, Jane's 2 kiwis, $50!"

# get the total length of the string
l = len(test_str)

# count the number characters
nums = len([c for c in test_str if c.isnumeric()])

# count the punctuation characters
punct = len([c for c in test_str if c in string.punctuation])

# use a set to count the unique letters
unique = "".join({c for c in test_str if c.isalpha()})

# print the data
str_data = {
    "Length:": l,
    "Digits:": nums,
    "Punctuation": punct,
    "Unique Letters": unique,
    "Unique Count": len(unique)
}
pprint.pp(str_data)
print(str_data)
