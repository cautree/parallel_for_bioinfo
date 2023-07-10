import pandas as pd
import numpy as np
def get_mean(x,y,z):
    m = np.mean([x,y,z])
    return(m)
##four different ways to create a new columns using three columns in the data  
df_for_avg["avg"] = df_for_avg.apply( lambda x: (x.math + x.physics+ x.chemistry)/3, axis =1)

df_for_avg["avg2"] = df_for_avg.apply( lambda x: np.mean([x.math , x.physics, x.chemistry]), axis =1)


df_for_avg["avg3"] = [(x+y+z)/3 for x,y,z in  zip(df_for_avg.math, df_for_avg.physics, df_for_avg.chemistry)  ]

df_for_avg["avg4"] = [get_mean(x,y,z) for x,y,z in  zip(df_for_avg.math, df_for_avg.physics, df_for_avg.chemistry)  ]



import pathlib2 as pl2
ps = pl2.Path('data/sp3')

dfs = (
    pd.read_csv(p, encoding='utf8') for p in ps.glob('*.csv')
)

res = pd.concat(dfs)
res


## create two new columns using one current column
df_test["first_last"] = df_test.apply( lambda x: x.Name.split(', '), axis =1 )
df_test[['first', 'last']] = pd.DataFrame(df_test["first_last"].tolist(), index=df_test.index)


## create a flag column
df_test["age_na"] = [ 1 if pd.isna(x) else 0 for x in df_test.Age]

def get_age(x,y,z): # x Age, y age_na, z Pclass
    if y== 0:
        x = x
    else:
        if z == 1:
            x = 37
        elif z ==2:
            x= 29
        else:
            x = 24
    return x
## impute based on two column info
df_test["age_imputated"] = df_test.apply( lambda x: get_age(x.Age, x.age_na, x.Pclass), axis =1)



## create a category group when the categories are more than 2
def create_age_groups( age):
    if age <=12:
        age_cat = "child"
    elif age <=18:
        age_cat = "Teen"
    elif age <=60:
        age_cat = "Adult"
    else:
        age_cat = "Senior"
    return age_cat


df_test["age_group"] = df_test.apply(lambda x: create_age_groups(x.age_imputated) , axis =1)
df_test.age_group.value_counts( dropna = False )


## another to it is to use pd.cut
bins=[0, 13, 19, 61, sys.maxsize]
labels=['<12', 'Teen', 'Adult', 'Older']
df_test[ "age_cat"] =  pd.cut(df_test.Age, bins = bins, labels = labels)
df_test


# A lambda expression for Standardization.
standardization = lambda x: (x - x.mean()) / x.std()
mtcars_s = mtcars.select_dtypes(include = 'float64')
mtcars_sd = mtcars_s.apply( standardization)
mtcars_sd = mtcars_s.transform( standardization)











