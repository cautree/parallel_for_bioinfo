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












