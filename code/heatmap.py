# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



df1 = pd.read_csv('.csv')
pd.set_option('max_columns', 10) #显示最多列数，超出该数以省略号表示
# 打印数据信息前5行的具体内容
#print(df1.head())
 
# 画热力图
plt.figure(figsize=(16, 9),dpi=300)
corr = df1.corr()
sns.heatmap(corr,cmap='Reds')
plt.show() 