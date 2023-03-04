# 数据导入 
import numpy as np
import pandas as pd

data = pd.read_excel("rat_eye.xlsx",header = None).to_numpy()
label = np.asarray(data[:,0],dtype=str)
mat = np.asarray(data[:,1:],dtype=float)
hyb = "1389163_at"
genpos = np.where(label == hyb)

# 转置与切分
x = np.delete(mat.T,genpos,axis=1)
x = np.asarray(x,dtype=float)

y = mat[genpos[0]].T
y = np.asarray(y,dtype=float)
label = np.delete(label,np.where(label == "1389163_at"))

data = np.column_stack((x,y))

# 标签调整
full_label = np.column_stack((label.reshape(1,-1),np.asarray("1389163_at")))
full_label = full_label.reshape(-1,1)


# KS检验
import matplotlib.pyplot as plt
sumrat = np.sum(data,axis=1)
# 散点图
plt.scatter(sumrat,range(1,121))
plt.show()

print(np.mean(sumrat)) # 均值
print(np.var(sumrat)) # 方差
print(np.std(sumrat)) # 标准差
mu = np.mean(sumrat)
sig = np.std(sumrat)


ratdf = pd.DataFrame(sumrat,columns = ['value'])

fig = plt.figure(figsize = (10,6))
ax1 = fig.add_subplot(2,1,1) # 创建子图1
ax1.scatter(ratdf.index, ratdf.values)
plt.grid() # # 绘制数据分布图
ax2 = fig.add_subplot(2,1,2) # 创建子图2
ratdf.hist(bins = 30,alpha = 0.5,ax = ax2)
ratdf.plot(kind = 'kde', secondary_y=True,ax = ax2)
plt.grid()

from scipy import stats
stats.kstest(ratdf['value'],'norm', (ratdf['value'].mean(), ratdf['value'].std()))

ks = []
for i in range(data.shape[1]):
    ratdf = pd.DataFrame(data[:,i],columns = ['value'])

    ksret = stats.kstest(ratdf['value'],'norm', (ratdf['value'].mean(), ratdf['value'].std()))
    ks.append(ksret)
    if i%1000==0:
        print(i)
    # .kstest方法：KS检验，参数分别是：待检验的数据，检验方法（这里设置成norm正态分布），均值与标准差
    # 结果返回两个值：statistic → D值，pvalue → P值# p值大于0.05，为正态分布
allgen_ks = [it[1] for it in ks]
allgen_ks = np.array(allgen_ks)
ks_pass_gens = np.where(allgen_ks>0.05)[0].tolist()
ks_fail_gens = list(set(range(0,18976)) ^ set(ks_pass_gens))

# 互信息
from sklearn.feature_selection import SelectKBest,mutual_info_regression

mut_selector = SelectKBest(score_func=mutual_info_regression, k=1000)
mut_selector.fit_transform(data[:,:-1], data[:,-1])

print('scores_:\n',mut_selector.scores_)
print('pvalues_:',mut_selector.pvalues_)
print('selected index:',mut_selector.get_support(True))
print('after transform:\n',mut_selector.transform(data[:,:-1]))

mut_selector.get_support(True)[:20]
data[:,mut_selector.get_support(True)] == mut_selector.transform(data[:,:-1])
mut_selected_gen = mut_selector.get_support(True)
mut_selector.scores_.shape

# Spearman
dim = data[1]
sample = data[0]

df = pd.DataFrame(data)
SP = df.corr("spearman")
print(SP.shape)
SPmatrix = SP.to_numpy()
print(SPmatrix[-1])

# MIC
import numpy as np
from minepy import MINE
 
MIC = []
mine = MINE(alpha=0.6, c=15)

for i in range(18976):
    mine.compute_score(data[:,i], data[:,-1])
    MIC.append(mine.mic())
mic_copy = MIC.copy()
mic_opt = []
for i in range(2000):
    pos = mic_copy.index(max(mic_copy))
    mic_opt.append(pos)
    del mic_copy[pos]

# 筛选
SPopt = np.where(np.abs(SPmatrix[-1])>0.5)
SPopt[0].shape
two_corr_ret = list(set(SPopt[0].tolist()) & set(mut_selected_gen))#& set(mic_opt)) #& set(f_selected_gen))# & set(Ropt[0].tolist())) 
two_corr_ret = sorted(two_corr_ret)
two_corr_ret.append(18975)
#print(two_corr_ret)
print(len(two_corr_ret))
len(set(two_corr_ret)&set(ks_fail_gens))

        # 根据相关性筛选
two_corr_ret = sorted(two_corr_ret)
decomposed = data[:,two_corr_ret]
decomposed.shape

decomposed_label = full_label[two_corr_ret]
decomposed_label.shape