# 第二轮筛选 
xgb_above_aver_pos.append(18975)
xgb_above_aver_pos = sorted(xgb_above_aver_pos)
xgb_decomposed = data[:,xgb_above_aver_pos]

xgb_decomposed_label = full_label[xgb_above_aver_pos]

xgb_decomposed_data = pd.DataFrame(xgb_decomposed)
xgb_decomposed_data

xgb_decomposed_data.to_csv("26genes.csv")

xgbx2, xgby2 = xgb_decomposed_data.iloc[:,:-1], xgb_decomposed_data.iloc[:,-1]

# KPCA降维
from KPCA import *
from sklearn.decomposition import KernelPCA
from sklearn.metrics import mean_squared_error,r2_score

kernels = ["rbf", "cosine"]
mse = list()  # 同一核函数不同主成分下的均方误差
n = 120

for kernel in kernels:
    for i in range(1,25):
        kpca = KernelPCA(n_components=int(i), kernel=kernel, fit_inverse_transform=True)
        X_reduced = kpca.fit_transform(xgbx2)
        X_back = kpca.inverse_transform(X_reduced)
        R_squ = 1-r2_score(X_back, xgbx2)
        #R_mean = mean_squared_error(X_back, xgbx2)

        R_adj = 1 -(1-R_squ)*(n-1)/(n-i-1)
        mse.append(R_adj)

print(len(mse))
mse_linear = mse[:int(len(mse)/2)]  #线性核函数的均方误差
mse_cosine = mse[int(len(mse)/2):]  #余弦核函数的均方误差

plt.bar(range(int(len(mse)/2)),mse_linear)
plt.ylabel("Adjusted_R2")
plt.xlabel("number of principal components")
plt.title("kPCA -- selected gene numbers")
plt.savefig("kpca1")
plt.show()

step = (np.array(mse_linear[:-1]) - np.array(mse_linear[1:]))
plt.bar(range(len(step)), step)
plt.ylabel("change rate of Adjusted_R2")
plt.xlabel("number of principal components")
plt.title("kPCA -- selected gene numbers")
plt.savefig("kpca2")
plt.show()


# RFE
from sklearn.feature_selection import RFE

xgbx2, xgby2 = xgb_decomposed_data.iloc[:,:-1], xgb_decomposed_data.iloc[:,-1]
x_slt = xgbx2
y_slt = xgby2

M = xgbx2.shape[1]

for n in range(1,M-10):
    print(x_slt.shape[1])
    selector = RFE(xgbmodel, n_features_to_select=M-n,step=1).fit(x_slt,y_slt)
    slt = selector.support_

    F_slt = np.where(slt==True)
    x_slt = x_slt[x_slt.columns[F_slt]]

x_slt.columns.to_list()

rfe_decomposed_label = xgb_decomposed_label[x_slt.columns.to_list()]
rfe_decomposed_label = np.array(rfe_decomposed_label,dtype=str)
rfe_decomposed_dict = {}
for it in rfe_decomposed_label:
    rfe_decomposed_dict.update({it[0]:xgb_dict[str(it)]})

rfe_decomposed_sorted = sorted(rfe_decomposed_dict.items(),key=lambda x:x[1],reverse=True)

rfe_plot_x, rfe_plot_y = [],[]
for f in range(10):
    rfe_plot_x.append(rfe_decomposed_sorted[f][0])
    rfe_plot_y.append(rfe_decomposed_sorted[f][1])
rfe_plot_x.reverse()
rfe_plot_y.reverse()
plt.figure(figsize=(10,6))
plt.barh(rfe_plot_x, rfe_plot_y)
plt.ylabel("genes")
plt.xlabel("relative importance")
plt.title("RFE -- Genes Importance")
plt.savefig("rfe")
plt.show()

# RFECV
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

xgbx2, xgby2 = xgb_decomposed_data.iloc[:,:-1], xgb_decomposed_data.iloc[:,-1]
# estimator = SVR(kernel="linear")

# 5折交叉
selector = RFECV(estimator=xgbmodel, step=1, cv=5)
selector = selector.fit(xgbx2, xgby2)

# 哪些特征入选最后特征，true表示入选
print(selector.support_)

# 每个特征的得分排名，特征得分越低（1最好），表示特征越好
print(selector.ranking_)

#  挑选了几个特征
print(selector.n_features_)
# 每次交叉迭代各个特征得分
print(selector.grid_scores_)

plt.plot(range(25),selectorr.grid_scores_)
plt.xlabel("numbers of features selected")
plt.ylabel("cross validation score")
plt.title("RFECV -- cross validation score")
plt.savefig("cv_score")

rfecv_decomposed_label = xgb_decomposed_label[:-1][selector.support_]
rfecv_decomposed_label = np.array(rfecv_decomposed_label,dtype=str).tolist()

rfecv_decomposed_dict = {}
for it in rfecv_decomposed_label:
    rfecv_decomposed_dict.update({it[0]:xgb_dict[str(it)]})

rfecv_decomposed_sorted = sorted(rfecv_decomposed_dict.items(),key=lambda x:x[1],reverse=True)

rfecv_plot_x, rfecv_plot_y = [],[]
for f in range(10):
    rfecv_plot_x.append(rfecv_decomposed_sorted[f][0])
    rfecv_plot_y.append(rfecv_decomposed_sorted[f][1])
rfecv_plot_x.reverse()
rfecv_plot_y.reverse()
plt.figure(figsize=(10,6))
plt.barh(rfecv_plot_x, rfecv_plot_y)
plt.ylabel("genes")
plt.title("RFECV -- Genes Importance")
plt.xlabel("feature importance")
plt.savefig("rfecv")
plt.show()


# 热力图
import seaborn as sns

heat_data = np.array(data[:,-1])

for gene in rfecv_plot_x:
    print(gene)
    heat_data = np.column_stack((heat_data,data[:,np.where(full_label==gene)[0]]))

heat_data.shape
heat_df = pd.DataFrame(heat_data)
heat_label = ["hyb"]+list(range(0,14))


plt.figure(figsize=(10,6))
sns.heatmap(heat_df.corr("spearman"),cmap='Reds')
plt.title("heat map for Spearman rank correlation")
plt.savefig("heatmap_sp")
plt.show()


#　表格
spline = heat_df.corr("spearman").to_numpy()[0][1:]
spline = spline.reshape(1,-1)

mutline = []
for gene in rfecv_plot_x[:-1]:
    mutline.append(mut_selector.scores_[np.where(full_label==gene)[0]])

mutline = np.array(mutline).reshape(1,-1)

micline = []
for gene in rfecv_plot_x[:-1]:
    micline.append(MIC[int(np.where(full_label==gene)[0])])

micline = np.array(micline).reshape(1,-1)

table_four_line = pd.DataFrame()
table_four_line = table_four_line.append(pd.DataFrame(spline))
table_four_line = table_four_line.append(pd.DataFrame(mutline))
table_four_line = table_four_line.append(pd.DataFrame(micline))
table_four_line.columns = rfecv_plot_x[:-1]
table_four_line.index = ["Spearman","MI","MIC"]

table_four_line.to_csv("four.csv")