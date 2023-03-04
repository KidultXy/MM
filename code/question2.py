# 数据筛选 
decomposeddf = pd.DataFrame(decomposed_label.T)
decomposedata = pd.DataFrame(decomposed)

decomposeddf = decomposeddf.append(decomposedata)
decomposedata.columns = decomposed_label.T.tolist()

# GBDT
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor

bfx, bfy = decomposedata.iloc[:,:-1], decomposedata.iloc[:,-1]
bx_train, bx_test, by_train, by_test = train_test_split(bfx, bfy, test_size = 0.2, random_state = 0)
bfeat_labels = decomposedata.columns[:-1]
boost = GradientBoostingRegressor(n_estimators=500)
boost.fit(bfx, bfy)

bimportances = boost.feature_importances_
bindices = np.argsort(bimportances)[::-1]
for f in range(bfx.shape[1]):
    print("%2d) %-*s %f" % (f + 1, 30, bfeat_labels[bindices[f]][0], bimportances[bindices[f]]))

bfeat_labels_list = bfeat_labels[bindices].to_list()
bfeat_labels_array = np.array([it[0] for it in bfeat_labels_list])
bfeat_labels_chosen = bfeat_labels_array[:int(len(bfeat_labels_array)/10)]
bfeat_labels_chosen = bfeat_labels_chosen.tolist()
bfeat_labels_chosen.append(hyb)
bfeat_labels_chosen = np.array(bfeat_labels_chosen)
bfeat_labels_chosen.shape

ret = []
for it in bfeat_labels_chosen:
    ret.append(data[:,np.where(full_label==it)[0]])
print(ret)
gbdt28 = np.asarray(ret).reshape(int(len(bfeat_labels_array)/10)+1,-1).T
gbdt28.shape

dfgbdt28 = pd.DataFrame(gbdt28)
dfgbdt28.columns = bfeat_labels_chosen
dfgbdt28

decomposeddf.to_csv("gbdt28.csv",index=False,header=None,sep=",")

# Random Forest
from sklearn.ensemble import RandomForestRegressor

rfx, rfy = decomposedata.iloc[:,:-1], decomposedata.iloc[:,-1]
bx_train, bx_test, by_train, by_test = train_test_split(bfx, bfy, test_size = 0.2, random_state = 0)
rfeat_labels = decomposedata.columns[:-1]
rfmodel = RandomForestRegressor(n_estimators=1000)
rfmodel.fit(rfx, rfy)
rimportances = rfmodel.feature_importances_
rindices = np.argsort(rimportances)[::-1]
for f in range(bx_train.shape[1]):
    print("%2d) %-*s %f" % (f + 1, 30, rfeat_labels[rindices[f]], rimportances[rindices[f]]))


# XGB
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
params = {'learning_rate': 0.08, 'max_depth': 5, 'min_child_weight': 6, 'n_estimators': 100}

xgbx, xgby = decomposedata.iloc[:,:-1], decomposedata.iloc[:,-1]
xx_train, xx_test, xy_train, xy_test = train_test_split(xgbx, xgby, test_size = 0.2, random_state = 0)
# dtrain = xgb.DMatrix(xx_train, label = xy_train)
# dtest = xgb.DMatrix(xx_test, label = xy_test)
xgb_feat_labels = decomposedata.columns[:-1]

xgbmodel = xgb.XGBRegressor(learning_rate=0.08, max_depth=5, min_child_weight=6, n_estimators=100)
xgbmodel.fit(xgbx, xgby)

xgb_importances = xgbmodel.feature_importances_
xgb_indices  = np.argsort(xgb_importances)[::-1]
xgb_above_aver = []
xgb_dict = {}
for f in range(xx_train.shape[1]):
    if xgb_importances[xgb_indices[f]]>0.01:
        xgb_above_aver.append(decomposed_label[xgb_indices[f]][0])
    xgb_dict.update({ str(decomposed_label[xgb_indices[f]]):xgb_importances[xgb_indices[f]]})
    print("%2d) %-*s %f" % (f + 1, 30, decomposed_label[xgb_indices[f]], xgb_importances[xgb_indices[f]]))


plt.figure(figsize=(8,4))
plt.bar(range(len(xgb_importances)), xgb_importances)
plt.xlabel("Genes")
plt.ylabel("Feature Score")
plt.title("XGBoost -- Genes Importance")
plt.savefig("xbg1")
plt.show()

xgb_plot_x, xgb_plot_y = [],[]
for f in range(20):
    xgb_plot_x.append(decomposed_label[xgb_indices[f]][0])
    xgb_plot_y.append(xgb_importances[xgb_indices[f]])
xgb_plot_x.reverse()
xgb_plot_y.reverse()
plt.figure(figsize=(15,10))
plt.barh(xgb_plot_x, xgb_plot_y)
plt.ylabel("genes")
plt.xlabel("relative importance")
plt.title("XGBoost -- Genes Importance")
plt.savefig("xbg2")
plt.show()

# xgboost网格搜索调参
gsCv = GridSearchCV(xgbmodel,
                {'max_depth':list(range(4, 10, 1)),
                 'learning_rate':[0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2],
                 'min_child_weight':list(range(4, 8, 2)),
                 'n_estimators':list(range(10, 101, 20))}
gsCv.fit(xgbx, xgby)
print(gsCv.best_params_)
cv_results = pd.DataFrame(gsCv.cv_results_)
print(cv_results)

#full_label = np.array(full_label,dtype=str)
xgb_above_aver_pos = []
for sig in xgb_above_aver:
    print(sig)
    xgb_above_aver_pos.append(np.where(full_label==sig)[0][0])
len(xgb_above_aver_pos)

