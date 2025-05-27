import time
import os
import pandas as pd
from sklearn.metrics import accuracy_score, precision_recall_curve, precision_score, recall_score, f1_score, roc_auc_score, classification_report, roc_curve, confusion_matrix
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def print_result(y_test, y_pred):
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
    print('accuracy_score:', accuracy_score(y_test, y_pred))
    print('precision_score:', precision_score(y_test, y_pred))
    print('recall_score:', recall_score(y_test, y_pred))
    print('f1_score:', f1_score(y_test, y_pred))
    print('roc_auc_score:', roc_auc_score(y_test, y_pred))

# 打印混淆矩阵
def print_confusion_matrix(y_test, y_pred):
    from sklearn.metrics import confusion_matrix
    print(confusion_matrix(y_test, y_pred))

# 打印分类报告
def print_classification_report(y_test, y_pred):
    from sklearn.metrics import classification_report
    print(classification_report(y_test, y_pred))

def get_scores(y_test,y_pred,y_pred_prob):
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
    scores = []
    scores.append(accuracy_score(y_test, y_pred))
    scores.append(precision_score(y_test, y_pred))
    scores.append(recall_score(y_test, y_pred))
    scores.append(f1_score(y_test, y_pred))
    try:
        scores.append(roc_auc_score(y_test, y_pred_prob[:,1]))
    except:
        scores.append(roc_auc_score(y_test, y_pred))
    return scores

# 新增函数：绘制混淆矩阵
def plot_confusion_matrix(y_test, y_pred, save_path):
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# 新增函数：绘制ROC曲线
def plot_roc_curve(y_test, y_pred_prob, model_names, save_path):
    plt.figure(figsize=(10, 8))
    
    for i, (y_prob, name) in enumerate(zip(y_pred_prob, model_names)):
        if isinstance(y_prob, np.ndarray) and y_prob.ndim > 1:
            # 对于返回概率矩阵的模型，取正类的概率
            y_prob = y_prob[:, 1]
        
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        auc = roc_auc_score(y_test, y_prob)
        plt.plot(fpr, tpr, lw=2, label=f'{name} (AUC = {auc:.3f})')
    
    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def test_xgboost(X_train, X_test, y_train, y_test, sign):
    import xgboost as xgb
    
    # 使用更优的参数配置
    xgb_clf = xgb.XGBClassifier(
        n_estimators=100,           # 增加树的数量
        learning_rate=0.1,          # 较小的学习率
        max_depth=5,                # 控制树的深度
        min_child_weight=1,         # 避免过拟合
        gamma=0,                    # 最小损失减少
        subsample=0.8,              # 随机采样训练样本
        colsample_bytree=0.8,       # 随机采样特征
        reg_alpha=0.01,             # L1正则化
        reg_lambda=1,               # L2正则化
        random_state=42,            # 随机种子
        objective='binary:logistic' # 二分类目标函数
    )
    
    # 使用早停策略
    xgb_clf.fit(
        X_train, y_train,
        eval_set=[(X_test, y_test)],
        early_stopping_rounds=10,
        verbose=False
    )
    
    y_pred = xgb_clf.predict(X_test)
    y_pred_prob = xgb_clf.predict_proba(X_test)
    scores = get_scores(y_test, y_pred, y_pred_prob)
    
    return xgb_clf, y_pred, y_pred_prob, scores

def test_random_forest(X_train, X_test, y_train, y_test,sign):
    from sklearn.ensemble import RandomForestClassifier
    rf_clf = RandomForestClassifier(n_estimators=100, random_state=42,n_jobs=-1)
    rf_clf.fit(X_train, y_train)
    y_pred = rf_clf.predict(X_test)
    y_pred_prob = rf_clf.predict_proba(X_test)
    scores = get_scores(y_test,y_pred,y_pred_prob)
    
    return rf_clf, y_pred, y_pred_prob, scores

def test_MLP(X_train, X_test, y_train, y_test,sign):
    from sklearn.neural_network import MLPClassifier
    mlp_clf = MLPClassifier(random_state=42)
    mlp_clf.fit(X_train, y_train)
    y_pred = mlp_clf.predict(X_test)
    y_pred_prob = mlp_clf.predict_proba(X_test)
    scores = get_scores(y_test,y_pred,y_pred_prob)

    return mlp_clf, y_pred, y_pred_prob, scores

def test_SVM(X_train, X_test, y_train, y_test,sign):
    from sklearn.svm import SVC
    svm_clf = SVC(random_state=42, probability=True)  # 添加probability=True以获取概率
    svm_clf.fit(X_train, y_train)
    y_pred = svm_clf.predict(X_test)
    y_pred_prob = svm_clf.predict_proba(X_test)  # 使用predict_proba替代decision_function
    scores = get_scores(y_test,y_pred,y_pred_prob)
    return svm_clf, y_pred, y_pred_prob, scores

def test_logistic_regression(X_train, X_test, y_train, y_test,sign):
    from sklearn.linear_model import LogisticRegression
    lr_clf = LogisticRegression(random_state=42)
    lr_clf.fit(X_train, y_train)
    y_pred = lr_clf.predict(X_test)
    y_pred_prob = lr_clf.predict_proba(X_test)
    scores = get_scores(y_test,y_pred,y_pred_prob)
    
    return lr_clf, y_pred, y_pred_prob, scores



def run_all_test(features_i,labels_all, train_index, test_index,sign):
    X_train = features_i[train_index]
    X_test = features_i[test_index]
    y_train = labels_all[train_index]
    y_test = labels_all[test_index]
    
    xgb_clf, y_pred_xgb, y_pred_prob_xgb, scores_xgb = test_xgboost(X_train, X_test, y_train, y_test,sign)
    rf_clf, y_pred_rf, y_pred_prob_rf, scores_rf = test_random_forest(X_train, X_test, y_train, y_test,sign)
    mlp_clf, y_pred_mlp, y_pred_prob_mlp, scores_mlp = test_MLP(X_train, X_test, y_train, y_test,sign)
    svm_clf, y_pred_svm, y_pred_prob_svm, scores_svm = test_SVM(X_train, X_test, y_train, y_test,sign)
    lr_clf, y_pred_lr, y_pred_prob_lr, scores_lr = test_logistic_regression(X_train, X_test, y_train, y_test,sign)
    
    scores = pd.DataFrame([scores_xgb, scores_rf, scores_mlp, scores_svm, scores_lr], columns=['accuracy_score', 'precision_score', 'recall_score', 'f1_score', 'roc_auc_score'], index=['xgboost', 'random_forest', 'MLP', 'SVM', 'logistic_regression'])
    
    # 以当前时间戳为文件名
    t = time.time()
    t = int(t)
    
    test_dir = './'+sign+'_'+str(t)

    os.makedirs(test_dir, exist_ok=True)
    
    # 绘制并保存混淆矩阵
    model_names = ['xgboost', 'random_forest', 'MLP', 'SVM', 'logistic_regression']
    y_preds = [y_pred_xgb, y_pred_rf, y_pred_mlp, y_pred_svm, y_pred_lr]
    
    for name, y_pred in zip(model_names, y_preds):
        plot_confusion_matrix(
            y_test, 
            y_pred, 
            save_path=f"{test_dir}/{sign}_{name}_confusion_matrix.png"
        )
    
    # 绘制并保存ROC曲线
    y_pred_probs = [y_pred_prob_xgb, y_pred_prob_rf, y_pred_prob_mlp, y_pred_prob_svm, y_pred_prob_lr]
    plot_roc_curve(
        y_test, 
        y_pred_probs, 
        model_names, 
        save_path=f"{test_dir}/{sign}_roc_curves.png"
    )
    
    # 保存模型与y_pred,y_test
    import joblib
    joblib.dump(xgb_clf, test_dir + '/' + sign+'_xgb_clf.pkl')
    joblib.dump(rf_clf, test_dir + '/' + sign+'_rf_clf.pkl')
    joblib.dump(mlp_clf, test_dir + '/' + sign+'_mlp_clf.pkl')
    joblib.dump(svm_clf, test_dir + '/' + sign+'_svm_clf.pkl')
    joblib.dump(lr_clf, test_dir + '/' + sign+'_lr_clf.pkl')
    
    np.save(test_dir + '/' + sign+'_y_pred_xgb.npy', y_pred_xgb)
    np.save(test_dir + '/' + sign+'_y_pred_rf.npy', y_pred_rf)
    np.save(test_dir + '/' + sign+'_y_pred_mlp.npy', y_pred_mlp)
    np.save(test_dir + '/' + sign+'_y_pred_svm.npy', y_pred_svm)
    np.save(test_dir + '/' + sign+'_y_pred_lr.npy', y_pred_lr)
    
    # 保存概率预测结果
    np.save(test_dir + '/' + sign+'_y_pred_prob_xgb.npy', y_pred_prob_xgb)
    np.save(test_dir + '/' + sign+'_y_pred_prob_rf.npy', y_pred_prob_rf)
    np.save(test_dir + '/' + sign+'_y_pred_prob_mlp.npy', y_pred_prob_mlp)
    np.save(test_dir + '/' + sign+'_y_pred_prob_svm.npy', y_pred_prob_svm)
    np.save(test_dir + '/' + sign+'_y_pred_prob_lr.npy', y_pred_prob_lr)
    
    np.save(test_dir + '/' + sign+'_y_test.npy', y_test)
    
    # 保存X_train, X_test, y_train, y_test
    np.save(test_dir + '/' + sign+'_X_train.npy', X_train)
    np.save(test_dir + '/' + sign+'_X_test.npy', X_test)
    np.save(test_dir + '/' + sign+'_y_train.npy', y_train)
    np.save(test_dir + '/' + sign+'_y_test.npy', y_test)
    
    scores.to_csv(test_dir + '/' + sign+'_scores.csv')
    
    return scores
