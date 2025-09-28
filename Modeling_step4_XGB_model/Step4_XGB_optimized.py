# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 11:41:27 2025

@author: PC
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import GroupShuffleSplit, StratifiedKFold
from xgboost import XGBClassifier
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score)
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import matplotlib.pyplot as plt
import joblib
import os

os.chdir("E:\\Thyroid_complex_XGB\\")
output_dir = "xgb_model_best_auprc_202506012"
os.makedirs(output_dir, exist_ok=True)

# 1. Process
data = pd.read_csv("TPC_912135pairs_featureswithlabel_20250608.tsv", sep="\t")
data['group_id'] = data['perturbation_A'].str.split('_').str[0] # ID have already been sorted when generating pair

# 2. Group the ID (original+perturbation) - train 80%, valid 10%, test 10%
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
train_idx, temp_idx = next(gss.split(data, groups=data['group_id']))
train_data, temp_data = data.iloc[train_idx], data.iloc[temp_idx]
gss_temp = GroupShuffleSplit(n_splits=1, test_size=0.5, random_state=42)
valid_idx, test_idx = next(gss_temp.split(temp_data, groups=temp_data['group_id']))
valid_data, test_data = temp_data.iloc[valid_idx], temp_data.iloc[test_idx]

# 3. Feature generation
features = ['pearson_R_smoothed', 'euclidean_distance', 'co_peak', 'pearson_R_raw','WCC_score'] # wih WCC
X_train, y_train = train_data[features], train_data['label']
X_valid, y_valid = valid_data[features], valid_data['label']
X_test, y_test = test_data[features], test_data['label']
scale_pos_weight = len(y_train[y_train==0]) / len(y_train[y_train==1]) # calculate the feature weight

# def xgb_auprc_score(y_true, y_pred):
#     return 'auprc', average_precision_score(y_true, y_pred)

# 4. Optimized para
def train_optimized_xgb(X_train, y_train, X_valid, y_valid, n_folds=3):
    best_auprc = 0.5 # need test initial might = 0
    # best_model = None
    models = []
    # auprc_scores = []    
    # setting para
    param_grid = {
        'max_depth': [4, 6, 8],
        'learning_rate': [0.01, 0.05, 0.1],
        'subsample': [0.6, 0.8, 1.0],
        'colsample_bytree': [0.6, 0.8, 1.0],
        'min_child_weight': [1, 3, 5],
        'gamma': [0, 0.1, 0.2]}    
    # grid search
    for max_depth in param_grid['max_depth']:
        for lr in param_grid['learning_rate']:
            for subsample in param_grid['subsample']:
                for colsample in param_grid['colsample_bytree']:
                    for min_child in param_grid['min_child_weight']:
                        for gamma in param_grid['gamma']:
                            params = {
                                'max_depth': max_depth,
                                'learning_rate': lr,
                                'subsample': subsample,
                                'colsample_bytree': colsample,
                                'min_child_weight': min_child,
                                'gamma': gamma,
                                'scale_pos_weight': scale_pos_weight,
                                'objective': 'binary:logistic',
                                'eval_metric': 'aucpr',
                                'random_state': 42,
                                'n_jobs': -1}                            
                            print(f"\nTraining with params: {params}") # print the processing                            
                            # use cross validation
                            kfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42) # n_fold default is 3
                            fold_auprc = []                            
                            for train_idx, val_idx in kfold.split(X_train, y_train):
                                X_tr, X_val = X_train.iloc[train_idx], X_train.iloc[val_idx]
                                y_tr, y_val = y_train.iloc[train_idx], y_train.iloc[val_idx]                                
                                model = XGBClassifier(**params)
                                model.fit(
                                    X_tr, y_tr,
                                    eval_set=[(X_val, y_val)],
                                    early_stopping_rounds=20,
                                    verbose=False)                                
                                y_pred = model.predict_proba(X_val)[:, 1]
                                auprc = average_precision_score(y_val, y_pred)
                                fold_auprc.append(auprc)
                                models.append((model, auprc))                            
                            mean_auprc = np.mean(fold_auprc)
                            print(f"Mean AUPRC: {mean_auprc:.4f}") # calculate the mean of AUPRC                            
                            if mean_auprc > best_auprc: # here best_auprc = 0.5
                                best_auprc = mean_auprc
                                best_params = params
                                print(f"New best AUPRC: {best_auprc:.4f}") # print the AUPRC>0.5    
    # use the best param in all training dataset
    print("\nTraining final model with best parameters...")
    final_model = XGBClassifier(**best_params)
    final_model.fit(
        X_train, y_train,
        eval_set=[(X_valid, y_valid)],
        early_stopping_rounds=20, # early stopping setting
        verbose=50)    
    # save the best AUPRC model
    models.sort(key=lambda x: x[1], reverse=True)
    best_auprc_model = models[0][0]    
    # Compare the final and best model
    final_model_auprc = average_precision_score(y_valid, final_model.predict_proba(X_valid)[:, 1])
    best_auprc_model_auprc = average_precision_score(y_valid, best_auprc_model.predict_proba(X_valid)[:, 1])    
    if best_auprc_model_auprc > final_model_auprc:
        print("Using best AUPRC model from cross-validation") # if refers to the grid search is working
        final_model = best_auprc_model    
    return final_model

# 6. Train the model
optimized_model = train_optimized_xgb(X_train, y_train, X_valid, y_valid)

# 7. Evaluate the model
def evaluate_model(model, X, y, set_name):
    y_pred = model.predict(X)
    y_prob = model.predict_proba(X)[:, 1]    
    metrics = {
        'Accuracy': accuracy_score(y, y_pred),
        'Precision': precision_score(y, y_pred),
        'Recall': recall_score(y, y_pred),
        'F1': f1_score(y, y_pred),
        'AUC-ROC': roc_auc_score(y, y_prob),
        'AUPRC': average_precision_score(y, y_prob)}    
    print(f"\n{set_name} Metrics:")
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")    
    importance = model.get_booster().get_score(importance_type='gain') # feature importance
    importance_df = pd.DataFrame({
        'Feature': features,
        'Importance': [importance.get(f, 0) for f in features]}).sort_values('Importance', ascending=False)
    plot_curves(y, y_prob, set_name)    
    return metrics, importance_df

# 8. Draw the curve
def plot_curves(y_true, y_prob, set_name):
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = auc(fpr, tpr) # ROC    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve ({set_name})')
    plt.legend(loc="lower right")
    plt.savefig(f"{output_dir}/roc_curve_{set_name}.png")
    plt.close()
    
    precision, recall, _ = precision_recall_curve(y_true, y_prob) # PR
    auprc = average_precision_score(y_true, y_prob)    
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, color='blue', lw=2, label=f'PR curve (AUPRC = {auprc:.2f})')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve ({set_name})')
    plt.legend(loc="upper right")
    plt.savefig(f"{output_dir}/pr_curve_{set_name}.png")
    plt.close()

# 9. Save results
print("\n== Final Evaluation ==")
train_metrics, train_importance = evaluate_model(optimized_model, X_train, y_train, "Training")
valid_metrics, valid_importance = evaluate_model(optimized_model, X_valid, y_valid, "Validation")
test_metrics, test_importance = evaluate_model(optimized_model, X_test, y_test, "Test")
# model
joblib.dump(optimized_model, f"{output_dir}/optimized_xgb_model_auprc.joblib")
# results - auc recall ...
metrics_df = pd.DataFrame({
    'Dataset': ['Training', 'Validation', 'Test'],
    'AUPRC': [train_metrics['AUPRC'], valid_metrics['AUPRC'], test_metrics['AUPRC']],
    'AUC-ROC': [train_metrics['AUC-ROC'], valid_metrics['AUC-ROC'], test_metrics['AUC-ROC']],
    'Recall': [train_metrics['Recall'], valid_metrics['Recall'], test_metrics['Recall']]})
    'Dataset': ['Training', 'Validation', 'Test'],
    'AUPRC': [train_metrics['AUPRC'], valid_metrics['AUPRC'], test_metrics['AUPRC']],
    'AUC-ROC': [train_metrics['AUC-ROC'], valid_metrics['AUC-ROC'], test_metrics['AUC-ROC']],
    'Recall': [train_metrics['Recall'], valid_metrics['Recall'], test_metrics['Recall']]
})
metrics_df.to_csv(f"{output_dir}/final_results_metrix.csv", index=False)
# feature importance
test_importance.to_csv(f"{output_dir}/test_feature_importance.csv", index=False)
# finish summary
print("\n== Complete ==")
print(f"Best model saved with AUPRC: {test_metrics['AUPRC']:.4f} on test set")
print(f"Results saved to: {output_dir}")