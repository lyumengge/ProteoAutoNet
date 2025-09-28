# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 00:16:29 2025

@author: PC
"""


import os
import pandas as pd
import numpy as np
import joblib
from sklearn.calibration import calibration_curve
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, average_precision_score, precision_recall_curve, roc_curve, brier_score_loss)
import matplotlib.pyplot as plt

# 1.Load path and model
os.chdir("E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\xgb_model_best_auprc_20250615")
output_dir = "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\xgb_model_best_auprc_20250615\\"
model_path = os.path.join(output_dir, 'optimized_xgb_model_auprc.joblib')
ftc_file = "E:\\MenggeLYU\\NC_rebuttal_t3\\New_feature_RF\\ComplexXGB\\TPC_all_feature_20250627.tsv"  # external validation data with label and features
features = ['pearson_R_smoothed', 'euclidean_distance', 'co_peak', 'pearson_R_raw', 'WCC_score']
eval_dir = os.path.join(output_dir, "TPC_all_label_withprecision_weighted3")
os.makedirs(eval_dir, exist_ok=True)

# 修改后的加权累积精度计算函数
def calculate_weighted_cumulative_precision(y_true, y_prob):
    """
    计算加权累积精度（NaN样本用预测概率作为权重）
    参数:
        y_true: 包含真实标签的数组（1=正例，0=负例，NaN=未知）
        y_prob: 预测概率数组
    返回:
        累积精度数组
    """
    # 替换NaN为预测概率（表示它是阳性的概率）
    y_true_weighted = np.where(np.isnan(y_true), y_prob, y_true)
    
    # 按预测概率从高到低排序
    sorted_indices = np.argsort(y_prob)[::-1]
    y_true_sorted = y_true_weighted[sorted_indices]
    
    # 计算累积加权TP和样本数
    weighted_TP = np.cumsum(y_true_sorted)  # 加权真阳性
    total_samples = np.arange(1, len(y_true_sorted)+1)
    
    precision = weighted_TP / total_samples
    return precision

# 保留原有的简单精度计算函数用于比较
def calculate_naive_cumulative_precision(labels):
    """
    简单累积精度计算（NaN视为负例）
    """
    labels = np.where(np.isnan(labels), 0, labels)
    TPs = np.cumsum(labels)
    FPs = np.cumsum(1 - labels)
    precision = TPs / (TPs + FPs)
    return precision

# load model
def load_specific_model(model_path):
    if os.path.exists(model_path.replace("_final_", "_v2_")): # make sure use L1 if have L1
        model = joblib.load(model_path.replace("_final_", "_v2_"))
        print("Load L1 model (model_v2)")
    else:
        model = joblib.load(model_path)
        print("Load basic model (model_v1)") # use base model
    return model

# 2. Evaluate external validation datasets
def evaluate_ftc133(model, validation_file, features, output_dir):
    ind_data = pd.read_csv(validation_file, sep="\t") # load validation data
    print(f"\nLabel distribution in external vali:\n{ind_data['label'].value_counts(dropna=False)}")
    # split the label (1,0 and NaN)
    clear_data = ind_data[ind_data['label'].notna()].copy()
    ambiguous_data = ind_data[ind_data['label'].isna()].copy()
    
    # 2.1 evaluate the label (1 and 0)
    if len(clear_data) > 0:
        X_clear = clear_data[features]
        y_clear = clear_data['label'].astype(int) # check the label - int type
        y_pred = model.predict(X_clear)
        y_prob = clear_data['predicted_prob'] = model.predict_proba(X_clear)[:, 1]  # 添加预测概率到数据框
        
        # 保存所有有标签样本的完整结果
        clear_data['predicted_label'] = y_pred  # 添加预测标签
        clear_data.to_csv(os.path.join(output_dir, "all_labeled_pairs_results.tsv"), sep="\t", index=False)
        print(f"\nSaved all labeled pairs results ({len(clear_data)} samples) to all_labeled_pairs_results.tsv")
        
        # 分析有标签样本中预测概率>0.95的pairs
        high_prob_labeled = clear_data[clear_data['predicted_prob'] > 0.95]
        print(f"\nLabeled samples with predicted_prob > 0.95: {len(high_prob_labeled)}")
        
        high_prob_labeled2 = clear_data[clear_data['predicted_prob'] > 0.5]
        print(f"\nLabeled samples with predicted_prob > 0.5: {len(high_prob_labeled2)}")
        
        # 保存高概率有标签样本
        high_prob_labeled.to_csv(os.path.join(output_dir, "high_prob095_labeled_pairs.tsv"), sep="\t", index=False)
        high_prob_labeled2.to_csv(os.path.join(output_dir, "high_prob05_labeled_pairs.tsv"), sep="\t", index=False)
        
        # 原有指标计算
        metrics = {
            'n_samples': len(X_clear),
            'accuracy': accuracy_score(y_clear, y_pred),
            'precision': precision_score(y_clear, y_pred),
            'recall': recall_score(y_clear, y_pred),
            'f1': f1_score(y_clear, y_pred),
            'auc_roc': roc_auc_score(y_clear, y_prob),
            'auprc': average_precision_score(y_clear, y_prob),
            'n_high_prob(>0.95)': len(high_prob_labeled)
        }
        pd.DataFrame([metrics]).to_csv(os.path.join(output_dir, "metrics_clear_labels.csv"), index=False)

        # ROC Curve
        fpr, tpr, _ = roc_curve(y_clear, y_prob)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f'AUC = {metrics["auc_roc"]:.3f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve (L1-adjusted Model)')
        plt.legend()
        plt.savefig(os.path.join(output_dir, "roc_curve.pdf"))
        plt.close()

        # PR Curve
        precision, recall, _ = precision_recall_curve(y_clear, y_prob)
        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, label=f'AUPRC = {metrics["auprc"]:.3f}')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve (L1-adjusted Model)')
        plt.legend()
        plt.savefig(os.path.join(output_dir, "pr_curve.pdf"))
        plt.close()

        # Calibration Curve
        prob_true, prob_pred = calibration_curve(y_clear, y_prob, n_bins=10)
        brier = brier_score_loss(y_clear, y_prob)
        ece = np.sum(np.abs(prob_true - prob_pred) * len(prob_true)) / len(y_clear)
        plt.figure(figsize=(8, 6))
        plt.plot(prob_pred, prob_true, marker='o', label=f'ECE={ece:.4f}, Brier={brier:.4f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('Mean Predicted Probability')
        plt.ylabel('Fraction of Positives')
        plt.title('Calibration Curve')
        plt.legend()
        plt.savefig(os.path.join(output_dir, "calibration_curve.pdf"))
        plt.close()

        print("\nLabel results:")
        print(pd.Series(metrics))

    # 2.2 predict the label (NaN)
    if len(ambiguous_data) > 0:
        X_amb = ambiguous_data[features]
        y_prob_amb = ambiguous_data['predicted_prob'] = model.predict_proba(X_amb)[:, 1]
        
        # 分析无标签样本中预测概率>0.95的pairs
        high_prob_unlabeled = ambiguous_data[ambiguous_data['predicted_prob'] > 0.95]
        print(f"\nUnlabeled samples with predicted_prob > 0.95: {len(high_prob_unlabeled)}")
        
        high_prob_unlabeled2 = ambiguous_data[ambiguous_data['predicted_prob'] > 0.5]
        print(f"\nUnlabeled samples with predicted_prob > 0.5: {len(high_prob_unlabeled2)}")
        
        # 保存高概率无标签样本
        high_prob_unlabeled.to_csv(os.path.join(output_dir, "high_prob095_unlabeled_pairs.tsv"), sep="\t", index=False)
        high_prob_unlabeled2.to_csv(os.path.join(output_dir, "high_prob05_unlabeled_pairs.tsv"), sep="\t", index=False)
        
        # 原有处理
        amb_results = ambiguous_data.copy()
        amb_results['predicted_label'] = (y_prob_amb > 0.5).astype(int)
        amb_results.to_csv(os.path.join(output_dir, "ambiguous_predictions.tsv"), sep="\t", index=False)
        
        # 预测概率直方图
        plt.figure(figsize=(10, 6))
        plt.hist(y_prob_amb, bins=50, color='purple', alpha=0.7)
        plt.xlabel('Predicted Probability')
        plt.ylabel('Count')
        plt.title('Prediction Distribution for Ambiguous Samples (L1 Model)')
        plt.savefig(os.path.join(output_dir, "ambiguous_distribution.pdf"))
        plt.close()
        
        # 高置信度样本
        high_confidence_mask = y_prob_amb >= 0.9
        high_confidence_pairs = amb_results.loc[high_confidence_mask]
        high_confidence_pairs.to_csv(os.path.join(output_dir, "high_confidence_predictions.tsv"), sep="\t", index=False)
        print(f"\nThe high confidence (>0.9): {len(high_confidence_pairs)}")

    if len(clear_data) > 0 or len(ambiguous_data) > 0:
        all_data = pd.concat([clear_data, ambiguous_data])
        all_data = all_data.sort_values('predicted_prob', ascending=False)

        # 获取真实标签和预测概率
        true_labels = all_data['label'].values
        predicted_probs = all_data['predicted_prob'].values
    
        # 计算三种精度
        weighted_precision = calculate_weighted_cumulative_precision(true_labels, predicted_probs)
        naive_precision = calculate_naive_cumulative_precision(true_labels)
    
        # 创建一个 DataFrame 来保存精度信息
        precision_df = pd.DataFrame({
            'protein_pair_rank': range(1, len(all_data)+1),
            'predicted_prob': predicted_probs,
            'weighted_precision': weighted_precision,
            'naive_precision': naive_precision
        })
    
        # 添加 labeled only 精度
        labeled_mask = ~np.isnan(true_labels)
        labeled_true = true_labels[labeled_mask]
    
        if len(labeled_true) > 0:
            labeled_precision_values = calculate_naive_cumulative_precision(labeled_true)
            # 创建与完整数据集相同长度的数组
            labeled_precision_full = np.full(len(all_data), np.nan)
            labeled_indices = np.where(labeled_mask)[0]
            copy_length = min(len(labeled_precision_values), len(labeled_indices))
            labeled_precision_full[labeled_indices[:copy_length]] = labeled_precision_values[:copy_length]
        else:
            labeled_precision_full = np.full(len(all_data), np.nan)
    
        # 将 labeled only precision 添加进 df
        precision_df['labeled_only_precision'] = labeled_precision_full
    
        # ✅ 关键修改：将 all_data 原始数据合并进来
        # 注意：因为 all_data 已经被排序过了，所以直接拼接
        precision_df = pd.concat([all_data.reset_index(drop=True), precision_df.reset_index(drop=True)], axis=1)
    
        # 保存结果
        precision_df.to_csv(os.path.join(eval_dir, "cumulative_precision_all_samples_with_info.tsv"), sep="\t", index=False)
        
        
        # 计算precision>0.95和>0.5的蛋白对数目（使用加权精度）
        n_high_precision_95 = np.sum(weighted_precision > 0.95)
        n_high_precision_50 = np.sum(weighted_precision > 0.5)
        
        print(f"\nNumber of protein pairs with weighted precision > 0.95: {n_high_precision_95}")
        print(f"Number of protein pairs with weighted precision > 0.5: {n_high_precision_50}")
        
        # 将结果保存到metrics文件
        if len(clear_data) > 0:
            metrics_df = pd.read_csv(os.path.join(output_dir, "metrics_clear_labels.csv"))
            metrics_df['n_high_precision(>0.95)'] = n_high_precision_95
            metrics_df['n_high_precision(>0.5)'] = n_high_precision_50
            metrics_df.to_csv(os.path.join(output_dir, "metrics_clear_labels.csv"), index=False)
        
        # 绘制三种累积精度曲线对比图
        plt.figure(figsize=(12, 6))
        plt.plot(precision_df['protein_pair_rank'], precision_df['weighted_precision'], 
                color='blue', label='Weighted (NaN=predicted_prob)')
        plt.plot(precision_df['protein_pair_rank'], precision_df['naive_precision'], 
                color='red', linestyle='--', label='Naive (NaN=0)')
        plt.plot(precision_df['protein_pair_rank'], precision_df['labeled_only_precision'], 
                color='green', linestyle=':', label='Labeled Only')
        
        plt.axhline(y=0.95, color='k', linestyle='--', alpha=0.3)
        plt.axhline(y=0.5, color='k', linestyle='--', alpha=0.3)
        
        plt.xlabel('Protein Pairs (Sorted by Confidence)')
        plt.ylabel('Cumulative Precision')
        plt.title('Comparison of Cumulative Precision Calculation Methods')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(output_dir, "cumulative_precision_comparison.pdf"))
        plt.close()
        
        # 打印Top K的精度（加权）
        print("\nWeighted Precision at Top-K:")
        for k in [1000, 5000, 10000, 20000]:
            if k < len(weighted_precision):
                print(f"Top {k}: {weighted_precision[k-1]:.4f}")
    
    # 优化阈值分析（保持不变）
    if len(clear_data) > 0:
        precision_vals, recall_vals, thresholds = precision_recall_curve(y_clear, y_prob)
        f1_scores = 2 * (precision_vals * recall_vals) / (precision_vals + recall_vals + 1e-8)
        best_idx = np.argmax(f1_scores)
        best_threshold = thresholds[best_idx]
        print(f"\n Best cut-off（F1 score）: {best_threshold:.4f} | F1: {f1_scores[best_idx]:.4f}")
        
        y_pred_opt = (y_prob >= best_threshold).astype(int)
        opt_metrics = {
            'threshold': best_threshold,
            'precision_opt': precision_score(y_clear, y_pred_opt),
            'recall_opt': recall_score(y_clear, y_pred_opt),
            'f1_opt': f1_score(y_clear, y_pred_opt)
        }
        print("\nOptimized results:")
        print(opt_metrics)
    
    # Precision at K分析（保持不变）
    if len(ambiguous_data) > 0:
        prob_sorted = np.sort(y_prob_amb)[::-1]
        cumulative_precision = np.cumsum(prob_sorted) / (np.arange(len(prob_sorted)) + 1)
        
        plt.figure(figsize=(12, 6))
        plt.plot(np.arange(1, len(prob_sorted)+1), cumulative_precision, color='darkorange', lw=2)    
        plt.xlabel('Number of Protein Pairs (Sorted by Confidence)', fontsize=12)
        plt.ylabel('Estimated Precision', fontsize=12)
        plt.title('Precision at Top-K Predictions (FTC133 Unlabeled Data)', fontsize=14)    
        
        for k in [5000, 20000, 50000]:
            if k < len(prob_sorted):
                plt.scatter(k, cumulative_precision[k-1], color='red', zorder=5)
                plt.text(k, cumulative_precision[k-1]+0.02, 
                        f'Top {k}: {cumulative_precision[k-1]:.2%}',
                        ha='center', fontsize=10)    
        
        plt.grid(True, alpha=0.3)
        plt.xlim(0, len(prob_sorted))
        plt.ylim(0, 1.05)  
        
        plot_path = os.path.join(output_dir, "precision_unlabeled_allpoint.pdf")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nSave Precision-protein pairs in: {plot_path}")
        
        sorted_predictions = ambiguous_data.copy()
        sorted_predictions = sorted_predictions.sort_values('predicted_prob', ascending=False)
        sorted_predictions.to_csv(os.path.join(output_dir, "sorted_predictions_unlabeled.tsv"), 
                                 sep="\t", index=False)

# 3. Main function
if __name__ == "__main__":
    model = load_specific_model(model_path)
    evaluate_ftc133(model, ftc_file, features, eval_dir)
    
