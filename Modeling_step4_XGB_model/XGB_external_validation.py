# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 12:59:58 2025

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
os.chdir("E:\\XGB_model\\")
output_dir = "xgb_maxauprc_TPC_valiFTC133"
model_path = os.path.join(output_dir, 'optimized_xgb_model_auprc.joblib')
ftc_file = "FTC133_vali_data_20250608.tsv"  # external validation data with label and features
features = ['pearson_R_smoothed', 'euclidean_distance', 'co_peak', 'pearson_R_raw', 'WCC_score']
eval_dir = os.path.join(output_dir, "external_validation_results")
os.makedirs(eval_dir, exist_ok=True)

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
        y_prob = model.predict_proba(X_clear)[:, 1]
        # calculate the auc,f1 ...
        metrics = {
            'n_samples': len(X_clear),
            'accuracy': accuracy_score(y_clear, y_pred),
            'precision': precision_score(y_clear, y_pred),
            'recall': recall_score(y_clear, y_pred),
            'f1': f1_score(y_clear, y_pred),
            'auc_roc': roc_auc_score(y_clear, y_prob),
            'auprc': average_precision_score(y_clear, y_prob)}
        pd.DataFrame([metrics]).to_csv(os.path.join(output_dir, "metrics_clear_labels.csv"), index=False) # save results

        # ROC Curve
        fpr, tpr, _ = roc_curve(y_clear, y_prob)
        plt.figure(figsize=(8, 6))
        plt.plot(fpr, tpr, label=f'AUC = {metrics["auc_roc"]:.3f}')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve (L1-adjusted Model)')
        plt.legend()
        plt.savefig(os.path.join(output_dir, "roc_curve.png"))
        plt.close()

        # PR Curve
        precision, recall, _ = precision_recall_curve(y_clear, y_prob)
        plt.figure(figsize=(8, 6))
        plt.plot(recall, precision, label=f'AUPRC = {metrics["auprc"]:.3f}')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve (L1-adjusted Model)')
        plt.legend()
        plt.savefig(os.path.join(output_dir, "pr_curve.png"))
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
        plt.savefig(os.path.join(output_dir, "calibration_curve.png"))
        plt.close()

        print("\nLabel results:")
        print(pd.Series(metrics))

 # 2.2 predict the label (NaN)
    if len(ambiguous_data) > 0:
        X_amb = ambiguous_data[features]
        y_prob_amb = model.predict_proba(X_amb)[:, 1]
        amb_results = ambiguous_data.copy()
        amb_results['predicted_prob'] = y_prob_amb
        amb_results['predicted_label'] = (y_prob_amb > 0.5).astype(int)
        amb_results.to_csv(os.path.join(output_dir, "ambiguous_predictions.tsv"), sep="\t", index=False)
        # prediction prob histogram
        plt.figure(figsize=(10, 6))
        plt.hist(y_prob_amb, bins=50, color='purple', alpha=0.7)
        plt.xlabel('Predicted Probability')
        plt.ylabel('Count')
        plt.title('Prediction Distribution for Ambiguous Samples (L1 Model)')
        plt.savefig(os.path.join(output_dir, "ambiguous_distribution.png"))
        plt.close()
        # check the high confidence
        high_confidence_mask = y_prob_amb >= 0.9
        high_confidence_pairs = amb_results.loc[high_confidence_mask]
        high_confidence_pairs.to_csv(os.path.join(output_dir, "high_confidence_predictions.tsv"), sep="\t", index=False)
        print(f"\nThe high confidence (>0.9): {len(high_confidence_pairs)}")
    # optimized para（maximum F1 score）
    if len(clear_data) > 0:
        precision_vals, recall_vals, thresholds = precision_recall_curve(y_clear, y_prob)
        f1_scores = 2 * (precision_vals * recall_vals) / (precision_vals + recall_vals + 1e-8)
        best_idx = np.argmax(f1_scores)
        best_threshold = thresholds[best_idx]
        print(f"\n Best cut-off（F1 score）: {best_threshold:.4f} | F1: {f1_scores[best_idx]:.4f}")
        # use the threshold - optimized
        y_pred_opt = (y_prob >= best_threshold).astype(int)
        opt_metrics = {
            'threshold': best_threshold,
            'precision_opt': precision_score(y_clear, y_pred_opt),
            'recall_opt': recall_score(y_clear, y_pred_opt),
            'f1_opt': f1_score(y_clear, y_pred_opt)
        }
        print("\nOptimized results:")
        print(opt_metrics)
    # Precision of predicted label - draw figure
    if len(ambiguous_data) > 0:
        X_amb = ambiguous_data[features]
        y_prob_amb = model.predict_proba(X_amb)[:, 1]
        prob_sorted = np.sort(y_prob_amb)[::-1]
        cumulative_precision = np.cumsum(prob_sorted) / (np.arange(len(prob_sorted)) + 1) # use high prob replace the precision
        # Draw precision plot
        plt.figure(figsize=(12, 6))
        plt.plot(np.arange(1, len(prob_sorted)+1), cumulative_precision, 
             color='darkorange', lw=2)    
        plt.xlabel('Number of Protein Pairs (Sorted by Confidence)', fontsize=12)
        plt.ylabel('Estimated Precision', fontsize=12)
        plt.title('Precision at Top-K Predictions (FTC133 Unlabeled Data)', fontsize=14)    
    # Key point labeling
    for k in [5000, 20000, 50000]:
        if k < len(prob_sorted):
            plt.scatter(k, cumulative_precision[k-1], color='red', zorder=5)
            plt.text(k, cumulative_precision[k-1]+0.02, 
                    f'Top {k}: {cumulative_precision[k-1]:.2%}',
                    ha='center', fontsize=10)    
    plt.grid(True, alpha=0.3)
    plt.xlim(0, len(prob_sorted)) # limited 1k points in the fig - plt.xlim(0, min(100000, len(prob_sorted))) 
    plt.ylim(0, 1.05)  
    
    plot_path = os.path.join(output_dir, "precision_unlabeled_allpoint.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nSave Precision-protein pairs in: {plot_path}")
    
    # save the sorted pairs
    sorted_predictions = ambiguous_data.copy()
    sorted_predictions['predicted_prob'] = y_prob_amb
    sorted_predictions = sorted_predictions.sort_values('predicted_prob', ascending=False)
    sorted_predictions.to_csv(os.path.join(output_dir, "sorted_predictions_unlabeled.tsv"), 
                             sep="\t", index=False)

# 3. Main function
if __name__ == "__main__":
    model = load_specific_model(model_path)
    evaluate_ftc133(model, ftc_file, features, eval_dir)