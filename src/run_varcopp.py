"""Functions to predict feature vectors with the VarCoPP2.0 predictor"""
import numpy as np


def predict(run_d, test_d, varcopp_model):
    results = {}
    prob_classes = varcopp_model.predict_proba(run_d)
    prob_classes = np.array(prob_classes)
    prob_classes = prob_classes[:, 1]
    for i in range(len(test_d)):
        comb_id = test_d[i]
        prob = prob_classes[i]
        if prob > 0.5:
            results[comb_id] = [prob, 'Disease-causing']
        else:
            results[comb_id] = [prob, 'Neutral']
    return results


def run_varcopp(feat_values, comb_ids, varcopp):
    """Runs VarCoPP predictor in chunks"""
    results = {}
    test_d = []
    run_d = []
    count = 0
    for comb_id, comb_feat in zip(comb_ids, feat_values):
        if len(run_d) <= 1000:
            test_d.append(comb_id)
            run_d.append(comb_feat)
        else:
            test_d.append(comb_id)
            run_d.append(comb_feat)
            count += 1
            results.update(predict(run_d, test_d, varcopp))
            test_d = []
            run_d = []
    if len(run_d) == 1:
        # If only one prediction is left array needs to be reshaped
        results.update(predict(np.array(run_d).reshape(1, -1), test_d, varcopp))
    elif len(run_d) > 1:
        results.update(predict(run_d, test_d, varcopp))
    return results