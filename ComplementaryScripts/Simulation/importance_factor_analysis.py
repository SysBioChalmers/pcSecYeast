#!/usr/bin/env python
# -*- coding: utf-8 -*-

import shap
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from matplotlib import rc
from matplotlib import pyplot as plt


def ABS_SHAP(df_shap,df):
    #import matplotlib as plt
    # Make a copy of the input data
    shap_v = pd.DataFrame(df_shap)
    feature_list = df.columns
    shap_v.columns = feature_list
    df_v = df.copy().reset_index().drop('index',axis=1)
    
    # Determine the correlation in order to plot with different colors
    corr_list = list()
    for i in feature_list:
        b = np.corrcoef(shap_v[i],df_v[i])[1][0]
        corr_list.append(b)
    corr_df = pd.concat([pd.Series(feature_list),pd.Series(corr_list)],axis=1).fillna(0)
    # Make a data frame. Column 1 is the feature, and Column 2 is the correlation coefficient
    corr_df.columns  = ['Variable','Corr']
    corr_df['Sign'] = np.where(corr_df['Corr']>0,'red','blue')
    
    # Plot it
    shap_abs = np.abs(shap_v)
    k=pd.DataFrame(shap_abs.mean()).reset_index()
    k.columns = ['Variable','SHAP_abs']
    k2 = k.merge(corr_df,left_on = 'Variable',right_on='Variable',how='inner')
    k2 = k2.sort_values(by='SHAP_abs',ascending = True)
    colorlist = k2['Sign']
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42

    # plt.axes([0.12,0.12,0.83,0.83])
    
    # plt.tick_params(direction='in')
    # plt.tick_params(which='major',length=1.5)
    # plt.tick_params(which='major',width=0.4)
    # plt.tick_params(which='major',width=0.4)

    # ax = k2.plot.barh(x='Variable',y='SHAP_abs',color = colorlist, figsize=(1.5,1.5),legend=False)

    # ax.spines['bottom'].set_linewidth(0.5)
    # ax.spines['left'].set_linewidth(0.5)
    # ax.spines['top'].set_linewidth(0.5)
    # ax.spines['right'].set_linewidth(0.5)

    # plt.xticks(fontsize=7)
    # plt.yticks(fontsize=6)

    ax = k2.plot.barh(x='Variable',y='SHAP_abs',color = colorlist, figsize=(5,6),legend=False)
    ax.set_xlabel("SHAP Value (Red = Positive Impact)")
    
def main() :
    with open('./FakeTPresMore.txt', 'r') as infile :  # Open the file from a specific directory
        lines = infile.readlines()

    first_line = lines[0].strip().split('\t')
    # print(first_line)

    aas = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']  # 20 amino acids

    feature_names = first_line[2:-2] + aas
    # print(feature_names)
    # print(len(feature_names))  # 24

    cols = feature_names + [first_line[-1]]
    # print(cols)

    i = 0
    all_line = list()
    for line in lines[1:] :
        i += 1
        data = line.strip().split('\t')
        # print(data)

        sequence = data[-2]
        # print(list(sequence))
        s = list(sequence)

        compositions = list()
        for aa in aas :
            compositions.append(float(s.count(aa))/float(len(s)))
            # compositions.append(float(s.count(aa)))
        # print(compositions)
        PTMs = [float(d) for d in data[2:-2]]
        # print(PTMs)

        one_line = PTMs + compositions + [float(data[-1])]
        all_line.append(one_line)

        # if PTMs[:-1] == [0,0,0] :
        #     one_line = PTMs + compositions + [float(data[-1])]
        #     all_line.append(one_line)

    all_data = pd.DataFrame(all_line, columns=cols)

    # print(len(all_data))  # 374 entries
    # print(all_data.maxTP)

    feature_data = pd.DataFrame(all_data,columns=feature_names)
    # print(feature_data)
    # print(type(feature_data))

    labels = all_data.maxTP
    # print(type(labels))

    X_train, X_test, y_train, y_test = train_test_split(feature_data, labels, test_size=0.20, random_state=42)
    # print(len(X_train))
    # print(len(X_test))
    print(X_test)

    rf = RandomForestRegressor(n_estimators=10)
    rf.fit(X_train, y_train)

    # # f = plt.figure(figsize=(1.5, 1.5))

    # # rc('font',**{'family':'serif','serif':['Helvetica']})
    # # plt.rcParams['pdf.fonttype'] = 42

    # # plt.axes([0.12,0.12,0.83,0.83])

    # plt.tick_params(direction='in')
    # plt.tick_params(which='major',length=1.5)
    # plt.tick_params(which='major',width=0.4)
    # plt.tick_params(which='major',width=0.4)

    # # Method 1:
    # # print(rf.feature_importances_)
    # sorted_idx = rf.feature_importances_.argsort()
    # # print(sorted_idx)
    # print(np.array(feature_names)[sorted_idx])
    # print(rf.feature_importances_[sorted_idx])

    # plt.barh(np.array(feature_names)[sorted_idx], rf.feature_importances_[sorted_idx])
    # plt.xlabel('Feature importance')
    # plt.savefig('./feature_importance_analysis_rf.pdf', dpi=400)

    # Method 2:  # conda install -c conda-forge shap
    # https://mljar.com/blog/feature-importance-in-random-forest/
    # https://mljar.com/blog/feature-importance-in-random-forest/
    # https://towardsdatascience.com/explain-your-model-with-the-shap-values-bc36aac4de3d
    # https://towardsdatascience.com/a-novel-approach-to-feature-importance-shapley-additive-explanations-d18af30fc21b
    # https://shap.readthedocs.io/en/latest/example_notebooks/tabular_examples/tree_based_models/Explaining%20the%20Loss%20of%20a%20Model.html

    # The feature importance can be plotted with more details, showing the feature value:
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X_train)
    # print(shap_values)
    ABS_SHAP(shap_values,X_train) 

    # shap_values = explainer.shap_values(X_test)
    # ABS_SHAP(shap_values,X_test) 

    # shap.summary_plot(shap_values, X_train, show=False)
    # shap.summary_plot(shap_values, X_train, plot_type="bar")

    # plt.xticks(fontsize=7)
    # plt.yticks(fontsize=6)

    # plt.savefig('./feature_importance_analysis_aa_value.pdf', dpi=400)
    plt.savefig('./feature_importance_analysis_new.pdf', dpi=400)


if __name__== "__main__":
    main()
