# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import tensorflow as tf
import matplotlib.pyplot as plt
import seaborn as sns
import six
import sys
sys.modules['sklearn.externals.six'] = six
from mlrose import NeuralNetwork, GeomDecay
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA, TruncatedSVD
import matplotlib.patches as mpatches
import time
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, accuracy_score, \
    classification_report, make_scorer
from scipy.stats import norm
# Classifier Libraries
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
import collections

from sklearn.preprocessing import StandardScaler, RobustScaler
# Other Libraries
from sklearn.model_selection import train_test_split, GridSearchCV, cross_validate
from sklearn.pipeline import make_pipeline
from imblearn.pipeline import make_pipeline as imbalanced_make_pipeline
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import NearMiss
from imblearn.metrics import classification_report_imbalanced
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, accuracy_score, \
    classification_report
from collections import Counter
from sklearn.model_selection import KFold, StratifiedKFold
import warnings

warnings.filterwarnings("ignore")


def read__credit_card_csv():
    return pd.read_csv('data/creditcard.csv')

from imblearn.over_sampling import RandomOverSampler



def main(name):
    df = read__credit_card_csv()
    print(df.head())
    print(df.describe())
    print(df.isnull().sum().max())
    print(df.columns)

    print('No Frauds', round(df['Class'].value_counts()[0] / len(df) * 100, 2), '% of the dataset')
    print('Frauds', round(df['Class'].value_counts()[1] / len(df) * 100, 2), '% of the dataset')

    # Class Distribution
    colors = ["#0101DF", "#DF0101"]
    #
    sns.countplot('Class', data=df, palette=colors)
    plt.title('Class Distributions \n (0: No Fraud || 1: Fraud)', fontsize=14)

    # Distirbution of amount and time
    fig, ax = plt.subplots(1, 2, figsize=(18, 4))

    amount_val = df['Amount'].values
    time_val = df['Time'].values

    sns.distplot(amount_val, ax=ax[0], color='r')
    ax[0].set_title('Distribution of Transaction Amount', fontsize=14)
    ax[0].set_xlim([min(amount_val), max(amount_val)])

    sns.distplot(time_val, ax=ax[1], color='b')
    ax[1].set_title('Distribution of Transaction Time', fontsize=14)
    ax[1].set_xlim([min(time_val), max(time_val)])
    #
    #
    #
    #
    # plt.show()

    std_scaler = StandardScaler()
    rob_scaler = RobustScaler()

    df['scaled_amount'] = rob_scaler.fit_transform(df['Amount'].values.reshape(-1, 1))
    df['scaled_time'] = rob_scaler.fit_transform(df['Time'].values.reshape(-1, 1))

    df.drop(['Time', 'Amount'], axis=1, inplace=True)

    X = df.drop('Class', axis=1)
    y = df['Class']

    sss = StratifiedKFold(n_splits=5, random_state=None, shuffle=False)

    for train_index, test_index in sss.split(X, y):
        print("Train:", train_index, "Test:", test_index)
        original_Xtrain, original_Xtest = X.iloc[train_index], X.iloc[test_index]
        original_ytrain, original_ytest = y.iloc[train_index], y.iloc[test_index]

    original_Xtrain = original_Xtrain.values
    original_Xtest = original_Xtest.values
    original_ytrain = original_ytrain.values
    original_ytest = original_ytest.values

    # See if both the train and test label distribution are similarly distributed
    train_unique_label, train_counts_label = np.unique(original_ytrain, return_counts=True)
    test_unique_label, test_counts_label = np.unique(original_ytest, return_counts=True)
    print('-' * 100)

    print('Label Distributions: \n')
    print(train_counts_label / len(original_ytrain))
    print(test_counts_label / len(original_ytest))

    df = df.sample(frac=1)

    # amount of fraud classes 492 rows.
    fraud_df = df.loc[df['Class'] == 1]
    non_fraud_df = df.loc[df['Class'] == 0][:len(fraud_df)]

    normal_distributed_df = pd.concat([fraud_df, non_fraud_df])

    # Shuffle dataframe rows
    new_df = normal_distributed_df.sample(frac=1, random_state=42)

    print('Distribution of the Classes in the subsample dataset')
    print(new_df['Class'].value_counts() / len(new_df))

    sns.countplot('Class', data=new_df, palette=colors)
    plt.title('Equally Distributed Classes', fontsize=14)

    # Make sure we use the subsample in our correlation

    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 20))

    # Entire DataFrame
    corr = df.corr()
    sns.heatmap(corr, cmap='coolwarm_r', annot_kws={'size': 20}, ax=ax1)
    ax1.set_title("Imbalanced Correlation Matrix ", fontsize=14)

    sub_sample_corr = new_df.corr()
    sns.heatmap(sub_sample_corr, cmap='coolwarm_r', annot_kws={'size': 20}, ax=ax2)
    ax2.set_title('SubSample Correlation Matrix ', fontsize=14)

    plt.savefig('correl.png')

    # f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(20, 6))
    #
    # v14_fraud_dist = new_df['V14'].loc[new_df['Class'] == 1].values
    # sns.distplot(v14_fraud_dist, ax=ax1, fit=norm, color='#FB8861')
    # ax1.set_title('V14 Distribution \n (Fraud Transactions)', fontsize=14)
    #
    # v12_fraud_dist = new_df['V12'].loc[new_df['Class'] == 1].values
    # sns.distplot(v12_fraud_dist, ax=ax2, fit=norm, color='#56F9BB')
    # ax2.set_title('V12 Distribution \n (Fraud Transactions)', fontsize=14)
    #
    # v10_fraud_dist = new_df['V10'].loc[new_df['Class'] == 1].values
    # sns.distplot(v10_fraud_dist, ax=ax3, fit=norm, color='#C5B3F9')
    # ax3.set_title('V10 Distribution \n (Fraud Transactions)', fontsize=14)
    #
    #
    # v12_fraud_dist = new_df['V11'].loc[new_df['Class'] == 1].values
    # sns.distplot(v12_fraud_dist, ax=ax4, fit=norm, color='#56F9BB')
    # ax4.set_title('V11 Distribution \n (Fraud Transactions)', fontsize=14)
    #
    # v10_fraud_dist = new_df['V4'].loc[new_df['Class'] == 1].values
    # sns.distplot(v10_fraud_dist, ax=ax5, fit=norm, color='#C5B3F9')
    # ax5.set_title('V4 Distribution \n (Fraud Transactions)', fontsize=14)
    #
    # f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(20, 6))

    # v14_fraud_dist = new_df['V14'].loc[new_df['Class'] == 0].values
    # sns.distplot(v14_fraud_dist, ax=ax1, fit=norm, color='#FB8861')
    # ax1.set_title('V14 Distribution \n (Non Fraud Transactions)', fontsize=14)
    #
    # v12_fraud_dist = new_df['V12'].loc[new_df['Class'] == 0].values
    # sns.distplot(v12_fraud_dist, ax=ax2, fit=norm, color='#56F9BB')
    # ax2.set_title('V12 Distribution \n (Non Fraud Transactions)', fontsize=14)
    #
    # v10_fraud_dist = new_df['V10'].loc[new_df['Class'] == 0].values
    # sns.distplot(v10_fraud_dist, ax=ax3, fit=norm, color='#C5B3F9')
    # ax3.set_title('V10 Distribution \n (Non Fraud Transactions)', fontsize=14)
    #
    #
    # v12_fraud_dist = new_df['V11'].loc[new_df['Class'] == 0].values
    # sns.distplot(v12_fraud_dist, ax=ax4, fit=norm, color='#56F9BB')
    # ax4.set_title('V11 Distribution \n (Non Fraud Transactions)', fontsize=14)
    #
    # v10_fraud_dist = new_df['V4'].loc[new_df['Class'] == 0].values
    # sns.distplot(v10_fraud_dist, ax=ax5, fit=norm, color='#C5B3F9')
    # ax5.set_title('V4 Distribution \n (Non Fraud Transactions)', fontsize=14)
    plt.close()

    # outlier(new_df,'outlier')
    #
    new_df = remove_outlier(new_df)
    #
    # outlier(new_df_up, 'after_removing_outlier')
    # decisionTreeGridSearch(new_df)
    # decision_tree_training_sample(new_df)
    # decistiontree_pruning(new_df)
    # decistiontree_min_leaf(new_df)
    # decisionTreeGridSearch(new_df)
    # neural_network_grid_search(new_df)
    # neural_network_size(new_df)
    # neural_network_impact_ofweights(new_df)
    # knn_with_dim_reduction(new_df)
    # knn_test(new_df)
    # knn_test_k_size(new_df)
    # ada_boost_estimators(new_df)
    # ada_boost_size(new_df)
    # svm_kernels(new_df)
    #
    #
    # svn_learning_curve(new_df)
    # X,y = split_data(new_df)
    # models_evaluation(X, y,5)
    # confusion_matrix(new_df)
    test_fit_random_hill_climb(new_df)

    test_fit_random_sa(new_df)

    test_fit_random_ga(new_df)

    # plt.show()

from sklearn import svm, datasets

def svm_kernels(new_df, img_name="svm_kernel"):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    kernels = ['linear', 'sigmoid' , 'rbf', 'poly' ]
    C = 1.0
    for kernel in kernels:
        sizes.append(kernel)

        X, y = split_data(new_df)
        svc_classifier = svm.SVC(kernel=kernel, C=C)

        scr = cross_validate(svc_classifier, X, y, cv=5, return_train_score=True, n_jobs=-1)
        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())

    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format(img_name))
    # plt.show()
def svn_learning_curve(new_df, img_name="svm_learning"):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    C = 1.0
    for size in range(50, 900, 50):
        sizes.append(size)
        sample_size = new_df.sample(n=size)
        X, y = split_data(sample_size)
        # X, y = split_data(new_df)
        svc_classifier = svm.SVC(kernel='linear', C=C)

        scr = cross_validate(svc_classifier, X, y, cv=5, return_train_score=True, n_jobs=-1)

        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())
    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format(img_name))
    # plt.show()


from sklearn.ensemble import AdaBoostClassifier

def ada_boost_size(new_df, img_name="ada_boost_size"):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    for size in range(50, len(new_df), 50):
        sizes.append(size)
        sample_size = new_df.sample(n=size)
        X, y = split_data(sample_size)
        dt = DecisionTreeClassifier(criterion="gini",
                                    splitter="best",
                                    max_depth=3,
                                    min_samples_split=2,
                                    min_samples_leaf=2,
                                    min_weight_fraction_leaf=0.0,
                                    max_features=None,
                                    random_state=None,
                                    max_leaf_nodes=None,
                                    min_impurity_decrease=0.0,
                                    class_weight=None,
                                    ccp_alpha=0.0)
        model = AdaBoostClassifier(base_estimator=dt, n_estimators=100)

        scr = cross_validate(model, X, y, cv=5, return_train_score=True, n_jobs=-1, scoring="recall")
        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())
    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format(img_name))
    # plt.show()

def ada_boost_estimators(new_df):

    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    new_df = new_df.reset_index(drop=True)
    X, Y = split_data(new_df)
    X = scaled_data(X)

    dt = DecisionTreeClassifier(criterion="gini",
                                splitter="best",
                                max_depth=5,
                                min_samples_split=2,
                                min_samples_leaf=2,
                                min_weight_fraction_leaf=0.0,
                                max_features=None,
                                random_state=None,
                                max_leaf_nodes=None,
                                min_impurity_decrease=0.0,
                                class_weight=None,
                                ccp_alpha=0.0)
    model = AdaBoostClassifier(base_estimator=dt, n_estimators=10)
    viz = ValidationCurve(model, param_name="n_estimators", n_jobs=-1,
                          param_range=np.arange(1, 100), cv=5, scoring="recall")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="ada_boost_estimators.png")







from sklearn.neighbors import KNeighborsClassifier

# def knn_test(new_df):
#
#     model = KNeighborsClassifier(n_neighbors=3)

def knn_with_dim_reduction(new_df):
    new_df = pca_dimension_reduction(new_df)
    knn_test(new_df, "with_dim_reduction")

def knn_test_k_size(new_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    new_df = new_df.reset_index(drop=True)
    X, Y = split_data(new_df)
    X = scaled_data(X)
    clf = KNeighborsClassifier(n_neighbors=3)
    # Creating a visualizer object
    viz = ValidationCurve(clf, param_name="n_neighbors", n_jobs=-1,
                          param_range=np.arange(1, 50), cv=5, scoring="precision")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="knn_test_k_size.png")



def knn_test(new_df, img_name="without_dim_reduction"):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    for size in range(50, len(new_df), 50):
        sizes.append(size)
        sample_size = new_df.sample(n=size)
        X, y = split_data(sample_size)
        clf = KNeighborsClassifier(n_neighbors=20)

        scr = cross_validate(clf, X, y, cv=5, return_train_score=True, n_jobs=-1)
        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())
    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format(img_name))
    # plt.show()


def scaled_data(X):
    scaler = StandardScaler().fit(X)
    X_train_scaled_transformed = pd.DataFrame(
        scaler.transform(X), columns=[
            c for c in X.columns if "target" not in c
        ]
    )
    return X_train_scaled_transformed


def pca_dimension_reduction(new_df):
    # we apply a standard scaler:
    X, y = split_data(new_df)
    scaler = StandardScaler().fit(X)
    X_train_scaled_transformed = pd.DataFrame(
        scaler.transform(X), columns=[
            c for c in X.columns if "target" not in c
        ]
    )

    pca = PCA()
    pca.fit(X_train_scaled_transformed)
    comps = len([c for c in list(pca.explained_variance_ratio_) if c >= 0.05])
    pca_train = PCA(n_components=comps)

    df_pca_train = pd.DataFrame(pca_train.fit_transform(
        X_train_scaled_transformed), columns=[
        f"{co}_COMP" for co in range(0, comps)
    ]
    )
    y_train = pd.DataFrame(y, columns=["Class"]).reset_index(drop=True)


    df_pca_for_knn = pd.concat(
        [
            df_pca_train, y_train
        ], axis=1
    )
    return df_pca_for_knn



from sklearn.neural_network import MLPClassifier



def test_fit_random_ga(new_df):
    """Test fit method using the random hill climbing algorithm"""
    X, Y = split_data(new_df)
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    # X, X_test, y, y_test = train_test_split(X, y, test_size=0.33)
    # ensure all data are floating point values
    X = X.to_numpy()
    X = X.astype('float32')
    X = X / 255
    Y = Y.to_numpy()


    traiing_scores_cool = list()
    testing_score_cool = list()

    mutations_prob = [0.0,0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
    restarts = [10, 100, 250, 500, 1000, 2000]
    X_train, X_test, y_train, y_test = train_test_split(X, Y, random_state = 0)
    for mp in mutations_prob:
        traiing_scores = list()
        testing_score = list()
        for iter in restarts:
            network = NeuralNetwork(hidden_nodes=[7,7], activation='relu',
                                    algorithm='genetic_alg',
                                    bias=False, is_classifier=True,
                                    learning_rate=1, clip_max=1,
                                    mutation_prob=mp,
                                    max_iters=iter,
                                    max_attempts=100)



            network.fit(X_train, y_train)
            y_train_prd = network.predict(X_train)

            y_test_pred = network.predict(X_test)




            testing_score.append(accuracy_score(y_train, y_train_prd))
            traiing_scores.append(accuracy_score(y_test, y_test_pred))
            print("Random Restart - {}".format(mp))

        traiing_scores_cool.append(traiing_scores)
        testing_score_cool.append(testing_score)
    #
    # plt_df = pd.DataFrame(
    #     {'x_values': restarts, 'Training score': traiing_scores_cool})
    for training_error in traiing_scores_cool:
        plt_df = pd.DataFrame(
            {'x_values': restarts, 'Training score': training_error})
        plt.plot('x_values', 'Training score', data=plt_df, marker='o',  markersize=12,
                 linewidth=4)
    plt.legend(mutations_prob, title='Mutation Probability')
    plt.savefig('{}.png'.format('ga_neural_network_size_training'))

    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    for testing_score in testing_score_cool:
        plt_df = pd.DataFrame(
            {'x_values': restarts, 'Testing score': testing_score})
        plt.plot('x_values', 'Testing score', data=plt_df, marker='o',  markersize=12,
                 linewidth=4)
    plt.legend(mutations_prob, title='Mutation Probability')
    plt.savefig('{}.png'.format('ga_neural_network_size_testing'))


def test_fit_random_sa(new_df):
    """Test fit method using the random hill climbing algorithm"""
    X, Y = split_data(new_df)
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    # X, X_test, y, y_test = train_test_split(X, y, test_size=0.33)
    # ensure all data are floating point values
    X = X.to_numpy()
    X = X.astype('float32')
    X = X / 255
    Y = Y.to_numpy()


    traiing_scores_cool = list()
    testing_score_cool = list()

    cooling_exponents = list([0.1, 0.3, 0.5, 0.7, 0.9, 1.0])
    restarts = [10, 100, 250, 500, 1000, 2000]
    X_train, X_test, y_train, y_test = train_test_split(X, Y, random_state = 0)
    for c_exp in cooling_exponents:
        traiing_scores = list()
        testing_score = list()
        for iter in restarts:
            network = NeuralNetwork(hidden_nodes=[7,7], activation='relu',
                                    algorithm='simulated_annealing',
                                    bias=False, is_classifier=True,
                                    learning_rate=1, clip_max=1,
                                    restarts=c_exp,
                                    max_iters=iter,
                                    schedule=GeomDecay(1E11, c_exp, 0.0001),
                                    max_attempts=100)



            network.fit(X_train, y_train)
            y_train_prd = network.predict(X_train)

            y_test_pred = network.predict(X_test)




            testing_score.append(accuracy_score(y_train, y_train_prd))
            traiing_scores.append(accuracy_score(y_test, y_test_pred))
            print("Random Restart - {}".format(c_exp))

        traiing_scores_cool.append(traiing_scores)
        testing_score_cool.append(testing_score)
    #
    # plt_df = pd.DataFrame(
    #     {'x_values': restarts, 'Training score': traiing_scores_cool})
    for training_error in traiing_scores_cool:
        plt_df = pd.DataFrame(
            {'x_values': restarts, 'Training score': training_error})
        plt.plot('x_values', 'Training score', data=plt_df, marker='o',  markersize=12,
                 linewidth=4)
    plt.legend(cooling_exponents, title='Cooling Exponent')
    plt.savefig('{}.png'.format('sa_neural_network_size_training'))

    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    for testing_score in testing_score_cool:
        plt_df = pd.DataFrame(
            {'x_values': restarts, 'Testing score': testing_score})
        plt.plot('x_values', 'Testing score', data=plt_df, marker='o',  markersize=12,
                 linewidth=4)
    plt.legend(cooling_exponents, title='Cooling Exponent')
    plt.savefig('{}.png'.format('sa_neural_network_size_testing'))










def test_fit_random_hill_climb(new_df):
    """Test fit method using the random hill climbing algorithm"""
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(new_df)

    # X, X_test, y, y_test = train_test_split(X, y, test_size=0.33)
    # ensure all data are floating point values
    X = X.to_numpy()
    X = X.astype('float32')
    X = X / 255
    Y = Y.to_numpy()

    traiing_scores = list()
    testing_score = list()
    restarts = [10, 100, 250, 500, 1000, 2000]
    X_train, X_test, y_train, y_test = train_test_split(X, Y, random_state = 0)
    for random_restart in restarts:
        network = NeuralNetwork(hidden_nodes=[7,7], activation='relu',
                                algorithm='random_hill_climb',
                                bias=False, is_classifier=True,
                                learning_rate=1, clip_max=1,
                                restarts=random_restart,
                                max_attempts=100)



        network.fit(X_train, y_train)
        y_train_prd = network.predict(X_train)

        y_test_pred = network.predict(X_test)




        testing_score.append(accuracy_score(y_train, y_train_prd))
        traiing_scores.append(accuracy_score(y_test, y_test_pred))
        print("Random Restart - {}".format(random_restart))

    plt_df = pd.DataFrame(
        {'x_values': restarts, 'training score': traiing_scores, 'Testing Score': testing_score})
    print(plt_df)
    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Testing Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)
    plt.legend()
    plt.savefig('{}.png'.format('random_hill_neural_network_size'))




def neural_network_size(new_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    for size in range(50, len(new_df), 50):
        sizes.append(size)
        sample_size = new_df.sample(n=size)
        X, y = split_data(sample_size)
        clf = MLPClassifier(solver='sgd', alpha=1e-5, hidden_layer_sizes=(7, 7), random_state=1, activation='relu',
                            learning_rate='adaptive', max_iter=100000)

        scr = cross_validate(clf, X, y, cv=5, return_train_score=True, n_jobs=-1)
        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())
    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format('neural_network_size'))
    # plt.show()


import tensorflow as tf
import wandb
from wandb.keras import WandbCallback

CLASS_NAMES = ['Normal', 'Fraud']
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense
from keras.regularizers import l1

def neural_network_impact_ofweights(new_df):
    X, y = split_data(new_df)

    # X, X_test, y, y_test = train_test_split(X, y, test_size=0.33)
    # ensure all data are floating point values
    X = X.to_numpy()
    X = X.astype('float32')
    X = X / 255
    y = y.to_numpy()

    # encode strings to integer
    y = LabelEncoder().fit_transform(y)
    # # split into train and test datasets
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
    # Define per-fold score containers
    acc_per_fold = []
    loss_per_fold = []
    init_schemes = {'zeroes': tf.keras.initializers.Zeros()
    , 'ones' : tf.keras.initializers.Ones()
    , 'random_uniform' : tf.keras.initializers.RandomUniform(minval=0, maxval=1, seed=666)}
    n_features = X.shape[1]
    # Define the K-fold Cross Validator
    # kfold = KFold(n_splits=5, shuffle=True)

    for init_scheme in init_schemes.keys():
        f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
        model = get_training_model(n_features, init_schemes[init_scheme])
        history = model.fit(X, y, validation_split=0.40, epochs=300, batch_size=32, verbose=1, shuffle=True)

        # list all data in history
        print(history.history.keys())
        # summarize history for accuracy

        plt.plot(history.history['accuracy'], markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
        plt.plot(history.history['val_accuracy'], markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)
        plt.title('model accuracy with {} weights'.format(init_scheme))
        plt.ylabel('accuracy')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='lower right')
        # plt.show()
        plt.savefig('{}.png'.format(init_scheme))
        print("Saved file")
        # for train, test in kfold.split(X, y):
        #     model = get_training_model(n_features, init_scheme)
        #     # Fit data to model
        #     history = model.fit(X[train], y[train],
        #                         epochs=25,
        #                         verbose=1)
        #     print("Historical accuracy - %s", history.history)
        #     scores = model.evaluate(X[test], y[test], verbose=1)
        #     acc_per_fold.append(scores[1] * 100)
        #     loss_per_fold.append(scores[0])


def neural_network_wt_zero(new_df):
    X, y = split_data(new_df)

    # ensure all data are floating point values
    X = X.to_numpy()
    X = X.astype('float32')
    X = X / 255
    y = y.to_numpy()
    n_features = X.shape[1]
    X, X_test, y, y_test = train_test_split(X, y, test_size=0.33)

    model_w_zeros = get_training_model(n_features)

    # tb_callback = tf.keras.callbacks.TensorBoard(log_dir=run_dir, histogram_freq=4, write_images=True)
    #
    # model_w_zeros.fit(X, y,
    #                   validation_data=(X_test, y_test),
    #                   epochs=20, batch_size=128,
    #                   callbacks=[WandbCallback( labels=CLASS_NAMES,
    #                                            validation_data=(X_test, y_test)),
    #                              tb_callback])


def get_training_model(n_features, init_scheme):


    model = tf.keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(n_features,)),
        tf.keras.layers.Dense(10, activation='relu', kernel_initializer=init_scheme,
                              bias_initializer='zeros'),
        # tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(8, activation='relu', kernel_initializer=init_scheme,
                              bias_initializer='zeros'),
        # tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(2, activation='softmax', kernel_initializer=init_scheme,
                              bias_initializer='zeros')
    ])

    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])

    return model


def get_training_model_w_rule(n_features, init_scheme='zeros'):
    if isinstance(init_scheme, str):
        wandb.init(project='weight-initialization-tb', sync_tensorboard=True,
                   id=init_scheme)
    else:
        wandb.init(project='weight-initialization-tb', sync_tensorboard=True,
                   id=str(init_scheme))

    model = tf.keras.models.Sequential([
        tf.keras.layers.Flatten(input_shape=(n_features,)),
        tf.keras.layers.Dense(256, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(128, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(10, activation='softmax')
    ])

    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])

    return model, wandb.run.dir


def get_training_model1(n_features, init_scheme='zeros'):
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.Dense(10, activation='relu', kernel_initializer=init_scheme, input_shape=(n_features,)))
    model.add(tf.keras.layers.Dense(8, activation='relu', kernel_initializer=init_scheme))
    model.add(tf.keras.layers.Dense(2, activation='sigmoid', kernel_initializer=init_scheme))

    model.compile(optimizer='adam',
                  loss='sparse_categorical_crossentropy',
                  metrics=['accuracy'])

    return model


def neural_network(new_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(new_df)

    clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(7, 7), random_state=1)

    # Creating a visualizer object
    viz = ValidationCurve(clf, param_name="max_depth", n_jobs=-1,
                          param_range=np.arange(0, 10), cv=5, scoring="precision")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="pruning_depth.png")


def neural_network_learning_network(new_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(new_df)

    clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(7, 7), random_state=1)

    # Creating a visualizer object
    viz = ValidationCurve(clf, param_name="max_depth", n_jobs=-1,
                          param_range=np.arange(0, 10), cv=5, scoring="precision")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="pruning_depth.png")


def neural_network_grid_search(new_df):
    for i in range(10):
        X, Y = split_data(new_df)
        mlp_gs = MLPClassifier(max_iter=100000)
        parameter_space = {
            'hidden_layer_sizes': [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (6), (7), (7, 7), (8, 8), (9, 9),
                                   (10, 10)],
            'activation': ['tanh', 'relu'],
            'solver': ['sgd', 'adam'],
            'alpha': [0.0001, 0.05],
            'learning_rate': ['constant', 'adaptive'],
        }
        clf = GridSearchCV(mlp_gs, parameter_space, n_jobs=-1, cv=5)
        clf.fit(X, Y)
        print(clf.feature_names_in_)
        print(clf.best_params_)


def split_data(df):
    df = df.sample(frac=1)
    X = df.drop('Class', axis=1)
    y = df['Class']
    return X, y


from sklearn.model_selection import cross_val_score
from yellowbrick.model_selection import ValidationCurve


def decistiontree_pruning(res_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(res_df)

    dt = DecisionTreeClassifier(criterion="gini",
                                splitter="best",
                                max_depth=None,
                                min_samples_split=2,
                                min_samples_leaf=2,
                                min_weight_fraction_leaf=0.0,
                                max_features=None,
                                random_state=None,
                                max_leaf_nodes=None,
                                min_impurity_decrease=0.0,
                                class_weight=None,
                                ccp_alpha=0.0)

    # Creating a visualizer object
    viz = ValidationCurve(dt, param_name="max_depth", n_jobs=-1,
                          param_range=np.arange(0, 10), cv=5, scoring="precision")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="pruning_depth.png")
    # plt.savefig("pruning_depth.png")
    # print(res_df)


def decistiontree_min_leaf(res_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(res_df)

    dt = DecisionTreeClassifier(criterion="gini",
                                splitter="best",
                                max_depth=None,
                                min_samples_split=2,
                                min_samples_leaf=1,
                                min_weight_fraction_leaf=0.0,
                                max_features=None,
                                random_state=None,
                                max_leaf_nodes=None,
                                min_impurity_decrease=0.0,
                                class_weight=None,
                                ccp_alpha=0.0)

    # Creating a visualizer object
    viz = ValidationCurve(dt, param_name="min_samples_leaf", n_jobs=-1,
                          param_range=np.arange(1, 20), cv=5, scoring="precision")
    viz.fit(X, Y)
    viz.ax.get_lines()[1].set_label('Cross Validation Score')

    viz.show(outpath="dt_min_samples_leaf.png")
    # plt.savefig("pruning_depth.png")
    # print(res_df)


def decision_tree_training_sample(new_df):
    f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    sizes = list()
    traiing_scores = list()
    cross_val_score = list()
    for size in range(50, len(new_df), 50):
        sizes.append(size)
        sample_size = new_df.sample(n=size)
        X, y = split_data(sample_size)
        dt = DecisionTreeClassifier(criterion="gini",
                                    splitter="best",
                                    max_depth=3,
                                    min_samples_split=2,
                                    min_samples_leaf=3,
                                    min_weight_fraction_leaf=0.0,
                                    max_features=None,
                                    random_state=None,
                                    max_leaf_nodes=None,
                                    min_impurity_decrease=0.0,
                                    class_weight=None,
                                    ccp_alpha=0.0)

        scr = cross_validate(dt, X, y, cv=5, return_train_score=True, n_jobs=-1)
        traiing_scores.append(scr['train_score'].mean())
        cross_val_score.append(scr['test_score'].mean())
    plt_df = pd.DataFrame(
        {'x_values': sizes, 'training score': traiing_scores, 'Cross Validation Score': cross_val_score})

    plt.plot('x_values', 'training score', data=plt_df, marker='o', markerfacecolor='blue', markersize=12,
             color='skyblue',
             linewidth=4)
    plt.plot('x_values', 'Cross Validation Score', data=plt_df, marker='o', markerfacecolor='green', markersize=12,
             color='lightgreen',
             linewidth=4)

    # show legend
    plt.legend()
    plt.savefig('{}.png'.format('decision_tree_training_size_credit_card'))
    # plt.show()


def decisionTreeGridSearch(new_df):
    # f, (ax1) = plt.subplots(1, 1, figsize=(20, 6))
    X, Y = split_data(new_df)

    param_dict = {
        "criterion": ['gini', 'entropy'],
        "max_depth": range(1, 10),
        "min_samples_leaf": range(1, 20)
    }

    dt = DecisionTreeClassifier()

    grid = GridSearchCV(dt, param_grid=param_dict, cv=10, verbose=1, n_jobs=-1)
    grid.fit(X, Y)
    print(grid.best_params_)


def outlier(new_df, file_name):
    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(20, 6))
    ax11 = sns.swarmplot(x='Class', y='V4', data=new_df, ax=ax1, zorder=0)

    ax21 = sns.swarmplot(x='Class', y='V11', data=new_df, ax=ax2, zorder=0)

    ax31 = sns.swarmplot(x='Class', y='V10', data=new_df, ax=ax3, zorder=0)

    ax41 = sns.swarmplot(x='Class', y='V12', data=new_df, ax=ax4, zorder=0)

    ax51 = sns.swarmplot(x='Class', y='V14', data=new_df, ax=ax5, zorder=0, hue="Class")

    plt.legend(title='Transaction', loc='lower left', labels=['Normal', 'Fraud'])
    plt.savefig('{}.png'.format(file_name))


def percentile(values, low=25, high=75):
    return np.percentile(values, low), np.percentile(values, high)


def remove_outlier(new_df):
    # Normal Transactions
    new_df = remove_outlier_util(new_df, 'V4', 0)

    new_df = remove_outlier_util(new_df, 'V11', 0)

    new_df = remove_outlier_util(new_df, 'V10', 0)

    new_df = remove_outlier_util(new_df, 'V12', 0)

    # Fraud Transaction
    new_df = remove_outlier_util(new_df, 'V11', 1)

    new_df = remove_outlier_util(new_df, 'V14', 1)

    return new_df


def remove_outlier_util(new_df, feature, cat):
    V4_non_fraud = new_df[feature].loc[new_df['Class'] == cat].values
    V_25, V_75 = percentile(V4_non_fraud)
    V_iqr_cut_off = (V_75 - V_25) * 6
    v_lower, v_upper = V_25 - V_iqr_cut_off, V_75 + V_iqr_cut_off
    new_df = new_df.drop(new_df[(new_df[feature] > v_upper) | (new_df[feature] < v_lower)].index)
    return new_df


def models_evaluation(X, y, folds):
    '''
    X : data set features
    y : data set target
    folds : number of cross-validation folds

    '''

    scoring = {
               'precision': make_scorer(precision_score)}
    dtree_model = DecisionTreeClassifier(criterion='gini', max_depth=4)
    knn_model = KNeighborsClassifier(n_neighbors=5)
    hl = [5] * 2
    nn_model = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=tuple(hl), random_state=1)
    trick = 'linear'
    svm_model = SVC(kernel=trick, C=1.0, random_state=0)
    boost_model = AdaBoostClassifier(
        DecisionTreeClassifier(max_depth=2),
        n_estimators=6,
        learning_rate=1)

    # Perform cross-validation to each machine learning classifier
    start_time = time.time()
    dtree = cross_validate(dtree_model,X,y,cv=folds, scoring=scoring)
    dtree_time = time.time() - start_time
    start_time = time.time()
    knn = cross_validate(knn_model, X, y, cv=folds, scoring=scoring)
    knn_time = time.time() - start_time
    start_time = time.time()
    nn = cross_validate(nn_model, X, y, cv=folds, scoring=scoring)
    nn_time = time.time() - start_time
    start_time = time.time()
    svm = cross_validate(svm_model, X, y, cv=folds, scoring=scoring)
    svm_time = time.time() - start_time
    start_time = time.time()
    boost = cross_validate(boost_model, X, y, cv=folds, scoring=scoring)
    boost_time = time.time() - start_time

    # Create a data frame with the models perfoamnce metrics scores

    models_scores_table = pd.DataFrame( {'precision ': [dtree['test_precision'].mean(), svm['test_precision'].mean(),knn['test_precision'].mean(), boost['test_precision'].mean(), nn['test_precision'].mean()],
        }, index= ['Decision Tree', 'Support Vector Machine', 'K Nearest Neighbours', 'AdaBoost',  'Neural Network'])

    # Add 'Best Score' column
    # models_scores_table['Best Score'] = models_scores_table.idxmax(axis=0)

    # runtime_dt = pd.DataFrame(
    #                           index = ['Decision Tree','Support Vector Machine','K Nearest Neighbours','AdaBoost','Neural Network'], )

    # Return models performance metrics scores data frame
    plt.figure(figsize=(37, 8))
    models_scores_table.plot(kind="bar")
    plt.xticks(rotation=10, horizontalalignment="center")
    plt.ylabel("Precision Score")
    # plt.figure(figsize=(17, 8))
    # print(runtime_dt.head())
    # runtime_dt.plot.bar()
    plt.savefig('prediction_cc.png')


from sklearn.metrics import plot_confusion_matrix
def confusion_matrix(new_df):
    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(20, 6))
    knn_model = KNeighborsClassifier(n_neighbors=5)
    X, Y = split_data(new_df)
    X_train, X_test, y_train, y_test = train_test_split(X, Y, random_state = 0)
    knn_model.fit(X_train, y_train)
    plot_confusion_matrix(knn_model, X_test, y_test)
    plt.savefig("knn_confusion_matrix_cc.png")

    boost_model = AdaBoostClassifier(
        DecisionTreeClassifier(max_depth=2),
        n_estimators=6,
        learning_rate=1)
    boost_model.fit(X_train, y_train)
    plot_confusion_matrix(boost_model, X_test, y_test)
    plt.savefig("ada_boost_confusion_matrix_cc.png")








# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
