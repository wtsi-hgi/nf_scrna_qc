#!/usr/bin/env python

__date__ = '2020-04-24'
__version__ = '0.0.1'


import argparse
import os
import random
import numpy as np
import scipy as sp
import pandas as pd
import scanpy as sc
import csv
from distutils.version import LooseVersion

# import joblib  # for numpy matrix, joblib faster than pickle
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import plotnine as plt9

from sklearn import metrics
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.model_selection import GridSearchCV

import keras
from keras.utils import np_utils
from keras.models import Sequential
from keras.layers import Dense
from keras.regularizers import L1L2
from keras.wrappers.scikit_learn import KerasClassifier

from tensorflow.python.client import device_lib
import tensorflow as tf

# Check that we are working on GPU or CPU
# print(device_lib.list_local_devices())  # list of DeviceAttributes
# tf.config.list_physical_devices('GPU')

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)
# 3. Set the `tensorflow` pseudo-random generator at a fixed value
tf.random.set_seed(seed_value)

# Get compression opts for pandas
compression_opts = 'gzip'
if LooseVersion(pd.__version__) > '1.0.0':
    compression_opts = dict(method='gzip', compresslevel=9)


def _create_colors(classes):
    n_cts = len(classes)
    color_norm = colors.Normalize(vmin=-n_cts / 3, vmax=n_cts)
    ct_arr = np.arange(n_cts)
    ct_colors = cm.YlGnBu(color_norm(ct_arr))

    return ct_colors


def plot_roc(y_prob, y_test, classes):
    """Plot ROC curve. Based off of NaiveDE library."""
    ct_colors = _create_colors(classes)

    for i, cell_type in enumerate(classes):
        fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
        plt.plot(fpr, tpr, c=ct_colors[i], lw=2)

    plt.plot([0, 1], [0, 1], color='k', ls=':')
    plt.xlabel('FPR')
    plt.ylabel('TPR')


def class_report(y_true, y_pred, classes, y_pred_proba=None):
    """
    Build a text report showing the main classification metrics.

    Replaces sklearn.metrics.classification_report.

    Derived from:
    https://stackoverflow.com/questions/39685740/calculate-sklearn-roc-auc-score-for-multi-class
    """
    if y_true.shape != y_pred.shape:
        raise Exception(
            'Error! y_true {} is not the same shape as y_pred {}'.format(
                y_true.shape,
                y_pred.shape
            )
        )

    # NOTE: Y may not have predictions for all classes
    model_report = pd.DataFrame(classification_report(
        y_true,
        y_pred,
        classes,
        output_dict=True
    )).transpose()

    if not (y_pred_proba is None):
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        aupc = dict()
        mcc = dict()
        for label_it, label in enumerate(model_report.index):
            if label in classes:  # skip accuracy, macro avg, weighted avg
                fpr[label], tpr[label], _ = metrics.roc_curve(
                    (y_true == label).astype(int),
                    y_pred_proba[:, label_it]
                )
                roc_auc[label] = metrics.auc(fpr[label], tpr[label])
                aupc[label] = metrics.average_precision_score(
                    (y_true == label).astype(int),
                    y_pred_proba[:, label_it],
                    average=None  # No need since iter over labels
                )
                mcc[label] = metrics.matthews_corrcoef(
                    (y_true == label).astype(int),
                    (y_pred == label).astype(int)
                )
            else:
                fpr[label] = np.nan
                tpr[label] = np.nan
                roc_auc[label] = np.nan
                aupc[label] = np.nan
                mcc[label] = np.nan
        model_report['AUC'] = pd.Series(roc_auc)
        model_report['average_precision_score'] = pd.Series(aupc)
        model_report['MCC'] = pd.Series(mcc)

    # Catch the case where true label not predicted in lr, perhaps because
    # too few training cases.
    for i in np.unique(y_true):
        if i not in model_report.index:
            print(
                'Adding category ({}) from {}.'.format(
                    i,
                    'truth with no prediction to report'
                )
            )
            model_report = model_report.append(pd.Series(
                [np.nan]*len(model_report.columns),
                index=model_report.columns,
                name=i
            ))

    model_report = model_report.sort_index()

    return model_report


def keras_grid(
    model_function,
    encoder,
    X_std,
    y,
    n_epochs=100,
    batch_size=32
):
    # Run same proceedure on the test data
    y_encoded = encoder.transform(y)
    Y_onehot = np_utils.to_categorical(y_encoded)

    # Initial parameter sweep for different activation, optimizer, and loss.
    # NOTE: From 100k cells from 10x, best settings were:
    # 'activation': 'softmax',
    # 'loss': 'categorical_crossentropy',
    # 'optimizer': toss up between adam and sgd, though sgd generally better
    # 'sparsity_l1': 0.001
    # param_grid = dict(
    #     activation=['softmax', 'sigmoid'],
    #     optimizer=['sgd', 'adam'],
    #     loss=['categorical_crossentropy', 'mean_squared_error'],
    #     sparsity_l1=[0.1, 0.01, 0.001, 0.0005]
    # )
    # NOTE: sparse_categorical_crossentropy is for classes that are not one
    # hot encoded.
    # https://www.quora.com/What-is-the-difference-between-categorical_crossentropy-and-sparse_categorical-cross-entropy-when-we-do-multiclass-classification-using-convolution-neural-networks
    # param_grid = dict(
    #     activation=['softmax'],
    #     optimizer=['sgd'],
    #     loss=['categorical_crossentropy'],
    #     sparsity_l2__activity=[0.0, 1e-6],
    #     sparsity_l1__activity=[0.1, 1e-4, 1e-10, 0.0],
    #     sparsity_l2__kernel=[0.0, 1e-6],
    #     sparsity_l1__kernel=[0.1, 1e-4, 1e-10, 0.0],
    #     sparsity_l2__bias=[0.0, 1e-6],
    #     sparsity_l1__bias=[0.1, 1e-4, 1e-10, 0.0]
    # )
    param_grid = dict(
        activation=['softmax'],
        optimizer=['sgd'],
        loss=['categorical_crossentropy'],
        sparsity_l2__activity=[0.0],
        sparsity_l1__activity=[0.1, 1e-4],
        sparsity_l2__kernel=[0.0],
        sparsity_l1__kernel=[0.1, 1e-4],
        sparsity_l2__bias=[0.0],
        sparsity_l1__bias=[0.1, 1e-4]
    )
    n_splits = 5
    grid = GridSearchCV(
        estimator=KerasClassifier(build_fn=model_function),
        param_grid=param_grid,
        n_jobs=1,
        cv=n_splits  # Number of cross validation.
    )
    # NOTE: We could pass batch_size and epochs here, but we get results much
    # faster if we just use the defaults.
    grid_result = grid.fit(
        # batch_size=batch_size,
        # epochs=n_epochs,
        X=X_std,
        y=Y_onehot
    )

    # Make a dataframe of the results of all of the models.
    cv_results = grid_result.cv_results_.copy()
    del cv_results['param_activation']
    df_grid_result = pd.DataFrame(cv_results.pop('params'))
    # Rename so we know that these columns are parameters
    df_grid_result.columns = [
        'param__{}'.format(i) for i in df_grid_result.columns
    ]
    df_grid_result = pd.concat([
        df_grid_result,
        pd.DataFrame(cv_results)
    ], axis=1)

    print('Best: %f using %s' % (
        grid_result.best_score_,
        grid_result.best_params_
    ))

    return grid_result, df_grid_result


def fit_model_keras(
    model_function,
    encoder,
    X_std,
    y,
    sparsity_l1=0.01,
    sparsity_l2=0.0,
    n_epochs=100,
    batch_size=32,
    train_size_fraction=0.67,
    verbose=True
):
    # References:
    # https://machinelearningmastery.com/multi-class-classification-tutorial-keras-deep-learning-library/
    # https://stackoverflow.com/questions/59643062/scikit-learn-vs-keras-tensorflow-for-multinomial-logistic-regression
    # https://medium.com/@luwei.io/logistic-regression-with-keras-d75d640d175e

    # Make the training and test dataset
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X_std,
        y,
        stratify=y,
        random_state=61,
        train_size=train_size_fraction
    )
    if verbose:
        print(
            'Split X into training {} and test {} sets.'.format(
                X_train.shape,
                X_test.shape
            )
        )

    # One hot encode y (the cell type classes)
    # encode class values as integers
    # encoder = preprocessing.LabelEncoder()
    # encoder.fit(y_train)
    y_train_encoded = encoder.transform(y_train)
    # convert integers to dummy variables (i.e. one hot encoded)
    Y_train_onehot = np_utils.to_categorical(y_train_encoded)
    # Run same proceedure on the test data
    y_test_encoded = encoder.transform(y_test)
    Y_test_onehot = np_utils.to_categorical(y_test_encoded)

    # Training
    model = model_function(
        sparsity_l1__activity=sparsity_l1,
        sparsity_l2__activity=sparsity_l2,
        sparsity_l1__kernel=sparsity_l1,
        sparsity_l2__kernel=sparsity_l2,
        sparsity_l1__bias=sparsity_l1,
        sparsity_l2__bias=sparsity_l2
    )
    history = model.fit(
        X_train,
        Y_train_onehot,
        batch_size=batch_size,
        epochs=n_epochs,
        verbose=0,
        # use_multiprocessing=True,
        # validation_split=0.33  # Frac of the training used for validation.
        validation_data=(X_test, Y_test_onehot)
    )

    # Train using KFold validation
    # from keras.wrappers.scikit_learn import KerasClassifier
    # from sklearn.model_selection import KFold
    # from sklearn.model_selection import cross_val_score
    # estimator = KerasClassifier(
    #     build_fn=classification_model,
    #     epochs=200,
    #     # batch_size=5,
    #     verbose=1
    # )
    # kfold = KFold(n_splits=10, shuffle=True)
    # results = cross_val_score(estimator, X_std, y_onehot, cv=kfold)
    # print("Baseline: %.2f%% (%.2f%%)" % (
    #     results.mean()*100, results.std()*100)
    # )

    # Make a classifier report
    classes = np.argmax(model.predict(X_test), axis=1)
    y_test_pred = encoder.inverse_transform(classes)
    y_test_proba = model.predict_proba(X_test)
    model_report = class_report(
        y_test,
        y_test_pred,
        encoder.classes_,
        y_test_proba
    )
    # Add the number of cells in each class (index) in the
    # (a) full dataset and (b) training dataset.
    categories, counts = np.unique(y, return_counts=True)
    cat_counts = dict(zip(categories, counts))
    model_report['n_cells_full_dataset'] = model_report.index.map(cat_counts)
    categories, counts = np.unique(y_train, return_counts=True)
    cat_counts = dict(zip(categories, counts))
    model_report['n_cells_training_dataset'] = model_report.index.map(
        cat_counts
    )

    # Get a matrix of predictions on the test set
    y_prob_df = pd.DataFrame(
        y_test_proba,
        columns=['class__{}'.format(i) for i in encoder.classes_]
    )
    y_prob_df['cell_label_predicted'] = y_test_pred
    y_prob_df['cell_label_true'] = y_test
    for i in ['cell_label_predicted', 'cell_label_true']:
        y_prob_df[i] = 'class__' + y_prob_df[i].astype(str)

    score = model.evaluate(X_test, Y_test_onehot, verbose=0)
    print('Test score:', score[0])
    print('Test accuracy:', score[1])

    return model, model_report, y_prob_df, history


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Fits logistic regression to predict labels.'
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file where clusters have been saved to cluster slot.'
    )

    parser.add_argument(
        '-s', '--sparsity_l1',
        action='store',
        dest='sparsity_l1',
        default=0.0001,
        type=float,
        help='Smaller values specify stronger regularization.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-nepoch', '--number_epoch',
        action='store',
        dest='number_epoch',
        default=25,
        type=int,
        help='Number of epochs.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-bs', '--batch_size',
        action='store',
        dest='batch_size',
        default=32,
        type=int,
        help='Batch size. Divides the dataset into n batches and updates the\
            weights at the end of each one.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-tsc', '--train_size_cells',
        action='store',
        dest='train_size_cells',
        default=0,
        type=int,
        help='Number of cells to use for training set. If > 0 all\
            remaining cells not randomly selected for training will be used\
            for the test set. Overrides <train_size_fraction>.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-tsf', '--train_size_fraction',
        action='store',
        dest='train_size_fraction',
        default=0.67,
        type=float,
        help='Fraction of the data to use for training set.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--dict_add',
        action='store',
        dest='dict_add',
        default='',
        type=str,
        help='Additional information to add to output model_report.\
            Format: key::value:::key2::value2.\
            Example: method::leiden:::resolution::3.0\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--grid_search',
        action='store_true',
        dest='grid_search',
        default=False,
        help='Run a grid search of hyperparameters.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--memory_limit',
        action='store',
        dest='memory_limit',
        default=50,
        type=int,
        help='Memory limit in Gb.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: keras_model-<params>)'
    )
    options = parser.parse_args()

    verbose = True

    # Set GPU memory limits
    gpus = tf.config.list_physical_devices('GPU')
    print(gpus)
    if gpus:
        # For TF v1
        # config = tf.ConfigProto()
        # config.gpu_options.allow_growth = True
        # session = tf.Session(config=config)

        # For TF v2
        try:
            # Method 1:
            # Currently, memory growth needs to be the same across GPUs
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)

            # Method 2:
            # Restrict TensorFlow to only allocate 1GB of memory on the first
            # GPU
            # tf.config.experimental.set_virtual_device_configuration(
            #     gpus[0],
            #     [tf.config.experimental.VirtualDeviceConfiguration(
            #         memory_limit=options.memory_limit*1024
            #     )])
            # logical_gpus = tf.config.list_logical_devices('GPU')
            # print(
            #     len(gpus),
            #     "Physical GPUs,",
            #     len(logical_gpus),
            #     "Logical GPUs"
            # )
        except RuntimeError as e:
            # Virtual devices must be set before GPUs have been initialized
            print(e)
    else:
        raise Exception('ERROR: no GPUs detected.')

    # Get additional data we are going to append to the output model info
    dict_add = {}
    if options.dict_add != '':
        for item in options.dict_add.split(':::'):
            _tmp = item.split('::')
            if len(_tmp) != 2:
                raise Exception('ERROR: check dict_add.')
            else:
                dict_add[_tmp[0]] = _tmp[1]
    print(dict_add)

    # Load the AnnData file.
    # This file should already have clusters identified and saved to the
    # clusters slot.
    adata = sc.read_h5ad(filename=options.h5)

    # Set X to cp10k
    # adata.X = np.expm1(adata.layers['log1p_cp10k'])
    # Set X to ln(cp10k+1)
    # NOTE: Testing with 100k cells from 10x, we were able to achieve
    # higher accuracy with log1p_cp10k - likely becuase better spread in
    # distribution.
    adata.X = adata.layers['log1p_cp10k']
    # Set X to raw counts
    # adata.X = adata.layers['counts']

    # Add some info from adata to dict_add
    for key, value in adata.uns['neighbors']['params'].items():
        dict_add['neighbors__{}'.format(key)] = value
    for key, value in adata.uns['cluster']['params'].items():
        dict_add['cluster__{}'.format(key)] = value

    # If train_size_cells, override the fraction so that the total number of
    # cells in the training set will be equal to train_size_cells.
    train_size_fraction = options.train_size_fraction
    if options.train_size_cells > 0:
        if options.train_size_cells >= adata.n_obs:
            raise Exception('Invalid train_size_cells.')
        train_size_fraction = (
            1 - ((adata.n_obs-options.train_size_cells)/adata.n_obs)
        )
        if verbose:
            print('Set train_size_fraction to: {}.'.format(
                train_size_fraction
            ))
    if verbose:
        print('Number cells training ({}) and testing ({}).'.format(
            int(train_size_fraction*adata.n_obs),
            int((1-train_size_fraction)*adata.n_obs)
        ))

    # Set X and y
    X = adata.X
    y = adata.obs['cluster'].values

    # Set other variables
    sparsity_l1 = options.sparsity_l1
    n_epochs = options.number_epoch
    batch_size = options.batch_size

    # Center and scale the data
    if sp.sparse.issparse(X):
        X = X.todense()
    X_std = X
    scaler = preprocessing.StandardScaler(
        with_mean=True,
        with_std=True
    )
    X_std = scaler.fit_transform(X)
    if verbose:
        print('center={} scale={}'.format(
            True,
            True
        ))

    # One hot encode y (the cell type classes)
    # encode class values as integers
    encoder = preprocessing.LabelEncoder()
    encoder.fit(y)
    print('Found {} clusters'.format(len(encoder.classes_)))

    # Define the model
    # NOTE: Defaults determined via grid search of 160k TI single cells
    def classification_model(
        optimizer='sgd',
        activation='softmax',
        loss='categorical_crossentropy',
        sparsity_l1__activity=0.0001,
        sparsity_l2__activity=0.0,
        sparsity_l1__kernel=0.0,
        sparsity_l2__kernel=0.0,
        sparsity_l1__bias=0.0,
        sparsity_l2__bias=0.0
    ):
        # create model
        model = Sequential()
        # Use a “softmax” activation function in the output layer. This is to
        # ensure the output values are in the range of 0 and 1 and may be used
        # as predicted probabilities.
        #
        # https://developers.google.com/machine-learning/crash-course/multi-class-neural-networks/softmax
        # Softmax assigns decimal probabilities to each class in a multi-class
        # problem. Those decimal probabilities must add up to 1.0. This
        # additional constraint helps training converge more quickly than it
        # otherwise would. Softmax is implemented through a neural network
        # layer just before the output layer. The Softmax layer must have the
        # same number of nodes as the output layer.
        # Softmax assumes that each example is a member of exactly one class.
        #
        # Softmax should be used for multi-class prediction with single label
        # https://developers.google.com/machine-learning/crash-course/multi-class-neural-networks/video-lecture
        # NOTE: input dimension = number of features your data has
        model.add(Dense(
            len(encoder.classes_),  # output dim is number of classes
            use_bias=True,  # intercept
            activation=activation,  # softmax, sigmoid
            activity_regularizer=L1L2(
                l1=sparsity_l1__activity,
                l2=sparsity_l2__activity
            ),
            kernel_regularizer=L1L2(
                l1=sparsity_l1__kernel,
                l2=sparsity_l2__kernel
            ),
            bias_regularizer=L1L2(
                l1=sparsity_l1__bias,
                l2=sparsity_l2__bias
            ),
            input_dim=X.shape[1]
        ))
        # Example of adding additional layers
        # model.add(Dense(8, input_dim=4, activation='relu'))
        # model.add(Dense(3, activation='softmax'))

        # Metrics to check out over training epochs
        mets = [
            # loss,
            keras.metrics.CategoricalAccuracy(name='categorical_accuracy'),
            # keras.metrics.TruePositives(name='tp'),
            # keras.metrics.FalsePositives(name='fp'),
            # keras.metrics.TrueNegatives(name='tn'),
            # keras.metrics.FalseNegatives(name='fn'),
            # keras.metrics.Precision(name='precision'),
            # keras.metrics.Recall(name='recall'),
            # keras.metrics.AUC(name='auc'),
            keras.metrics.BinaryAccuracy(name='accuracy')
        ]
        # Use Adam gradient descent optimization algorithm with a logarithmic
        # loss function, which is called “categorical_crossentropy” in Keras.
        # UPDATE: sgd works better emperically.
        model.compile(
            optimizer=optimizer,  # adam, sgd
            loss=loss,
            metrics=mets
        )

        return model

    # Now, either call a grid search or specific model fit
    if options.grid_search:
        # Get the out file base.
        out_file_base = options.of
        if out_file_base == '':
            out_file_base = 'keras_model'
        out_file_base = '{}-grid_search'.format(out_file_base)

        # Call grid search of various parameters
        grid_result, df_grid_result = keras_grid(
            model_function=classification_model,
            encoder=encoder,
            X_std=X_std,
            y=y,
            n_epochs=n_epochs,
            batch_size=batch_size
        )

        # NOTE: This will fail because can't pickle KerasClassifier. This is
        # fine though becuase results are saved in tsv.gz format below.
        # Save the results
        # out_f = '{}-grid_result.gz'.format(out_file_base)
        # joblib.dump(
        #     grid_result,
        #     out_f,
        #     compress=('gzip', 3)
        # )
        # Load the model
        # lr = joblib.load(
        #     'test-lr_model.joblib.gz'
        # )
        # print(lr)

        # Save the results of our search to tsv
        out_f = '{}-grid_result.tsv.gz'.format(out_file_base)
        df_grid_result.to_csv(
            out_f,
            sep='\t',
            index=False,
            quoting=csv.QUOTE_NONNUMERIC,
            na_rep='',
            compression=compression_opts
        )

        # Add a single columns that summarizes params
        param_columns = [
            col for col in df_grid_result.columns if 'param__' in col
        ]
        df_grid_result['params'] = df_grid_result[
            param_columns
        ].astype(str).apply(lambda x: '-'.join(x), axis=1)

        # Plot the distribution of accuracy across folds
        split_columns = [
            col for col in df_grid_result.columns if 'split' in col
        ]
        split_columns = [
            col for col in split_columns if '_test_score' in col
        ]
        df_plt = pd.melt(
            df_grid_result,
            id_vars=['params'],
            value_vars=split_columns
        )
        gplt = plt9.ggplot(df_plt, plt9.aes(
            x='params',
            y='value'
        ))
        gplt = gplt + plt9.theme_bw()
        gplt = gplt + plt9.geom_boxplot(alpha=0.8)
        gplt = gplt + plt9.geom_jitter(alpha=0.75)
        gplt = gplt + plt9.scale_y_continuous(
            # trans='log10',
            # labels=comma_labels,
            minor_breaks=0
            # limits=[0, 1]
        )
        gplt = gplt + plt9.labs(
            x='Parameters',
            y='Score',
            title=''
        )
        gplt = gplt + plt9.theme(
            axis_text_x=plt9.element_text(angle=-45, hjust=0)
        )
        gplt.save(
            '{}-score.png'.format(out_file_base),
            # dpi=300,
            width=10,
            height=4,
            limitsize=False
        )

        # Plot the mean time and std err for fitting results
        gplt = plt9.ggplot(df_grid_result, plt9.aes(
            x='params',
            y='mean_fit_time'
        ))
        gplt = gplt + plt9.theme_bw()
        gplt = gplt + plt9.geom_point()
        gplt = gplt + plt9.geom_errorbar(
            plt9.aes(
                ymin='mean_fit_time-std_fit_time',
                ymax='mean_fit_time+std_fit_time'
            ),
            width=0.2,
            position=plt9.position_dodge(0.05)
        )
        gplt = gplt + plt9.scale_y_continuous(
            # trans='log10',
            # labels=comma_labels,
            minor_breaks=0
        )
        gplt = gplt + plt9.labs(
            x='Parameters',
            y='Mean fit time',
            title=''
        )
        gplt = gplt + plt9.theme(
            axis_text_x=plt9.element_text(angle=-45, hjust=0)
        )
        gplt.save(
            '{}-fit_time.png'.format(out_file_base),
            # dpi=300,
            width=10,
            height=4,
            limitsize=False
        )

    else:
        # Get the out file base.
        out_file_base = options.of
        if out_file_base == '':
            out_file_base = 'keras_model'
            # out_file_base = '{}-center={}-scale={}'.format(
            #     out_file_base,
            #     center,
            #     scale
            # )
            out_file_base = '{}-batch_size={}-epochs={}'.format(
                out_file_base,
                batch_size,
                n_epochs
            )
            out_file_base = '{}-sparsity_l1={}-train_size_fraction={}'.format(
                out_file_base,
                str(sparsity_l1).replace('.', 'pt'),
                str(train_size_fraction).replace('.', 'pt')
            )

        # Fit the specific model and save the results
        model, model_report, y_prob_df, history = fit_model_keras(
            model_function=classification_model,
            encoder=encoder,
            X_std=X_std,
            y=y,
            sparsity_l1=sparsity_l1,
            sparsity_l2=0.0,
            n_epochs=n_epochs,
            batch_size=batch_size,
            train_size_fraction=train_size_fraction
        )

        # Save the model, weights (coefficients), and bias (intercept)
        model.save(
            '{}.h5'.format(out_file_base),
            overwrite=True,
            include_optimizer=True
        )

        # Save the model and weights (coefficients) seperately
        # open('{}.json'.format(out_file_base), 'w').write(model.to_json())
        open('{}.yml'.format(out_file_base), 'w').write(model.to_yaml())
        model.save_weights('{}-weights.h5'.format(out_file_base))
        # Example read functions
        # model = model_from_yaml(open('my_model_architecture.yaml').read())
        # model.load_weights('my_model_weights.h5')

        # Save the model report
        # Add column telling us if this is cluster or summary value
        is_cluster = []
        for i in model_report.index:
            if i in encoder.classes_:
                is_cluster.append(True)
            else:
                is_cluster.append(False)
        model_report['is_cluster'] = is_cluster
        # Add in extra data
        model_report['sparsity_l1'] = sparsity_l1
        if dict_add:
            for key, value in dict_add.items():
                model_report[key] = value
        print(model_report)
        out_f = '{}-model_report.tsv.gz'.format(out_file_base)
        model_report.to_csv(
            out_f,
            sep='\t',
            index=True,
            index_label='cell_label',
            quoting=csv.QUOTE_NONNUMERIC,
            na_rep='',
            compression=compression_opts
        )
        if verbose:
            print('Completed: save {}.'.format(out_f))

        # Save the test results - each row is a cell and the columns are the
        # prob of that cell belonging to a particular class.
        # Add in extra data
        y_prob_df['sparsity_l1'] = sparsity_l1
        if dict_add:
            for key, value in dict_add.items():
                y_prob_df[key] = value
        out_f = '{}-test_result.tsv.gz'.format(out_file_base)
        y_prob_df.to_csv(
            out_f,
            sep='\t',
            index=False,   # NOTE: Not adding the label to test_result index.
            # index_label='cell_label',
            quoting=csv.QUOTE_NONNUMERIC,
            na_rep='',
            compression=compression_opts
        )
        if verbose:
            print('Completed: save {}.'.format(out_f))

        # Make a matrix of weights per gene
        # Columns = genes tested and rows = cell type label
        weight, bias = model.layers[-1].get_weights()
        # weight, bias = model.get_layer("output").get_weights()
        df_weights = pd.DataFrame.from_records(
            weight,
            index=adata.var.index,  # index is gene
            columns=encoder.classes_
        )
        # Save the weights dataframe.
        out_f = '{}-weights.tsv.gz'.format(out_file_base)
        df_weights.to_csv(
            out_f,
            sep='\t',
            index=True,
            index_label='ensembl_gene_id',
            quoting=csv.QUOTE_NONNUMERIC,
            na_rep='',
            compression=compression_opts
        )
        if verbose:
            print('Completed: save {}.'.format(out_f))

        # Plot the number of features with non-zero coefficients in each
        # cluster.
        out_f = '{}-n_features.png'.format(out_file_base)
        df_plt = pd.DataFrame({
            'classes': df_weights.columns,
            'features': (df_weights != 0).sum(axis=0)
        })
        df_plt = df_plt.set_index('classes')
        # print(df_plt)
        # Add in catgories with no predictive model (e.g., becuase they were
        # too few in training).
        for i in adata.obs['cluster'].cat.categories:
            if i not in df_plt.index:
                df_plt = df_plt.append(pd.Series(
                    [0],
                    index=df_plt.columns,
                    name=i
                ))
        fig = plt.figure(figsize=(max(0.5*len(df_plt.index), 5), 4))
        # plt.bar(lr.classes_, n_features)
        plt.bar(df_plt.index, df_plt['features'])
        plt.xlabel('Cluster')
        plt.ylabel('Features with coefficient != 0')
        plt.xticks(rotation=90)
        for i in df_plt.index:
            plt.annotate(
                str(df_plt.loc[i, 'features']),
                xy=(i, df_plt.loc[i, 'features'])
            )
        fig.savefig(
            out_f,
            # dpi=300,
            bbox_inches='tight'
        )
        plt.close(fig)

        # Plot ROC of the test and truth.
        out_f = '{}-roc.png'.format(out_file_base)
        fig = plt.figure()
        cell_label_true = y_prob_df.pop('cell_label_true')
        # Drop columns that are not cell type labels
        for i in y_prob_df.columns:
            if 'class__' not in i:
                del y_prob_df[i]
        plot_roc(y_prob_df.values, cell_label_true.values, y_prob_df.columns)
        fig.savefig(
            out_f,
            # dpi=300,
            bbox_inches='tight'
        )
        plt.close(fig)
        if verbose:
            print('Completed: save {}.'.format(out_f))

        # Plot metrics vs cluster size to see if smaller clusters have poorer
        # metric measures.
        df_plt = model_report.fillna(0)
        for i in df_plt.index:
            if i not in encoder.classes_:
                df_plt = df_plt.drop(i)
        for i in ['AUC', 'f1-score', 'average_precision_score', 'MCC']:
            out_f = '{}-cluster_size_{}.png'.format(out_file_base, i)
            fig = plt.figure()
            plt.scatter(df_plt['n_cells_full_dataset'], df_plt[i], alpha=0.5)
            plt.xlabel('Number of cells in cluster (full dataset)')
            plt.ylabel(i)
            if i in ['AUC', 'f1-score', 'average_precision_score']:
                plt.ylim(0, 1)
            elif i == 'MCC':
                plt.ylim(-1, 1)
            # Add annotation of the cluster
            for index, row in df_plt.iterrows():
                if row['n_cells_full_dataset'] == 0:
                    print(
                        'ERROP: n_cells_full_dataset = 0 for {}.'.format(index)
                    )
                plt.annotate(
                    index,  # this is the text
                    (row['n_cells_full_dataset'], row[i]),  # point to label
                    textcoords='offset points',  # how to position the text
                    xytext=(0, 10),  # distance from text to points (x,y)
                    ha='center'  # horiz alignment can be left, right, center
                )
            fig.savefig(
                out_f,
                # dpi=300,
                bbox_inches='tight'
            )
            plt.xscale('log', basex=10)
            fig.savefig(
                '{}-cluster_size_{}_log10.png'.format(out_file_base, i),
                # dpi=300,
                bbox_inches='tight'
            )
            plt.close(fig)
            if verbose:
                print('Completed: save {}.'.format(out_f))

        # Plot history of metrics over epochs
        for dat_i in history.history.keys():
            fig = plt.figure()
            plt.plot(history.history[dat_i])
            plt.ylabel(dat_i)
            plt.xlabel('Epoch')
            fig.savefig(
                '{}-model_iter_{}.png'.format(out_file_base, dat_i),
                # dpi=300,
                bbox_inches='tight'
            )
            plt.close(fig)


if __name__ == '__main__':
    main()
