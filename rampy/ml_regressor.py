import sklearn
from sklearn import model_selection
from sklearn.neural_network import MLPRegressor
import pandas as pd
import numpy as np

def chemical_splitting(Pandas_DataFrame, target,split_fraction):
    """split datasets depending on their chemistry

        Parameters
        ==========
        Pandas_DataFrame: A Pandas DataFrame
            The input DataFrame with in the first row the names of the different data compositions
        target: string
            The target in the DataFrame according to which we will split the dataset
        split_fraction: a float number between 0 and 1
            This is the amount of splitting you want, in reference to the second output dataset (see OUTPUTS).

        Returns
        =======
            frame1 : A Pandas DataFrame
                A DataSet with (1-split_fraction) datas from the input dataset with unique chemical composition / names
            frame2 : A Pandas DataFrame
                A DataSet with split_fraction datas from the input dataset with unique chemical composition / names
            frame1_idx :
                A numpy array containing the indexes of the data picked in Pandas_DataFrame to construct frame1
            frame2_idx :
                A numpy array containing the indexes of the data picked in Pandas_DataFrame to construct frame2
        Note
        ====
        This function avoids the same chemical dataset to be found in different training/testing/validating datasets that are used in ML.

        Indeed, it is worthless to put data from the same original dataset / with the same chemical composition
        in the training / testing / validating datasets. This creates a initial bias in the splitting process...
    """
    names = Pandas_DataFrame[target].unique()
    names_idx = np.arange(len(names))

    # getting index for the frames with the help of scikitlearn
    frame1_idx, frame2_idx = model_selection.train_test_split(names_idx, test_size = split_fraction,random_state=42)

    # and now grabbing the relevant pandas dataframes
    ttt = np.in1d(Pandas_DataFrame.logfo2,names[frame1_idx])
    frame1 = Pandas_DataFrame[ttt == True]
    frame1_idx = np.where(ttt == True)

    ttt2 = np.in1d(Pandas_DataFrame.logfo2,names[frame2_idx])
    frame2 = Pandas_DataFrame[ttt2 == True]
    frame2_idx = np.where(ttt2 == True)

    return frame1, frame2, frame1_idx, frame2_idx

def mlregressor(x, y, algorithm="SVM",**kwargs):
    """use machine learning algorithms from scikit learn to perform regression between spectra and a observed variable.

    Parameters
    ==========
    x: array{Float64}
            the spectra organised in rows (1 row = one spectrum). The spectra should share a common X axis.
    y: Array{Float64}
            the targets. Only a single target is possible for now.
    algorithm: String,
            "KernelRidge", "SVM", "LinearRegression", "Lasso", "ElasticNet", "NeuralNet", default = "SVM"

    Returns
    =======
    prediction_train: Array{Float64}
            the predicted target values for the training y dataset.
    prediction_test: Array{Float64}
            the predicted target values for the testing y_test dataset.
    model: Scikit learn model
            A Scikit Learn object model, see scikit learn library documentation.
    X_scaler:
            A Scikit Learn scaler object for the x values.
    Y_scaler
            A Scikit Learn scaler object for the y values.

    OPTIONS
    =======

    X_test: Array{Float64}
            spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
    y_test: Array{Float64}
            the target that you want to use as a testing dataset. Those targets should not be present in the y (training) dataset.
    scaler: String
            the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see http://scikit-learn.org/stable/modules/preprocessing.html for details. Default = "MinMaxScaler".
    rand_state: Float64
            the random seed that is used for reproductibility of the results. Default = 42.
    param_grid_kr: Dictionary
            contain the values of the hyperparameters that should be checked by gridsearch for the Kernel Ridge regression algorithm.
    param_grid_svm: Dictionary
            containg the values of the hyperparameters that should be checked by gridsearch for the Support Vector regression algorithm.

    For the last two parameters, the user is refered to the documentation of SciKit Learn. See the pages:

    http://scikit-learn.org/stable/modules/kernel_ridge.html

    http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html

    NOTES
    =====

    For Support Vector and Kernel Ridge regressions, mlregressor performs a cross_validation search with using 5 KFold cross validators.

    If the results are poor with Support Vector and Kernel Ridge regressions, you will have to tune the param_grid_kr or param_grid_svm dictionnary that records the hyperparameter space to investigate during the cross validation.

    """

    #
    # Kwargs extractions
    #

    X_test = kwargs.get("X_test",[0.0])
    y_test = kwargs.get("y_test",[0.0])
    test_sz = kwargs.get("test_sz",0.3)
    scaling = kwargs.get("scaling","yes")
    scaler = kwargs.get("scaler","MinMaxScaler")
    rand_state = kwargs.get("rand_state",42)
    param_grid_kr = kwargs.get("param_grid_kr",dict(alpha=[1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3],gamma=np.logspace(-4, 4, 9)))
    param_grid_svm= kwargs.get("param_grid_svm",dict(C= [1e0, 2e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5], gamma= np.logspace(-4, 4, 9)))
    user_kernel = kwargs.get("user_kernel","rbf")

    if len(X_test) == 1:
        X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
            x, y.reshape(-1, 1), test_size=test_sz, random_state=rand_state)
    elif X_test.shape[1] == x.shape[1]:
        X_train = np.copy(x)
        y_train = np.copy(y)
    else:
        ValueError("You tried to provide a testing dataset that has a different number of features (in columns) than the training set. Please correct this.")

    # to avoid any problem with SciKit Learn quite annoying demands for the shape of arrays...
    y_train = y_train.reshape(-1,1)
    y_test = y_test.reshape(-1, 1)

    # initialising the preprocessor scaler
    if scaler == "StandardScaler":
        X_scaler = sklearn.preprocessing.StandardScaler()
        Y_scaler = sklearn.preprocessing.StandardScaler()
    elif scaler == "MinMaxScaler":
        X_scaler = sklearn.preprocessing.MinMaxScaler()
        Y_scaler = sklearn.preprocessing.MinMaxScaler()
    else:
        InputError("Choose the scaler between MinMaxScaler and StandardScaler")

    X_scaler.fit(X_train)
    Y_scaler.fit(y_train)

    # scaling the data
    X_train_sc = X_scaler.transform(X_train)
    y_train_sc = Y_scaler.transform(y_train)

    X_test_sc = X_scaler.transform(X_test)
    y_test_sc = Y_scaler.transform(y_test)

    if algorithm == "KernelRidge":
        clf_kr = KernelRidge(kernel=user_kernel, gamma=0.1)
        model = sklearn.model_selection.GridSearchCV(
            clf_kr, cv=5, param_grid=param_grid_kr)
    elif algorithm == "SVM":
        clf_svm = SVR(kernel=user_kernel, gamma=0.1)
        model = sklearn.model_selection.GridSearchCV(
            clf_svm, cv=5, param_grid=param_grid_svm)
    elif algorithm == "Lasso":
        clf_lasso = sklearn.linear_model.Lasso(alpha=0.1)
        model = sklearn.model_selection.GridSearchCV(clf_lasso, cv=5, param_grid=dict(alpha=[1e-3,1e-2,1e-1,1.,1e1,1e2,1e3,1e4]))
    elif algorithm == "ElasticNet":
        clf_ElasticNet = sklearn.linear_model.ElasticNet(alpha=0.1, l1_ratio=0.5)
        model = sklearn.model_selection.GridSearchCV(clf_ElasticNet, cv=5, param_grid=dict(alpha=[1e-3, 1e-2, 1e-1, 1., 1e1, 1e2, 1e3, 1e4]))
    elif algorithm == "LinearRegression":
        model = sklearn.linear_model.LinearRegression()
    elif algorithm == "NeuralNet":
        model = MLPRegressor(hidden_layer_sizes=(2,), activation='relu',solver='lbfgs', random_state=rand_state)

    if scaling == "yes":
        model.fit(X_train_sc, y_train_sc)
        predict_train_sc = model.predict(X_train_sc)
        prediction_train = Y_scaler.inverse_transform(predict_train_sc)
        predict_test_sc = model.predict(X_test_sc)
        prediction_test = Y_scaler.inverse_transform(predict_test_sc)
        return prediction_train, prediction_test, model, X_scaler, Y_scaler
    else:
        model.fit(X_train, vec(y_train))
        prediction_train = model.predict(X_train)
        prediction_test = model.predict(X_test)
        return prediction_train, prediction_test, model
