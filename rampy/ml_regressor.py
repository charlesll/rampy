import sklearn
from sklearn import model_selection
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.ensemble import BaggingRegressor

import pandas as pd
import numpy as np

def chemical_splitting(Pandas_DataFrame, target,split_fraction=0.30, rand_state=42):
    """split datasets depending on their chemistry

        Parameters
        ==========
        Pandas_DataFrame: A Pandas DataFrame
            The input DataFrame with in the first row the names of the different data compositions
        target: string
            The target in the DataFrame according to which we will split the dataset
        split_fraction: a float number between 0 and 1
            This is the amount of splitting you want, in reference to the second output dataset (see OUTPUTS).
        rand_state : Float64
            the random seed that is used for reproductibility of the results. Default = 42.

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
    frame1_idx, frame2_idx = model_selection.train_test_split(names_idx, test_size = split_fraction,random_state=rand_state)

    # and now grabbing the relevant pandas dataframes
    ttt = np.in1d(Pandas_DataFrame[target],names[frame1_idx])
    frame1 = Pandas_DataFrame[ttt == True]
    frame1_idx = np.where(ttt == True)

    ttt2 = np.in1d(Pandas_DataFrame[target],names[frame2_idx])
    frame2 = Pandas_DataFrame[ttt2 == True]
    frame2_idx = np.where(ttt2 == True)

    return frame1, frame2, frame1_idx, frame2_idx

class mlregressor:
    """use machine learning algorithms from scikit learn to perform regression between spectra and an observed variable.


    Attributes
    ----------
    x : {array-like, sparse matrix}, shape = (n_samples, n_features)
        Spectra; n_features = n_frequencies.
    y : array, shape = (n_samples,)
        Returns predicted values.
    X_test : {array-like, sparse matrix}, shape = (n_samples, n_features)
        spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
    y_test : array, shape = (n_samples,)
        the target that you want to use as a testing dataset. Those targets should not be present in the y (training) dataset.
    algorithm : String,
        "KernelRidge", "SVM", "LinearRegression", "Lasso", "ElasticNet", "NeuralNet", "BaggingNeuralNet", default = "SVM"
    scaling : Bool
        True or False. If True, data will be scaled during fitting and prediction with the requested scaler (see below),
    scaler : String
        the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see http://scikit-learn.org/stable/modules/preprocessing.html for details. Default = "MinMaxScaler".
    test_size : float
        the fraction of the dataset to use as a testing dataset; only used if X_test and y_test are not provided.
    rand_state : Float64
        the random seed that is used for reproductibility of the results. Default = 42.
    param_kr : Dictionary
        contain the values of the hyperparameters that should be provided to KernelRidge and GridSearch for the Kernel Ridge regression algorithm.
    param_svm : Dictionary
        containg the values of the hyperparameters that should be provided to SVM and GridSearch for the Support Vector regression algorithm.
    param_neurons : Dictionary
        contains the parameters for the Neural Network (MLPregressor model in sklearn).
        Default= dict(layers=(3,),solver = 'lbfgs',funct='relu',early_stopping=True)
    param_bagging : Dictionary
        contains the parameters for the BaggingRegressor sklearn function that uses a MLPregressor base method.
        Default= dict(n_estimators=100, max_samples=1.0, max_features=1.0, bootstrap=True,
                      bootstrap_features=False, oob_score=False, warm_start=False, n_jobs=1, random_state=rand_state, verbose=0)
    prediction_train : Array{Float64}
        the predicted target values for the training y dataset.
    prediction_test : Array{Float64}
        the predicted target values for the testing y_test dataset.
    model : Scikit learn model
        A Scikit Learn object model, see scikit learn library documentation.
    X_scaler :
        A Scikit Learn scaler object for the x values.
    Y_scaler :
        A Scikit Learn scaler object for the y values.

    Remarks
    -------

    For details on hyperparameters of each algorithms, please directly consult the documentation of SciKit Learn at:

    http://scikit-learn.org/stable/

    For Support Vector and Kernel Ridge regressions, mlregressor performs a cross_validation search with using 5 KFold cross validators.

    If the results are poor with Support Vector and Kernel Ridge regressions, you will have to tune the param_grid_kr or param_grid_svm dictionnary that records the hyperparameter space to investigate during the cross validation.

    Results for machine learning algorithms can vary from run to run. A way to solve that is to fix the random_state. 
    For neural nets, results from multiple neural nets (bagging technique) may also generalise better, such that
    it may be better to use the BaggingNeuralNet function.

    """

    def __init__(self,x,y,**kwargs):
        """
        Parameters
        ----------
        x : array{Float64}
            the spectra organised in rows (1 row = one spectrum). The spectra should share a common X axis.
        y : Array{Float64}
            Target. Only a single target is possible for now.

        """
        self.x = x
        self.y = y
        #
        # Kwargs extractions
        #
        self.X_test = kwargs.get("X_test",[0.0])
        self.y_test = kwargs.get("y_test",[0.0])
        self.algorithm = kwargs.get("algorithm","SVM")
        self.test_sz = kwargs.get("test_size",0.3)
        self.scaling = kwargs.get("scaling",True)
        self.scaler = kwargs.get("scaler","MinMaxScaler")
        self.rand_state = kwargs.get("rand_state",42)

        # hyperparameters for the algorithms
        self.user_kernel = kwargs.get("kernel","rbf")
        self.param_kr = kwargs.get(
            "param_kr",dict(alpha=[1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3],gamma=np.logspace(-4, 4, 9)))

        self.param_svm= kwargs.get(
            "param_svm",dict(C= [1e0, 2e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5], gamma= np.logspace(-4, 4, 9)))

        self.param_neurons = kwargs.get("param_neurons",
                                        dict(hidden_layer_sizes=(3,),
                                             solver='lbfgs',
                                             activation='relu',
                                             early_stopping=True,
                                             random_state=self.rand_state))

        self.param_bag = kwargs.get("param_bagging",
                                    dict(n_estimators=100, 
                                         max_samples=1.0, 
                                         max_features=1.0, 
                                         bootstrap=True, 
                                         bootstrap_features=False, 
                                         oob_score=False, 
                                         warm_start=False, 
                                         n_jobs=1, verbose=0,
                                         random_state=self.rand_state))    

        if len(self.X_test) == 1:
            self.X_train, self.X_test, self.y_train, self.y_test = sklearn.model_selection.train_test_split(
            self.x, self.y.reshape(-1, 1), test_size=self.test_sz, random_state=self.rand_state)
        elif self.X_test.shape[1] == self.x.shape[1]:
            self.X_train = np.copy(self.x)
            self.y_train = np.copy(self.y)
        else:
            ValueError("You tried to provide a testing dataset that has a different number of features (in columns) than the training set. Please correct this.")

        # to avoid any problem with SciKit Learn quite annoying demands for the shape of arrays...
        self.y_train = self.y_train.reshape(-1,1)
        self.y_test = self.y_test.reshape(-1, 1)

        # initialising the preprocessor scaler
        if self.scaler == "StandardScaler":
            self.X_scaler = sklearn.preprocessing.StandardScaler()
            self.Y_scaler = sklearn.preprocessing.StandardScaler()
        elif self.scaler == "MinMaxScaler":
            self.X_scaler = sklearn.preprocessing.MinMaxScaler()
            self.Y_scaler = sklearn.preprocessing.MinMaxScaler()
        else:
            InputError("Choose the scaler between MinMaxScaler and StandardScaler")

    def fit(self):
        """Train the model with the indicated algorithm.

        Do not forget to tune the hyperparameters.

        Parameters
        ----------
        algorithm : String,
            "KernelRidge", "SVM", "LinearRegression", "Lasso", "ElasticNet", "NeuralNet", "BaggingNeuralNet", default = "SVM"

        """
        self.X_scaler.fit(self.X_train)
        self.Y_scaler.fit(self.y_train)

        # scaling the data in all cases, it may not be used during the fit later
        self.X_train_sc = self.X_scaler.transform(self.X_train)
        self.y_train_sc = self.Y_scaler.transform(self.y_train)

        self.X_test_sc = self.X_scaler.transform(self.X_test)
        self.y_test_sc = self.Y_scaler.transform(self.y_test)

        if self.algorithm == "KernelRidge":
            clf_kr = KernelRidge(kernel=self.user_kernel)
            self.model = sklearn.model_selection.GridSearchCV(clf_kr, cv=5, param_grid=self.param_kr)

        elif self.algorithm == "SVM":
            clf_svm = SVR(kernel=self.user_kernel)
            self.model = sklearn.model_selection.GridSearchCV(clf_svm, cv=5, param_grid=self.param_svm)

        elif self.algorithm == "Lasso":
            clf_lasso = sklearn.linear_model.Lasso(alpha=0.1,random_state=self.rand_state)
            self.model = sklearn.model_selection.GridSearchCV(clf_lasso, cv=5,
                                                              param_grid=dict(alpha=np.logspace(-5,5,30)))

        elif self.algorithm == "ElasticNet":
            clf_ElasticNet = sklearn.linear_model.ElasticNet(alpha=0.1, l1_ratio=0.5,random_state=self.rand_state)
            self.model = sklearn.model_selection.GridSearchCV(clf_ElasticNet,cv=5, 
                                                              param_grid=dict(alpha=np.logspace(-5,5,30)))

        elif self.algorithm == "LinearRegression":
            self.model = sklearn.linear_model.LinearRegression()

        elif self.algorithm == "NeuralNet":
            self.model = MLPRegressor(**self.param_neurons)
        elif self.algorithm == "BaggingNeuralNet":
            nn_m = MLPRegressor(**self.param_neurons)

            self.model = BaggingRegressor(base_estimator = nn_m, **self.param_bag)

        if self.scaling == True:
            self.model.fit(self.X_train_sc, self.y_train_sc.reshape(-1,))
            predict_train_sc = self.model.predict(self.X_train_sc)
            self.prediction_train = self.Y_scaler.inverse_transform(predict_train_sc.reshape(-1,1))
            predict_test_sc = self.model.predict(self.X_test_sc)
            self.prediction_test = self.Y_scaler.inverse_transform(predict_test_sc.reshape(-1,1))
        else:
            self.model.fit(self.X_train, self.y_train.reshape(-1,))
            self.prediction_train = self.model.predict(self.X_train)
            self.prediction_test = self.model.predict(self.X_test)

    def predict(self,X):
        """Predict using the model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = (n_samples, n_features)
            Samples.

        Returns
        -------
        C : array, shape = (n_samples,)
            Returns predicted values.

        Remark
        ------
        if self.scaling == "yes", scaling will be performed on the input X.
        """
        if self.scaling == True:
            X_sc = self.X_scaler.transform(X)
            pred_sc = self.model.predict(X_sc)
            return self.Y_scaler.inverse_transform(pred_sc.reshape(-1,1))
        else:
            return self.model.predict(self.X)
