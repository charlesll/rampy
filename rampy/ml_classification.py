import numpy as np

import sklearn
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import BaggingClassifier

class mlclassificator:
    """use machine learning algorithms from scikit learn to perform classification of spectra.
    
    Attributes
    ----------
    x : {array-like, sparse matrix}, shape = (n_samples, n_features)
        Spectra; n_features = n_frequencies.
    y : array, shape = (n_samples,)
        numeric labels.
    X_test : {array-like, sparse matrix}, shape = (n_samples, n_features)
        spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
    y_test : array, shape = (n_samples,)
        numeric labels that you want to use as a testing dataset. Those targets should not be present in the y (training) dataset.
    algorithm : String,
        "Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
         "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
         "Naive Bayes", "QDA"  
    scaling : Bool
        True or False. If True, data will be scaled during fitting and prediction with the requested scaler (see below),
    scaler : String
        the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see http://scikit-learn.org/stable/modules/preprocessing.html for details. Default = "MinMaxScaler".
    test_size : float
        the fraction of the dataset to use as a testing dataset; only used if X_test and y_test are not provided.
    rand_state : Float64
        the random seed that is used for reproductibility of the results. Default = 42.
    params_ : Dictionary
        contain the values of the hyperparameters that should be provided to the algorithm. See scikit-learn documentation for details for each algorithm.
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
    
    Example
    -------
    Given an array X of n samples by m frequencies, and Y an array of n x 1 concentrations
    >>> model = rampy.mlclassificator(X,y)
    >>> model.algorithm("SVC")
    >>> model.user_kernel = 'poly'
    >>> model.fit()
    >>> y_new = model.predict(X_new)
    
    Remarks
    -------
    For details on hyperparameters of each algorithms, please directly consult the documentation of SciKit Learn at:
    http://scikit-learn.org/stable/
    
    In progress
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
        self.algorithm = kwargs.get("algorithm","Nearest Neighbors")
        self.test_sz = kwargs.get("test_size",0.3)
        self.scaling = kwargs.get("scaling",True)
        self.scaler = kwargs.get("scaler","MinMaxScaler")
        self.rand_state = kwargs.get("rand_state",42)

        # hyperparameters for the algorithms
        self.user_kernel = kwargs.get("kernel","rbf")
        
        self.params_ = kwargs.get(
            "params_",None)

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
        elif self.scaler == "MinMaxScaler":
            self.X_scaler = sklearn.preprocessing.MinMaxScaler()
        else:
            InputError("Choose the scaler between MinMaxScaler and StandardScaler")
            
            
        # now defining the model functions in a safe way:
        self.dispatcher = {"Nearest Neighbors" : KNeighborsClassifier(3), 
                      "Linear SVM" : SVC(kernel="linear", C=0.025),
                      "RBF SVM" : SVC(gamma=2, C=1),
                      "Gaussian Process" : GaussianProcessClassifier(),
                      "Decision Tree" : DecisionTreeClassifier(max_depth=15), 
                      "Random Forest" : RandomForestClassifier(max_depth=15, n_estimators=5, max_features=2), 
                      "Neural Net": MLPClassifier(),
                      "AdaBoost": AdaBoostClassifier(),
                      "Naive Bayes": GaussianNB(), 
                      "QDA": QuadraticDiscriminantAnalysis()}
                      
    def fit(self):
        """Scale data and train the model with the indicated algorithm.
        Do not forget to tune the hyperparameters.
        Parameters
        ----------
        algorithm : String,
            "Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
         "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
         "Naive Bayes", "QDA"
        """
        
        # scaling the data in all cases, it may not be used during the fit later
        self.X_scaler.fit(self.X_train)
        self.X_train_sc = self.X_scaler.transform(self.X_train)
        self.X_test_sc = self.X_scaler.transform(self.X_test)

        self.model = self.dispatcher[self.algorithm]
        
        if self.params_ != None:
            self.model(**self.params)
            #self.model = BaggingRegressor(base_estimator = nn_m, **self.param_bag)

        if self.scaling == True:
            self.model.fit(self.X_train_sc, self.y_train.reshape(-1,))
            self.prediction_train = self.model.predict(self.X_train_sc)
            self.prediction_test = self.model.predict(self.X_test_sc)
        else:
            self.model.fit(self.X_train, self.y_train.reshape(-1,))
            self.prediction_train = self.model.predict(self.X_train)
            self.prediction_test = self.model.predict(self.X_test)

    def refit(self):
        """Re-train a model previously trained with fit()
        """
        if self.scaling == True:
            self.model.fit(self.X_train_sc, self.y_train.reshape(-1,))
            self.prediction_train = self.model.predict(self.X_train_sc)
            self.prediction_test = self.model.predict(self.X_test_sc)
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
            return self.model.predict(X_sc)
        else:
            return self.model.predict(self.X)