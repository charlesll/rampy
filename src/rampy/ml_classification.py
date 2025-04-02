# -*- coding: utf-8 -*-
#############################################################################
#Copyright (c) 2018-2025 Charles Le Losq
#
# Licence GNU-GPL
#
#
#############################################################################
import numpy as np

import sklearn

class mlclassificator:
    """Perform automatic classification of spectral data using scikit-learn machine learning algorithms.

    This class supports various classification algorithms and allows customization of hyperparameters. 
    It also handles scaling and splitting of training and testing datasets.

    Attributes:
        x (np.ndarray): Training spectra organized in rows (1 row = one spectrum).
        y (np.ndarray): Target labels for training data.
        X_test (np.ndarray): Testing spectra organized in rows.
        y_test (np.ndarray): Target labels for testing data.
        algorithm (str): Machine learning algorithm to use. Options:
            "Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
            "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
            "Naive Bayes", "QDA".
        scaling (bool): Whether to scale the data during fitting and prediction.
        scaler (str): Type of scaler to use ("MinMaxScaler" or "StandardScaler").
        test_size (float): Fraction of the dataset to use as a testing dataset if X_test and y_test are not provided.
        rand_state (int): Random seed for reproducibility. Default is 42.
        params_ (dict): Hyperparameters for the selected algorithm.
        model: Scikit-learn model instance.
        X_scaler: Scikit-learn scaler instance for X values.
    """

    def __init__(self,x,y,**kwargs):
        """Initialize the mlclassificator.

        Args:
            x (np.ndarray): Training spectra organized in rows.
            y (np.ndarray): Target labels for training data.

        Keyword Args:
            X_test (np.ndarray, optional): Testing spectra. Default is None.
            y_test (np.ndarray, optional): Testing labels. Default is None.
            algorithm (str, optional): Machine learning algorithm to use. Default is "Nearest Neighbors".
            test_size (float, optional): Fraction of data used for testing if X_test and y_test are not provided. Default is 0.3.
            scaling (bool, optional): Whether to scale data. Default is True.
            scaler (str, optional): Type of scaler ("MinMaxScaler" or "StandardScaler"). Default is "MinMaxScaler".
            rand_state (int, optional): Random seed for reproducibility. Default is 42.
            params_ (dict, optional): Hyperparameters for the selected algorithm. Default is None.

        Raises:
            ValueError: If X_test has a different number of features than x.
        """
        self.x = x
        self.y = y
        self.X_test = kwargs.get("X_test", None)
        self.y_test = kwargs.get("y_test", None)
        self.algorithm = kwargs.get("algorithm", "Nearest Neighbors")
        self.test_size = kwargs.get("test_size", 0.3)
        self.scaling = kwargs.get("scaling", True)
        self.scaler = kwargs.get("scaler", "MinMaxScaler")
        self.rand_state = kwargs.get("rand_state", 42)
        self.params_ = kwargs.get("params_", None)

        # Split data if no test set is provided
        if self.X_test is None or self.y_test is None:
            from sklearn.model_selection import train_test_split
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
                self.x, self.y, test_size=self.test_size, random_state=self.rand_state
            )
        else:
            if self.X_test.shape[1] != self.x.shape[1]:
                raise ValueError("X_test must have the same number of features as x.")
            self.X_train = np.copy(self.x)
            self.y_train = np.copy(self.y)

        # Initialize scaler
        if self.scaler == "StandardScaler":
            from sklearn.preprocessing import StandardScaler
            self.X_scaler = StandardScaler()
        elif self.scaler == "MinMaxScaler":
            from sklearn.preprocessing import MinMaxScaler
            self.X_scaler = MinMaxScaler()
        else:
            raise ValueError("Invalid scaler type. Choose 'MinMaxScaler' or 'StandardScaler'.")
        
        # Define dispatcher for algorithms
        from sklearn.neighbors import KNeighborsClassifier
        from sklearn.svm import SVC
        from sklearn.gaussian_process import GaussianProcessClassifier
        from sklearn.tree import DecisionTreeClassifier
        from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
        from sklearn.naive_bayes import GaussianNB
        from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
        from sklearn.neural_network import MLPClassifier

        self.dispatcher = {
            "Nearest Neighbors": KNeighborsClassifier,
            "SVC": SVC,
            "Gaussian Process": GaussianProcessClassifier,
            "Decision Tree": DecisionTreeClassifier,
            "Random Forest": RandomForestClassifier,
            "Neural Net": MLPClassifier,
            "AdaBoost": AdaBoostClassifier,
            "Naive Bayes": GaussianNB,
            "QDA": QuadraticDiscriminantAnalysis,
        }

    def scale_data(self):
        """Scale training and testing data."""
        if not hasattr(self, 'X_scaler'):
            raise AttributeError("Scaler has not been initialized.")
        
        # Fit and transform training data; transform testing data
        self.X_train_sc = self.X_scaler.fit_transform(self.X_train)
        self.X_test_sc = self.X_scaler.transform(self.X_test)

    def fit(self, params_: dict = None):
        """Scale data and train or re-train the model with the specified algorithm.

        This method initializes and trains the model if it hasn't been trained yet. If a model
        already exists (from a previous fit), it reuses the existing model and optionally updates
        its hyperparameters.

        Args:
            params_ (dict, optional): Hyperparameters for the selected algorithm. If provided,
                these parameters will override any previously set parameters.

        Raises:
            ValueError: If an invalid algorithm is specified or if scaling is inconsistent.
        """
        # Update hyperparameters if provided
        if params_ is not None:
            self.params_ = params_

        # Ensure params_ is always a valid dictionary
        self.params_ = self.params_ or {}

        # Initialize or reuse the model
        ModelClass = self.dispatcher[self.algorithm]
        self.model = ModelClass(**self.params_)

        # Scale data if required
        if self.scaling:
            self.X_scaler.fit(self.X_train)
            self.X_train_sc = self.X_scaler.transform(self.X_train)
            self.X_test_sc = self.X_scaler.transform(self.X_test)
            X_train_used = self.X_train_sc
            X_test_used = self.X_test_sc
        else:
            X_train_used = self.X_train
            X_test_used = self.X_test

        # Fit the model on training data and make predictions
        self.model.fit(X_train_used, self.y_train.reshape(-1,))
        self.prediction_train = self.model.predict(X_train_used)
        self.prediction_test = self.model.predict(X_test_used)


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
        """Predict target values using the trained model.

        Args:
            X (np.ndarray): Samples to predict with shape `(n_samples, n_features)`.

        Returns:
            np.ndarray: Predicted target values with shape `(n_samples,)`.

        Notes:
            - If `scaling` is enabled, input samples will be scaled before prediction.

        Raises:
           ValueError: If the model has not been fitted yet.
        """
        if self.scaling == True:
            X_sc = self.X_scaler.transform(X)
            return self.model.predict(X_sc)
        else:
            return self.model.predict(X)
