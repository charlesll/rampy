import sklearn
from sklearn.decomposition import PCA, NMF

class mlexplorer:
    """use machine learning algorithms from scikit learn to explore spectroscopic datasets

    Performs automatic scaling and train/test split before NMF or PCA fit.

    Attributes
    ----------
    x : {array-like, sparse matrix}, shape = (n_samples, n_features)
        Spectra; n_features = n_frequencies.
    X_test : {array-like, sparse matrix}, shape = (n_samples, n_features)
        spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
    algorithm : String,
        "PCA", "NMF", default = "PCA"
    scaling : Bool
        True or False. If True, data will be scaled prior to fitting (see below),
    scaler : String
        the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see http://scikit-learn.org/stable/modules/preprocessing.html for details. Default = "MinMaxScaler".
    test_size : float
        the fraction of the dataset to use as a testing dataset; only used if X_test and y_test are not provided.
    rand_state : Float64
        the random seed that is used for reproductibility of the results. Default = 42.
    model : Scikit learn model
        A Scikit Learn object model, see scikit learn library documentation.

    Remarks
    -------

    For details on hyperparameters of each algorithms, please directly consult the documentation of SciKit Learn at:

    http://scikit-learn.org/stable/

    Results for machine learning algorithms can vary from run to run. A way to solve that is to fix the random_state.

    Example
    -------

    Given an array X of n samples by m frequencies, and Y an array of n x 1 concentrations

    >>> explo = rampy.mlexplorer(X) # X is an array of signals built by mixing two partial components
    >>> explo.algorithm = 'NMF' # using Non-Negative Matrix factorization
    >>> explo.nb_compo = 2 # number of components to use
    >>> explo.test_size = 0.3 # size of test set
    >>> explo.scaler = "MinMax" # scaler
    >>> explo.fit() # fitting!
    >>> W = explo.model.transform(explo.X_train_sc) # getting the mixture array
    >>> H = explo.X_scaler.inverse_transform(explo.model.components_) # components in the original space
    >>> plt.plot(X,H.T) # plot the two components

    """

    def __init__(self,x,**kwargs):
        """
        Parameters
        ----------
        x : array{Float64}
            the spectra organised in rows (1 row = one spectrum). The spectra should share a common X axis.

        """
        self.x = x
        #
        # Kwargs extractions
        #
        self.X_test = kwargs.get("X_test",[0.0])
        self.algorithm = kwargs.get("algorithm","PCA")
        self.test_size = kwargs.get("test_size",0.3)
        self.scaling = kwargs.get("scaling",True)
        self.scaler = kwargs.get("scaler","MinMaxScaler")
        self.rand_state = kwargs.get("rand_state",42)
        self.nb_compo = kwargs.get("n_components",2)

        if len(self.X_test) == 1:
            self.X_train, self.X_test = sklearn.model_selection.train_test_split(
            self.x, test_size=self.test_size, random_state=self.rand_state)
        elif self.X_test.shape[1] == self.x.shape[1]:
            self.X_train = np.copy(self.x)
        else:
            ValueError("You tried to provide a testing dataset that has a different number of features (in columns) than the training set. Please correct this.")

        # initialising the preprocessor scaler
        if self.scaler == "StandardScaler":
            self.X_scaler = sklearn.preprocessing.StandardScaler()
        elif self.scaler == "MinMaxScaler":
            self.X_scaler = sklearn.preprocessing.MinMaxScaler()
        else:
            InputError("Choose the scaler between MinMaxScaler and StandardScaler")

        # fitting scaler
        self.X_scaler.fit(self.X_train)

        # scaling the data in all cases, it may not be used during the fit later
        self.X_train_sc = self.X_scaler.transform(self.X_train)
        self.X_test_sc = self.X_scaler.transform(self.X_test)

    def fit(self):
        """Train the model with the indicated algorithm.

        Do not forget to tune the hyperparameters.

        """
        if self.algorithm == "PCA":
            self.model = PCA(n_components=self.nb_compo)
        elif self.algorithm == "NMF":
            self.model = NMF(n_components=self.nb_compo,init = "nndsvd")

        if self.scaling == True:
            self.model.fit(self.X_train_sc)
        else:
            self.model.fit(self.X_train)

    def refit(self):
        """Train the model with the indicated algorithm.

        Do not forget to tune the hyperparameters.

        """
        if self.scaling == True:
            self.model.fit(self.X_train_sc)
        else:
            self.model.fit(self.X_train)

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
