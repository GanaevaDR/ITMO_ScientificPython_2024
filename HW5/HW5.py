'''
IMPORT PART
'''

#import necessary libraries

! pip -q install rdkit

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.feature_selection import SelectKBest, mutual_info_regression
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer as Imputer
from sklearn.pipeline import Pipeline
import numpy as np
import pandas as pd
import urllib.request
import json

# copy the function of your choice from regression.py and the necessary imports for it -> it will perform hyperparameters tuning on selected regression model

def fit_ElasticNet(X_train, X_test, y_train, y_test):

    a = Imputer(missing_values=np.nan, strategy='median')
    b = StandardScaler()
    c = SelectKBest(score_func=mutual_info_regression)
    clf = ElasticNetCV(cv=10)
    model = Pipeline([('impute', a), ('scaling', b), ('anova', c), ('rf', clf)])

    # Grid Search CV
    parameters = {'anova__k': [5, 10, 20, 40]}

    grid = GridSearchCV(model, parameters)
    grid.fit(X_train, y_train)
    y_pred = grid.predict(X_test)

    # Metrics
    metric = [grid.score(X_test, y_test),
               metrics.explained_variance_score(y_test, y_pred),
               metrics.mean_absolute_error(y_test, y_pred),
               metrics.mean_squared_error(y_test, y_pred),
               metrics.median_absolute_error(y_test, y_pred),
               metrics.r2_score(y_test, y_pred)]

    return grid, y_pred, metric

# find and import the corresponding model from sklearn
from sklearn.linear_model import  ElasticNetCV

# find and import the function that computes all molecular descriptors "from .molecular_descriptors.py import ..."
from molecular_descriptors import getAllDescriptors


'''
DESCRIPTORS PART
'''

def desc_calc(data, mode) -> pd.DataFrame: # fill the necessary arguments
    return getAllDescriptors(data, mode)# compute all descriptors

def sar_model_evaluation(descriptors: pd.DataFrame):

    '''
    Function for model evaluation with fit_LinearSVR().
    Saves best model.
    '''

    y = descriptors.iloc[:, -1]
    X = descriptors.iloc[:, :-1] 

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    model1, y_pred1, metrics1 = fit_ElasticNet(X_train, X_test, y_train, y_test) # hyperparameters tuning with selected function
    return model1, y_pred1, metrics1

def sar_model_train(descriptors_train: pd.DataFrame, indices):

    '''
    Function for training the model with best paramaters.
    Don't forget to add here the input arguments required by the selected model.
    '''

    y_train = descriptors_train.iloc[:, -1] 
    X_train = descriptors_train.iloc[:, :-1]

    X_train = X_train[X_train.columns[indices]] # keeping only necessary descriptors according to ANOVA evaluation

    # reproducing the pipeline from GridSearchCV but for one selected model
    a = Imputer(missing_values=np.nan, strategy='median')
    b = StandardScaler()
    clf = ElasticNetCV() # selected model from sklearn
    model = Pipeline([('impute', a), ('scaling', b), ('rf', clf)]) # without ANOVA now

    model.fit(X_train, y_train)
    return model

def sar_model_predict(model, descriptors_pred, indices):

    '''
    Function for casting predictions on unseen data
    '''

    X_pred = descriptors_pred
    X_pred = X_pred[X_pred.columns[indices]]
    return model.predict(X_pred)


'''
PUBCHEM PART
'''

def pubchem_parsing(url):

    '''
    Function for Pubchem request
    '''

    req = urllib.request.Request(url)
    res = urllib.request.urlopen(req).read()
    fin = json.loads(res.decode())
    return fin


def get_similar_cids(compound_smiles, threshold=95, maxentries=10):

    '''
    Function for finding similar CIDS in PubChem with fastsimilarity_2d
    '''

    pubchem_pug_rest_api_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    pubchem_pug_rest_api_link+= "fastsimilarity_2d/smiles/%(smiles)s/cids/JSON?Threshold=%(threshold)s&MaxRecords=%(maxentries)s" % {
        "smiles": compound_smiles, "threshold": threshold, "maxentries": maxentries}
    similar_cids = pubchem_parsing(pubchem_pug_rest_api_link)['IdentifierList']['CID']
    return similar_cids


def get_xlogp(compound_cid):

    '''
    Function for parsing XLogP from the response
    '''

    pubchem_pug_rest_api_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    pubchem_pug_rest_api_link += "compound/cid/%s/property/XLogP/JSON" % compound_cid

    try:
        xlogp = pubchem_parsing(pubchem_pug_rest_api_link)['PropertyTable']['Properties'][0]['XLogP']
        return xlogp
    except KeyError:
        return None


'''
MAIN PART
'''

if __name__ == "__main__":

    pd.set_option.use_inf_as_na = True

    # loading data
    train_data = pd.read_csv('logp_100.csv')
    pred_data = pd.read_csv('logp_inputs.csv')
    cpds = [row for row in pred_data.loc[:, 'SMILES']]

    # calculating descriptors
    print("Calculating descriptors for training data...")
    train_descriptors = desc_calc(train_data, 'train')
    print("Calculating descriptors for prediction data...")
    pred_descriptors = desc_calc(pred_data, 'pred')

    # finding best estimator
    print("Evaluating regression model parameters...")
    model = sar_model_evaluation(train_descriptors)
    print('Best parameters are:', model[0].best_params_)
    cols = model[0].best_estimator_.named_steps['anova'].get_support(indices=True) # this are indices from ANOVA

    # train the best estimator and predict values
    print("Training the model with the best parameters...")
    final_model = sar_model_train(train_descriptors, cols)

    for cpd in cpds:
        cpd_descriptors = pred_descriptors[pred_descriptors['SMILES']==cpd]
        pred = sar_model_predict(final_model, cpd_descriptors, cols)
        print(f"Predicted LogP value for compound {cpd}:", pred)

        result = []

        print("Searching for similar compounds...")
        similarity = get_similar_cids(cpd) # related pubchem function

        print("Filtering logP...")
        for cid in similarity:
            xlogp =  get_xlogp(cid) # related pubchem function
            if xlogp:
                if xlogp <= pred*1.1 and xlogp >=pred*0.9:
                    result.append((cid, xlogp))

        print(f"Request for compound {cpd} completed. I found the following CIDs in PubChem with XLogP in the range of {pred}+- 10%: {result}")
