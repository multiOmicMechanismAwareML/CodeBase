import csv
import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
import numpy as np
import pandas as pa
import os.path
import matplotlib.pyplot as plts
import seaborn as sns
from keras.models import Sequential, Model
from keras.layers.core import Dense, Dropout, Activation
from keras.layers import Concatenate, Input
from keras.optimizers import SGD
from keras.constraints import max_norm
from sklearn import preprocessing
from keras.utils import plot_model

def garson(A, B):
    """
    Computes Garson's algorithm
    A = matrix of weights of input-hidden layer (rows=input & cols=hidden)
    B = vector of weights of hidden-output layer
    """
    # B = np.diag(B)
    # connection weight through the different hidden node
    cw = np.dot(A, B)
    # weight through node (axis=0 is column; sum per input feature)
    cw_h = abs(cw).sum(axis= 0)

    # relative contribution of input neuron to outgoing signal of each hidden neuron
    # sum to find relative contribution of input neuron
    rc = np.divide(abs(cw), abs(cw_h))
    rc = rc.sum(axis=1)

    # normalize to 100% for relative importance
    ri = rc / rc.sum()
    return(ri)

#Checks that the reading of the data has worked
def check_file_read_ok(data, name):
    if data is None:
        print("error with " + name + " data not read")
    else:
        print("Success in loading " + name)
        print(data.shape)

def remove_zero_entry_columns(data):
    return data.loc[:, (data != 0).any(axis = 0)]

#Name of the row we are trying to predict (growth rate)
TARGET_NAME = 'log2relT'


with open('testing_index.csv', 'r') as csvfile:
    testing_index = []
    for row in csv.reader(csvfile, delimiter=';'):
        testing_index.append(row[0]) # careful here with [0]
#Load the data removing any
full_data = pa.read_csv('data/completeDataset.csv')
check_file_read_ok(full_data, "full data")

# Here we extract the features highlighted by the genetic algorithm
def get_the_genetic_algo_data():
    gens = []
    features = pa.read_csv("data/genetic_feature_selection_features.csv")
    for x in range(features.shape[0]):
        fet = features.iloc[[0]].values[1:]
        gens.append(full_data[fet])
    return gens

gens = get_the_genetic_algo_data()

full_data = remove_zero_entry_columns(full_data)
expression_data = pa.read_csv('data/expressionOnly.csv')
check_file_read_ok(expression_data, "expression data")
expression_data = remove_zero_entry_columns(expression_data)

metabolic_expression_data = pa.read_csv('data/metabolic_gene_data.csv')
check_file_read_ok(metabolic_expression_data, "met_expression data")

#Extract the target and drop target column from main data
target_data = full_data[TARGET_NAME]
full_data = full_data.drop(columns=TARGET_NAME)
expression_data = expression_data.drop(columns=TARGET_NAME)
#metabolic_expression = metabolic_expression.drop(columns=TARGET_NAME)

iRF = pa.read_csv('data/Features_Extracted_Using_iRF.csv', header = None)
iRF.columns = ['Genes']
check_file_read_ok(iRF, "iRF data")
iRF = remove_zero_entry_columns(iRF)
iRF = full_data[iRF.Genes]
print(iRF.shape, full_data.shape)


sgl = pa.read_csv('data/Features_Extracted_Using_SGL.csv')
sgl.columns = ['Genes']
check_file_read_ok(sgl, "sgl data")
sgl_list= []
for x in sgl.Genes:
    if x in full_data.columns:
        sgl_list.append(x)   #We do this because some of the zero features have been removed from full_data

sgl = full_data[sgl_list]
genes = full_data['Row']


full_data = full_data.drop(columns = 'Row')
reaction_data = full_data.drop(columns = expression_data.columns.values)

def init_model_3layer(input_dim, learning_rate, epochs, momentum, neurons):
    model = Sequential()
    model.add(Dense(neurons, input_dim = input_dim, kernel_constraint=max_norm(3)))
    model.add(Activation('sigmoid'))
    model.add(Dropout(0.6))
    model.add(Dense(neurons, input_dim = input_dim, kernel_constraint=max_norm(3)))
    model.add(Activation('sigmoid'))
    model.add(Dropout(0.6))
    model.add(Dense(neurons, input_dim = input_dim, kernel_constraint=max_norm(3)))
    model.add(Activation('sigmoid'))
    model.add(Dropout(0.6))
    model.add(Dense(1))
    model.add(Activation('linear'))
    rms = SGD(lr= learning_rate, decay= learning_rate / epochs, momentum=momentum)
    model.compile(loss='mean_absolute_error', optimizer=rms, metrics = ["mean_absolute_error", 'mean_squared_error'])
    return model

def init_model(input_dim, learning_rate, epochs, momentum, neurons, trainable = True):
    input = Input(shape = (input_dim,))
    layer = Dense(neurons, activation='sigmoid', kernel_constraint=max_norm(3), name = "expression_1") (input)
    layer = Dropout(rate=0.6) (layer)
    layer = Dense(neurons, activation='sigmoid', kernel_constraint=max_norm(3), name = "expression_2") (layer)
    layer = Dropout(rate=0.6) (layer)
    predictions = Dense(1, activation='linear') (layer)
    model = Model(inputs = input, outputs = predictions)
    rms = SGD(lr= learning_rate, decay= learning_rate / epochs, momentum=momentum)
    model.trainable = trainable
    if (trainable) :
        model.compile(loss='mean_squared_error', optimizer=rms, metrics = ["mean_absolute_error"])
    return model

def init_multi_model(input_dim,input_dim2, learning_rate, epochs, momentum, neurons, reaction_trained, expression_trained):
    reaction_input = Input(shape = (input_dim,))
    expression_input = Input(shape = (input_dim2,))
    comb_layer = Concatenate()([reaction_trained(reaction_input), expression_trained(expression_input)])
    comb_layer = Dense(neurons, activation='sigmoid', kernel_constraint=max_norm(3), name = "last_hidden") (comb_layer)
    predictions = Dense(1, activation='linear') (comb_layer)
    model = Model(inputs = [reaction_input,expression_input], outputs = predictions)
    rms = SGD(lr= learning_rate, decay= learning_rate / epochs, momentum=momentum)
    model.compile(loss='mean_squared_error', optimizer=rms, metrics = ["mean_absolute_error"])
    return model

SEED = 120
number_of_instances = len(target_data)
testing_index = list(map(int, testing_index))
training_index = np.setxor1d(range(1,number_of_instances), testing_index)
epochs = 6000
batches = 256
lrate = 0.005
validation = 0.1
print(training_index)
#Split the data 80:20
full_data_train, full_data_test = full_data.drop(full_data.index[testing_index]), full_data.iloc[testing_index, :]
expression_data_train, expression_data_test = expression_data.drop(expression_data.index[testing_index]), expression_data.iloc[testing_index, :]
metabolic_expression_data_train, metabolic_expression_data_test = metabolic_expression_data.drop(metabolic_expression_data.index[testing_index]), metabolic_expression_data.iloc[testing_index, :]
#metabolic_expression_data_train, metabolic_expression_data_test = metabolic_expression_data[training_index,:], metabolic_expression_data[testing_index,:]
reaction_data_train, reaction_data_test = reaction_data.drop(reaction_data.index[testing_index]), reaction_data.iloc[testing_index, :]
target_data_train, target_data_test = target_data.drop(target_data.index[testing_index]), target_data.iloc[testing_index]
iRF_train, iRF_test = iRF.drop(iRF.index[testing_index]), iRF.iloc[testing_index, :]
sgl_train, sgl_test = sgl.drop(sgl.index[testing_index]), sgl.iloc[testing_index, :]
print(sgl_train.shape, 'SGL')
#Preprocessing the data to mean of zero and unit variance - stopped using this for improved results
full_scaler = preprocessing.StandardScaler().fit(full_data_train)
full_data_scaled_train = full_scaler.transform(full_data_train).astype(np.float32)
full_data_scaled_test = full_scaler.transform(full_data_test).astype(np.float32)

expression_scaler = preprocessing.StandardScaler().fit(expression_data_train)
expression_data_scaled_train = expression_scaler.transform(expression_data_train).astype(np.float32)
expression_data_scaled_test = expression_scaler.transform(expression_data_test).astype(np.float32)

reaction_scaler = preprocessing.StandardScaler().fit(reaction_data_train)
reaction_data_scaled_train = reaction_scaler.transform(reaction_data_train).astype(np.float32)
reaction_data_scaled_test = reaction_scaler.transform(reaction_data_test).astype(np.float32)

metabolic_expression_scaler = preprocessing.StandardScaler().fit(metabolic_expression_data_train)
metabolic_expression_data_scaled_train = metabolic_expression_scaler.transform(metabolic_expression_data_train).astype(np.float32)
metabolic_expression_data_scaled_test = metabolic_expression_scaler.transform(metabolic_expression_data_test).astype(np.float32)

iRF_scaler = preprocessing.StandardScaler().fit(iRF_train)
iRF_scaled_train = iRF_scaler.transform(iRF_train).astype(np.float32)
iRF_scaled_test = iRF_scaler.transform(iRF_test).astype(np.float32)

sgl_scale = preprocessing.StandardScaler().fit(sgl_train)
sgl_scaled_train = sgl_scale.transform(sgl_train).astype(np.float32)
sgl_scaled_test = sgl_scale.transform(sgl_test).astype(np.float32)

target_train = target_data_train.astype(np.float32)
target_test = target_data_test.astype(np.float32)

result = {}
earlyStopping=EarlyStopping(monitor='val_loss', patience=15000, verbose=0, mode='auto')
print("STARTING")

model_metabolic_expression = init_model(metabolic_expression_data_scaled_train.shape[1], lrate, 3000, 0.75, 1000)
if os.path.exists('models/metabolic_expression_model.h5'):
    model_metabolic_expression.load_weights('models/metabolic_expression_model.h5')
else:
    model_metabolic_expression.fit(x = metabolic_expression_data_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation, callbacks=[earlyStopping])
    model_metabolic_expression.save_weights('models/metabolic_expression_model.h5')
score = model_metabolic_expression.evaluate(metabolic_expression_data_scaled_test, target_test, verbose=1)
print("GEM score" , score)
if not os.path.exists('predictions/GEM_DL_Predictions.csv'):
    prediction = model_metabolic_expression.predict_on_batch(metabolic_expression_data_scaled_test)
    np.savetxt(fname="predictions/GEM_DL_Predictions.csv",  X=prediction, delimiter=',')

model_full = init_model(full_data_scaled_train.shape[1], lrate,3000, 0.75, 1000)
if os.path.exists("models/concat_Flu_GE.h5"):
    model_full.load_weights("models/concat_Flu_GE.h5")
else:
    model_full.fit(x = full_data_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation, callbacks=[earlyStopping])
    model_full.save_weights("models/concat_Flu_GE.h5")
score = model_full.evaluate(full_data_scaled_test, target_test, verbose=1)
print("concate_Flu_GE score ", score)
if not os.path.exists("predictions/concate_Flu_GE_DL_predictions.csv"):
    prediction = model_full.predict_on_batch(full_data_scaled_test)
    np.savetxt(fname="predictions/concate_Flu_GE_DL_predictions.csv",  X=prediction, delimiter=',')

model_expression = init_model(expression_data_train.shape[1], lrate,3000, 0.75, 1000)
if os.path.exists('models/expression_model.h5'):
    model_expression.load_weights('models/expression_model.h5')
else:
    model_expression.fit(x = expression_data_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation, callbacks=[earlyStopping])
    model_expression.save_weights('models/expression_model.h5')
score = model_expression.evaluate(expression_data_scaled_test, target_test, verbose=1)
print("GE score ", score)
if not os.path.exists("predictions/GE_DL_predictions"):
    predictions = model_expression.predict(expression_data_scaled_test)
    np.savetxt(fname="predictions/GE_DL_predictions.csv", X=predictions, delimiter=',')

model_reaction = init_model(reaction_data_scaled_train.shape[1], lrate, 3000, 0.75,1000)
if os.path.exists("models/reaction_model.h5"):
    model_reaction.load_weights("models/reaction_model.h5")
else:
    model_reaction.fit(x = reaction_data_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation, callbacks=[earlyStopping])
    model_reaction.save_weights("models/reaction_model.h5")
score = model_reaction.evaluate(reaction_data_scaled_test, target_test, verbose=1)
print("Flu score " , score)
if not os.path.exists("predictions/Flu_DL_predictions.csv"):
    predictions = model_reaction.predict(reaction_data_scaled_test)
    np.savetxt(fname="predictions/Flu_DL_predictions.csv", X=predictions, delimiter=',')

if not os.path.exists("models/NSGA-II_model.h5"):
    min_score = [9999999,0]
    best_fs_model = None
    for x in range(9):
        next = gens[x]
        print(next.shape)
        next_train, next_test = next.drop(next.index[testing_index]), next.iloc[testing_index, :]
        scale = preprocessing.StandardScaler().fit(next_train)
        next_scale_train = scale.transform(next_train).astype(np.float32)
        next_scale_test = scale.transform(next_test).astype(np.float32)
        model = init_model(next_scale_train.shape[1], lrate, 3000, 0.75,1000)
        model.fit(x = next_scale_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation)
        score = model.evaluate(next_scale_test, target_test, verbose=1)
        if(score[0] < min_score[0]):
            min_score = score
            best_fs_model = model
            predictions = model.predict(next_scale_test)
    print("NSGA-II score ", min_score)
    np.savetxt(fname = "predictions/NSGA-II_DL_Predictions.csv", X = predictions, delimiter=',')
    best_fs_model.save_weights("models/NSGA-II_model.h5")
else:
    print("we have the NSGA-II scores")

model_iRF = init_model(iRF_scaled_train.shape[1], lrate, 3000, 0.75, 1000)
if not os.path.exists("models/iRF_model.h5"):
        model_iRF.fit(x = iRF_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation)
        model_iRF.save_weights("models/iRF_model.h5")
else :
    model_iRF.load_weights("models/iRF_model.h5")
score = model_iRF.evaluate(iRF_scaled_test, target_test, verbose=1)
print("iRF score ", score)
if not os.path.exists("predictions/iRF_DL_predictions.csv"):
    predictions  =model_iRF.predict(iRF_scaled_test)
    np.savetxt(fname="predictions/iRF_DL_predictions.csv", X = predictions, delimiter=',')

model_SGL = init_model(sgl_scaled_train.shape[1], lrate, 3000, 0.75, 1000)
if not os.path.exists("models/SGL_model.h5"):
        print('IN HERE')
        model_SGL.fit(x = sgl_scaled_train, y = target_train, epochs=epochs, batch_size=batches, validation_split=validation)
        model_SGL.save_weights("models/SGL_model.h5")
else :
    model_SGL.load_weights("models/SGL_model.h5")
score = model_SGL.evaluate(sgl_scaled_test, target_test, verbose=1)
print("SGL score ", score)
if not os.path.exists("predictions/SGL_DL_predictions.csv"):
    predictions  =model_SGL.predict(sgl_scaled_test)
    np.savetxt(fname="predictions/SGL_DL_predictions.csv", X = predictions, delimiter=',')

#Next we remove the  last layers of the pretrained to build the multi_model models

lrate = 0.05
epochs = 500
model_metabolic_expression = init_model(metabolic_expression_data_scaled_train.shape[1], lrate, 3000, 0.75, 1000, False)
model_expression = init_model(expression_data_train.shape[1], lrate, 3000, 0.75, 1000, False)
model_reaction = init_model(reaction_data_scaled_train.shape[1], lrate, 3000, 0.75,1000, False)
rms = SGD(lr= lrate , decay= lrate / epochs, momentum=0.75)
model_metabolic_expression.load_weights('models/metabolic_expression_model.h5')
model_metabolic_expression.trainable = True
model_metabolic_expression.compile(loss='mean_squared_error', optimizer=rms, metrics = ["mean_absolute_error"])
model_reaction.load_weights("models/reaction_model.h5")
model_reaction.trainable = True
model_reaction.compile(loss='mean_squared_error', optimizer=rms, metrics = ["mean_absolute_error"])
model_expression.trainable = True
model_expression.load_weights('models/expression_model.h5')
model_expression.compile(loss='mean_squared_error', optimizer=rms, metrics = ["mean_absolute_error"])

model_expression.layers.pop()
model_expression.layers.pop()
model_expression.outputs = [model_expression.layers[-1].output]
model_expression.layers[-1].outbound_nodes = []

model_metabolic_expression.layers.pop()
model_metabolic_expression.layers.pop()
model_metabolic_expression.outputs = [model_metabolic_expression.layers[-1].output]
model_metabolic_expression.layers[-1].outbound_nodes = []

model_reaction.layers.pop()
model_reaction.layers.pop()
model_reaction.outputs = [model_reaction.layers[-1].output]
model_reaction.layers[-1].outbound_nodes = []


multi_model_full_expression = init_multi_model(reaction_data_scaled_train.shape[1], expression_data_scaled_train.shape[1], lrate, epochs, 0.75, 15, model_reaction, model_expression)
multi_model_metabolic_expression = init_multi_model(reaction_data_scaled_train.shape[1], metabolic_expression_data_scaled_train.shape[1], lrate, epochs, 0.75, 10, model_reaction, model_metabolic_expression)


if not os.path.exists("models/MM-Flu_GE.h5"):
        multi_model_full_expression.fit(x = [reaction_data_scaled_train, expression_data_scaled_train],
         y = target_train,
         epochs=epochs,
         batch_size=batches,
         validation_split=validation)
        multi_model_full_expression.save_weights("models/MM-Flu_GE.h5")
else :
    multi_model_full_expression.load_weights("models/MM-Flu_GE.h5")
score = multi_model_full_expression.evaluate([reaction_data_scaled_test, expression_data_scaled_test], target_test, verbose=1)
print("MM-Full-Expression score ", score)
if not os.path.exists("predictions/MM-GE-Flu_DL_Predictions.csv"):
    predictions  =multi_model_full_expression.predict([reaction_data_scaled_test, expression_data_scaled_test])
    np.savetxt(fname="predictions/MM-GE-Flu_DL_Predictions.csv", X = predictions, delimiter=',')



if not os.path.exists("models/MM-Flu_GEM.h5"):
        multi_model_metabolic_expression.fit(x = [reaction_data_scaled_train, metabolic_expression_data_scaled_train],
         y = target_train,
         epochs=epochs,
         batch_size=batches,
         validation_split=validation)
        multi_model_metabolic_expression.save_weights("models/MM-Flu_GEM.h5")
else :
    multi_model_metabolic_expression.load_weights("models/MM-Flu_GEM.h5")
score = multi_model_metabolic_expression.evaluate([reaction_data_scaled_test, metabolic_expression_data_scaled_test], target_test, verbose=1)

print("MM-Meta-Expression score ", score)
if not os.path.exists("predictions/MM-Flu_GEM_DL_Predictions.csv"):
    predictions  =multi_model_metabolic_expression.predict([reaction_data_scaled_test, metabolic_expression_data_scaled_test])
    np.savetxt(fname="predictions/MM-Flu_GEM_DL_Predictions.csv", X = predictions, delimiter=',')


if not os.path.exists("layer_examine_full.pdf"):
    plts.figure()
    weights = multi_model_full_expression.get_layer(name = "last_hidden").get_weights()[0]
    output = multi_model_full_expression.layers[-1].get_weights()
    np.set_printoptions(threshold=20000)
    output = np.asarray(output[0])
    res = garson(weights, output)
    sns.set(color_codes=True)
    sns.set_style('white')
    ax = sns.distplot(res[1:1000], bins=30, hist=True, rug=False, label='MF network', hist_kws={'range': [0, 0.003]})
    ax = sns.distplot(res[1001:], bins=30, hist=True, rug=False, label="GE network", hist_kws={'range': [0, 0.003]})
    ax.grid(False)
    ax.set(xticks=np.arange(0,0.003,0.0005))
    ax.set_xlabel('Weight', fontsize=13, fontweight='bold')
    ax.set_ylabel('Absolute frequency', fontsize=13, fontweight='bold')
    ax.set_title('MMDNN weight distribution - GE and MF', fontdict={'fontsize': 15, 'fontweight': 'bold'})
#    ax2.grid(False)
    plts.legend(fontsize=13)
    fig = ax.get_figure()
    fig.savefig("layer_examine_full.pdf", format='pdf', dpi=600)
    plts.clf()

if not os.path.exists("layer_examine_met.pdf"):
    plts.figure()
    weights = multi_model_metabolic_expression.get_layer(name = "last_hidden").get_weights()[0]
    output = multi_model_metabolic_expression.layers[-1].get_weights()
    np.set_printoptions(threshold=20000)
    output = np.asarray(output[0])
    res = garson(weights, output)
    sns.set(color_codes=True)
    sns.set_style('white')
    ax = sns.distplot(res[1:1000], bins=35, hist=True, rug=False, label='MF network', hist_kws={'range': [0, 0.0035]})
    ax = sns.distplot(res[1001:], bins=35, hist=True, rug=False, label="MGE network", hist_kws={'range': [0, 0.0035]})
    ax.grid(False)
    ax.set(xticks=np.arange(0,0.004,0.0005))
    ax.set_xlabel('Weight', fontsize=13, fontweight='bold')
    ax.set_ylabel('Absolute frequency', fontsize=13, fontweight='bold')
    ax.set_title('MMDNN weight distribution - MGE and MF', fontdict={'fontsize': 15, 'fontweight': 'bold'})
#    ax2.grid(False)
    plts.legend(fontsize=13)
    fig = ax.get_figure()
    fig.savefig("layer_examine_met.pdf", format='pdf', dpi=600)
    plts.clf()



plot_model(
    multi_model_full_expression,
    to_file='multi_model.png',
    show_shapes=True,
    show_layer_names=True,
    rankdir='TB'
)
plot_model(
    multi_model_metabolic_expression,
    to_file='multi_model_metabolic_expression.png',
    show_shapes=True,
    show_layer_names=True,
    rankdir='TB'
)
