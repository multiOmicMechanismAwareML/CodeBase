import csv
from keras.callbacks import EarlyStopping
import tensorflow as tf
import numpy as np
import pandas as pa
import os.path
from keras.utils import np_utils
import keras.callbacks
import matplotlib.pyplot as plts
import seaborn as sns
from keras.models import Sequential, Model
from keras.layers.core import Dense, Dropout, Activation
import matplotlib.pyplot
import keras.backend as k
from keras.layers.advanced_activations import LeakyReLU
from keras.layers import Concatenate, Input
from keras.optimizers import SGD
from keras.constraints import max_norm
from sklearn.model_selection import train_test_split
from keras.callbacks import Callback
from sklearn import preprocessing
import matplotlib.pyplot
from keras.utils import plot_model

SEED = 250

#Checks that the reading of the data has worked
def check_file_read_ok(data, name):
    if data is None:
        print("error with " + name + " data not read")
    else:
        print("Success in loading " + name)
        print(data.shape)

def remove_zero_entry_columns(data):
    return data.loc[:, (data != 0).any(axis = 0)]

def remove_tiny_columns(data, threshold):
    cols_to_keep = abs(data.sum(axis = 0)) > threshold
    print(data.columns[cols_to_keep])
    return data[data.columns[cols_to_keep]]
#Name of the row we are trying to predict (growth rate)
TARGET_NAME = 'log2relT'

with open('data/testing_index.csv', 'r') as csvfile:
    testing_index = []
    for row in csv.reader(csvfile, delimiter=';'):
        testing_index.append(row[0]) # careful here with [0]

#Load the data removing any
full_data = pa.read_csv('data/completeDataset.csv')
check_file_read_ok(full_data, "full data")
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
#metabolic_expression_data = metabolic_expression_data.drop(columns=TARGET_NAME)
genes = full_data['Row']
full_data = full_data.drop(columns = 'Row')
reaction_data = full_data.drop(columns = expression_data.columns.values)
#reaction_data = remove_tiny_columns(reaction_data, 30)

def add_reactions_to_list(lis, reactions):
    for x in reactions:
        lis.append(x)
    return lis

reactions_to_examine = pa.read_csv('fluxes_to_examine.csv')
#start with growth and any reactions that mention ATP, then add the addtional pathways
#general = ['r_4041']
glycolyses = pa.read_csv("data/glycolyses.csv", skiprows = 1).values
citrate = pa.read_csv("data/citrate.csv", skiprows = 1).values
fatty_acid = pa.read_csv("data/fatty_acid.csv", skiprows = 1).values
amino_acid = pa.read_csv("data/amino_acid.csv", skiprows = 1).values
oxi_phos = pa.read_csv("data/oxi_phos.csv", skiprows = 1).values
pentose = pa.read_csv("data/pentose.csv", skiprows = 1).values
pyruvate = pa.read_csv("data/pyruvate.csv", skiprows = 1).values
pathways = []


pathways = {"glycolyses" : glycolyses,
 "citrate" : citrate,
  "fatty_acid" : fatty_acid,
  "amino_acid" : amino_acid,
   "oxi_phos" : oxi_phos,
    "pentose" : pentose,
     "pyruvate" : pyruvate}

def flatten(x):
    return [t for s in x for t in s]

pathway_names = [ "glycolyses", "citrate", "fatty_acid", "amino_acid", "oxi_phos", "pentose", "pyruvate"]
#pathway_names = ["citrate", "fatty_acid", "amino_acid", "oxi_phos", "pentose", "pyruvate"]

joined_reactions = (glycolyses, citrate, fatty_acid, amino_acid, oxi_phos, pentose, pyruvate)

#joined_reactions = (citrate, fatty_acid, amino_acid, oxi_phos, pentose, pyruvate)
joined_reactions = flatten(np.concatenate(joined_reactions))
joined_i = np.unique(joined_reactions, return_index = True) [1]
joined_reactions = [joined_reactions[index] for index in sorted(joined_i)]

print(joined_reactions)

reaction_data = reaction_data[joined_reactions]
print(reaction_data.columns, "col names ")



def init_model(input_dim, learning_rate, epochs, momentum,initial_neurons, pathway_neurons,  neurons, pathway_dic, pathway_names, trainable = True):
    outputs = []
    added = set()
    input = Input(shape = (input_dim,))
    initial = Dense(initial_neurons, activation='sigmoid', kernel_constraint=max_norm(4), name = "intial") (input)
    initial = Dropout(rate=0.6) (initial)
    #initial = Dense(initial_neurons, activation='sigmoid', kernel_constraint=max_norm(4), name = "intial2") (initial)
    for x in range(len(pathway_names)):
        layer_f = Dense(pathway_neurons, activation='sigmoid', kernel_constraint=max_norm(3), name = pathway_names[x] + "_1") (initial)
        #layerf = Dropout(rate=0.6) (layer_f)
        # layer = Dense(pathway_neurons, activation='sigmoid', kernel_constraint=max_norm(3), name = pathway_names[x] + "_2") (layerf)
        fluxes_in_pathway = flatten(pathway_dic[pathway_names[x]])
        layer_f = Dropout(rate=0.6) (layer_f)
        for i in range(len(fluxes_in_pathway)):
            if(fluxes_in_pathway[i] not in added):
                added.add(fluxes_in_pathway[i])
                layer = Dense(neurons, activation='sigmoid', kernel_constraint=max_norm(4), name = fluxes_in_pathway[i] + "_dense") (layer_f)
            #    test = Dense(neurons, activation='relu') (input)
            #    test = Dropout(rate = 0.6) (input)
            #    test = Dense(neurons, activation='relu') (test)
            #    final = Concatenate()([layer, test])
                outputs.append(Dense(1, activation = 'linear', name = fluxes_in_pathway[i]) (layer))

    model = Model(inputs = input, outputs = outputs)
    rms = SGD(lr= learning_rate, decay= learning_rate / epochs, momentum=momentum)
    model.trainable = trainable
    if (trainable) :
        model.compile(loss='mean_squared_error', optimizer=rms)
    return model



number_of_instances = len(target_data)
testing_index = list(map(int, testing_index))
training_index = np.setxor1d(range(1,number_of_instances), testing_index)
epochs = 10000
batches = 100
lrate = 0.0005
validation = 0.1
print(training_index)
#Split the data 80:20
full_data_train, full_data_test = full_data.drop(full_data.index[testing_index]), full_data.iloc[testing_index, :]
expression_data_train, expression_data_test = expression_data.drop(expression_data.index[testing_index]), expression_data.iloc[testing_index, :]
metabolic_expression_data_train, metabolic_expression_data_test = expression_data_train, expression_data_test
# metabolic_expression_data_train, metabolic_expression_data_test = metabolic_expression_data.drop(metabolic_expression_data.index[testing_index]), metabolic_expression_data.iloc[testing_index, :]
#metabolic_expression_data_train, metabolic_expression_data_test = metabolic_expression_data[training_index,:], metabolic_expression_data[testing_index,:]


print("reaction data size ", reaction_data.shape)
#NUMBER_OF_TARGETS = 500
#reaction_data = reaction_data.sample(NUMBER_OF_TARGETS, axis = 1, random_state = SEED)
#print("new shape ", reaction_data.shape)
#print(reaction_data)
reaction_data_train, reaction_data_test = reaction_data.drop(reaction_data.index[testing_index]), reaction_data.iloc[testing_index, :]


expression_scaler = preprocessing.StandardScaler().fit(expression_data_train)
expression_data_scaled_train = expression_scaler.transform(expression_data_train).astype(np.float32)
expression_data_scaled_test = expression_scaler.transform(expression_data_test).astype(np.float32)

metabolic_expression_scaler = preprocessing.StandardScaler().fit(metabolic_expression_data_train)
metabolic_expression_data_scaled_train = metabolic_expression_scaler.transform(metabolic_expression_data_train).astype(np.float32)
metabolic_expression_data_scaled_test = metabolic_expression_scaler.transform(metabolic_expression_data_test).astype(np.float32)

result = {}
earlyStopping=EarlyStopping(monitor='val_loss', patience=15000, verbose=0, mode='auto')
print("STARTING")

#Make the target dataframe into a format able to be used in tensorflow
def create_the_target_dataset(target_train_data, target_test_data):
    targets_train = []
    targets_test = []
    for x in target_train_data.columns.values:
        targets_train.append(target_train_data[x])
        targets_test.append(target_test_data[x])
    return(targets_train, targets_test)

targets, targets_test = create_the_target_dataset(reaction_data_train, reaction_data_test)

full_fba = init_model(metabolic_expression_data_scaled_train.shape[1],
 lrate,
  5000,
   0.75,
   8000,
   4000,
   500,
    pathways,
    pathway_names)

plot_model(
    full_fba,
    to_file='full_fba.png',
    show_shapes=True,
    show_layer_names=True,
    rankdir='TB'
)

if os.path.exists('Models/full_fba.h5'):
    full_fba.load_weights('Models/full_fba.h5')

else:
    print("STARTING TRAINING")
    full_fba.fit(x = metabolic_expression_data_scaled_train,
     y = targets,
     epochs=epochs,
     batch_size=batches,
     validation_split=validation,
     callbacks=[earlyStopping])
    full_fba.save_weights('Models/full_fba.h5')
score = full_fba.evaluate(metabolic_expression_data_scaled_test, targets_test, verbose=1)
plts.boxplot(score)
plts.show()
print("Set of scores", score)
print("Mean score ", np.mean(score))

predictions = full_fba.predict(metabolic_expression_data_scaled_test)
