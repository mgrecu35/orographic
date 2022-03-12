from keras.models import Sequential
from keras.layers.recurrent import LSTM
import tensorflow as tf
from tensorflow.keras.layers import *
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.optimizers import Adam, RMSprop


#first_model = Sequential()
#first_model.add(LSTM(output_dim, input_shape=(m, input_dim)))

#second_model = Sequential()
#second_model.add(LSTM(output_dim, input_shape=(n-m, input_dim)))
#
#model = Sequential()
#model.add(Merge([first_model, second_model], mode='concat'))
#model.add(Dense(1))
#model.add(Activation('sigmoid'))
#model.compile(optimizer='RMSprop', loss='binary_crossentropy')
#model.fit([X[:,:m,:], X[:,m:,:]], y)
from keras.layers import Input, Concatenate, Conv2D, Flatten, Dense
from keras.models import Model

def lstm_model(ndims=2,ninp=1):
    ntimes=None
    inp1 = tf.keras.layers.Input(shape=(ntimes,ndims,))
    inp2 = tf.keras.layers.Input(shape=(ninp,))
    out1 = tf.keras.layers.LSTM(12, return_sequences=True)(inp1)
    out1 = tf.keras.layers.LSTM(12, return_sequences=False)(out1)
    out1 = tf.keras.layers.Dense(1)(out1)
    concat_layer= Concatenate()([out1, inp2])
    out1 = tf.keras.layers.Dense(6,activation='relu')(concat_layer)
    out1 = tf.keras.layers.Dropout(0.1) (out1)
    out1 = tf.keras.layers.Dense(6,activation='relu')(out1)
    out1 = tf.keras.layers.Dropout(0.1) (out1)
    out = tf.keras.layers.Dense(1)(out1)
    model = tf.keras.Model(inputs=[inp1,inp2], outputs=out)
    return model

#mod=lstm_model(1,1)

#ndims=1
#ninp=1
#ntimes=2
#inp1 = tf.keras.layers.Input(shape=(ntimes,ndims,))
#inp2 = tf.keras.layers.Input(shape=(ninp,))
#out1 = tf.keras.layers.LSTM(12, return_sequences=True)(inp1)
#out1 = tf.keras.layers.LSTM(12, return_sequences=False)(out1)
#out1 = tf.keras.layers.Dense(1)(out1)
#out1 = Flatten()(out1)
#concat_layer= Concatenate()([out1, inp2])
#out1 = tf.keras.layers.Dense(6)(concat_layer)
# Defin1e two input layers
#image_input = Input((32, 32, 3))
#vector_input = Input((6,))

# Convolution + Flatten for the image
#conv_layer = Conv2D(32, (3,3))(image_input)
#flat_layer = Flatten()(conv_layer)

# Concatenate the convolutional features and the vector input
#concat_layer= Concatenate()([vector_input, flat_layer])
#output = Dense(3)(concat_layer)

# define a model with a list of two inputs
#model = Model(inputs=[image_input, vector_input], outputs=output)
