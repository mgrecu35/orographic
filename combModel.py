from keras.models import Sequential
from keras.layers import Merge, Activation, Dense
from keras.layers.recurrent import LSTM

n, m = 10, 3
input_dim = 10
output_dim = 20

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

# Define two input layers
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
