from Calcula_GF_Class import Green_Func
from Novo.Classe_Data_Handler import Data_Handler
import matplotlib.pyplot as plt
from tensorflow import keras
import pandas as pd
import numpy as np

N = 11


t_ = np.ones((N, N))
e0 = np.zeros(N)

obj_GF = Green_Func(100, 0.01, N, t_, e0)

_, GF = obj_GF.Calcula_GF()

DF = obj_GF.Armazena(GF ,r'/home/luan/Documents/', 'Teste1.csv')

#data = pd.read_csv(r'/home/luan/Documents/Teste1.csv', sep = ' ')

#DF = pd.DataFrame(data)

dados_form = np.array(DF)

obj_data = Data_Handler(dados_form)

tempo, GF1 = obj_data.handler()

lb = 200

(X_tr, y_tr), (X_te, y_te) = obj_data.prepross(lb, GF1[0][0], val_split = 0.5)

net = keras.Sequential()

net.add(keras.layers.Conv2D(40, 3, padding = 'same', input_shape = (lb, 1, 1)))
net.add(keras.layers.LeakyReLU(alpha = 0.01))
net.add(keras.layers.MaxPooling2D(pool_size = (5, 1)))

net.add(keras.layers.Conv2D(40, 3, padding = 'same'))
net.add(keras.layers.LeakyReLU(alpha = 0.01))
net.add(keras.layers.MaxPooling2D(pool_size = (5, 1)))

net.add(keras.layers.Flatten())

net.add(keras.layers.Dense(40))
net.add(keras.layers.LeakyReLU(alpha = 0.01))

net.add(keras.layers.Dense(1))

net.compile(optimizer = 'adam', loss = 'mse', metrics = ['mae'])

hist = net.fit(X_tr, y_tr, epochs = 300, validation_data = (X_te, y_te))

extrapol = []

seed_batch = y_tr[:lb].reshape((1,lb, 1))
current_batch = seed_batch
for i in range(lb+len(y_te)):
	predicted_value = net.predict(current_batch, verbose = 0)[0]
	extrapol.append(predicted_value) 
	current_batch = np.append(current_batch[:,1:,:],[[predicted_value]],axis=1)

  	
plt.plot(tempo, GF1[0][0])
plt.plot(tempo[int(0.5*(len(tempo) - lb)):]-0.01*lb, extrapol)

plt.show(block = True)

print("Rodei")
