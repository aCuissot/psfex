from psf import *
import tensorflow as tf

import sys
# to check we use the venv
print(sys.executable)

# TODO: create a dataloader, but we need a dataset previously.

model = tf.keras.models.Sequential([
    tf.keras.layers.Conv2D(32, kernel_size=3, activation='relu', input_shape=(15,15, 1), padding='same'),
    tf.keras.layers.MaxPool2D(),
    tf.keras.layers.Conv2D(64, kernel_size=3, activation='relu', padding='same'),
    tf.keras.layers.MaxPool2D(),
    tf.keras.layers.Conv2D(128, kernel_size=3, activation='relu', padding='same'),
    tf.keras.layers.MaxPool2D(),
 #   tf.keras.layers.Flatten(input_shape=(28, 28)),
    tf.keras.layers.Dense(2048, activation=tf.nn.relu),
   # tf.keras.layers.Dropout(0.2),
    tf.keras.layers.Dense(10, activation=tf.nn.softmax)
])


opt=tf.keras.optimizers.Adam(0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
model.compile(optimizer=opt,
              loss=tf.losses.huber_loss,
              metrics=['accuracy'])
print(model.summary())

#model.fit(x_train, y_train, epochs=5)
#model.evaluate(x_test, y_test) 

