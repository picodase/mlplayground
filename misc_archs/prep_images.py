
# %% [markdown]
# This tutorial shows how to load and preprocess an image dataset in three ways. First, you will use high-level Keras preprocessing [utilities](https://www.tensorflow.org/api_docs/python/tf/keras/preprocessing/image_dataset_from_directory) and [layers](https://www.tensorflow.org/api_docs/python/tf/keras/layers/experimental/preprocessing) to read a directory of images on disk. Next, you will write your own input pipeline from scratch using [tf.data](https://www.tensorflow.org/guide/data). Finally, you will download a dataset from the large [catalog](https://www.tensorflow.org/datasets/catalog/overview) available in [TensorFlow Datasets](https://www.tensorflow.org/datasets).
# %% [markdown]
# ## Setup

# %%
import numpy as np
import os
import PIL
import PIL.Image
import tensorflow as tf
import tensorflow_datasets as tfds

# %%
print(tf.__version__)

# %%
import pathlib
#dataset_url = "https://storage.googleapis.com/download.tensorflow.org/example_images/flower_photos.tgz"

# Set dataset directory
#data_dir = tf.keras.utils.get_file(origin=dataset_url, 
#                                   fname='flower_photos', 
#                                   untar=True)

data_dir = pathlib.Path(data_dir)       # convert path variable if necessary

# %%
# count the number of images, specifically those that are .png files
image_count = len(list(data_dir.glob('*/*.png')))
print(image_count)

# %% [markdown]
# ## Load using keras.preprocessing
# 
# Let's load these images off disk using [image_dataset_from_directory](https://www.tensorflow.org/api_docs/python/tf/keras/preprocessing/image_dataset_from_directory).
# %% [markdown]
# Note: The Keras Preprocesing utilities and layers introduced in this section are currently experimental and may change.
# %% [markdown]
# ### Create a dataset
# %% [markdown]
# Define some parameters for the loader:

# %%
batch_size = 32
img_height = 128
img_width = 128

# %% [markdown]
# It's good practice to use a validation split when developing your model. We will use 80% of the images for training, and 20% for validation.

# %%
# Load the training dataset
train_ds = tf.keras.preprocessing.image_dataset_from_directory(
  data_dir,
  validation_split=0.2,
  subset="training",
  seed=123,
  image_size=(img_height, img_width),
  batch_size=batch_size)


# %%
# Load the validation dataset
val_ds = tf.keras.preprocessing.image_dataset_from_directory(
  data_dir,
  validation_split=0.2,
  subset="validation",
  seed=123,
  image_size=(img_height, img_width),
  batch_size=batch_size)

# %%
class_names = train_ds.class_names
print(class_names)

# %% [markdown]
# ### Visualize the data
# 
# Here are the first 9 images from the training dataset.

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 10))
for images, labels in train_ds.take(1):
  for i in range(9):
    ax = plt.subplot(3, 3, i + 1)
    plt.imshow(images[i].numpy().astype("uint8"))
    plt.title(class_names[labels[i]])
    plt.axis("off")

# %% [markdown]
# You can train a model using these datasets by passing them to `model.fit` (shown later in this tutorial). If you like, you can also manually iterate over the dataset and retrieve batches of images:

# %%
for image_batch, labels_batch in train_ds:
  print(image_batch.shape)
  print(labels_batch.shape)
  break

# %% [markdown]
# The `image_batch` is a tensor of the shape `(32, 180, 180, 3)`. This is a batch of 32 images of shape `180x180x3` (the last dimension referes to color channels RGB). The `label_batch` is a tensor of the shape `(32,)`, these are corresponding labels to the 32 images. 
# 
# %% [markdown]
# Note: you can call `.numpy()` on either of these tensors to convert them to a `numpy.ndarray`.
# %% [markdown]
# ### Standardize the data
# 
# %% [markdown]
# The RGB channel values are in the `[0, 255]` range. This is not ideal for a neural network; in general you should seek to make your input values small. Here, we will standardize values to be in the `[0, 1]` by using a Rescaling layer.

# %%
from tensorflow.keras import layers

normalization_layer = tf.keras.layers.experimental.preprocessing.Rescaling(1./255)

# %% [markdown]
# There are two ways to use this layer. You can apply it to the dataset by calling map:

# %%
normalized_ds = train_ds.map(lambda x, y: (normalization_layer(x), y))
image_batch, labels_batch = next(iter(normalized_ds))
first_image = image_batch[0]
# Notice the pixels values are now in `[0,1]`.
print(np.min(first_image), np.max(first_image)) 

# %% [markdown]
# Or, you can include the layer inside your model definition to simplify deployment. We will use the second approach here.
# %% [markdown]
# Note: If you would like to scale pixel values to `[-1,1]` you can instead write `Rescaling(1./127.5, offset=-1)`
# %% [markdown]
# Note: we previously resized images using the `image_size` argument of `image_dataset_from_directory`. If you want to include the resizing logic in your model, you can use the [Resizing](https://www.tensorflow.org/api_docs/python/tf/keras/layers/experimental/preprocessing/Resizing) layer instead.
# 
# %% [markdown]
# ### Configure the dataset for performance
# 
# Let's make sure to use buffered prefetching so we can yield data from disk without having I/O become blocking. These are two important methods you should use when loading data.
# 
# `.cache()` keeps the images in memory after they're loaded off disk during the first epoch. This will ensure the dataset does not become a bottleneck while training your model. If your dataset is too large to fit into memory, you can also use this method to create a performant on-disk cache.
# 
# `.prefetch()` overlaps data preprocessing and model execution while training. 
# 
# Interested readers can learn more about both methods, as well as how to cache data to disk in the [data performance guide](https://www.tensorflow.org/guide/data_performance#prefetching).

# %%
AUTOTUNE = tf.data.experimental.AUTOTUNE

train_ds = train_ds.cache().prefetch(buffer_size=AUTOTUNE)
val_ds = val_ds.cache().prefetch(buffer_size=AUTOTUNE)

# %% [markdown]
# ### Train a model
# 
# For completeness, we will show how to train a simple model using the datasets we just prepared. This model has not been tuned in any way - the goal is to show you the mechanics using the datasets you just created. To learn more about image classification, visit this [tutorial](https://www.tensorflow.org/tutorials/images/classification).

### APPEND YOUR DESIGNED MODEL HERE!!!

# %%
num_classes = 5

model = tf.keras.Sequential([
  layers.experimental.preprocessing.Rescaling(1./255),
  layers.Conv2D(32, 3, activation='relu'),
  layers.MaxPooling2D(),
  layers.Conv2D(32, 3, activation='relu'),
  layers.MaxPooling2D(),
  layers.Conv2D(32, 3, activation='relu'),
  layers.MaxPooling2D(),
  layers.Flatten(),
  layers.Dense(128, activation='relu'),
  layers.Dense(num_classes)
])


# %%
model.compile(
  optimizer='adam',
  loss=tf.losses.SparseCategoricalCrossentropy(from_logits=True),
  metrics=['accuracy'])

# %% [markdown]
# Note: we will only train for a few epochs so this tutorial runs quickly. 

# %%
model.fit(
  train_ds,
  validation_data=val_ds,
  epochs=3
)

# %% [markdown]
# Note: you can also write a custom training loop instead of using `model.fit`. To learn more, visit this [tutorial](https://www.tensorflow.org/guide/keras/writing_a_training_loop_from_scratch).
# %% [markdown]
# You may notice the validation accuracy is low to the compared to the training accuracy, indicating our model is overfitting. You can learn more about overfitting and how to reduce it in this [tutorial](https://www.tensorflow.org/tutorials/keras/overfit_and_underfit).
# %% [markdown]
# ## Using tf.data for finer control
# %% [markdown]
# The above keras.preprocessing utilities are a convenient way to create a `tf.data.Dataset` from a directory of images. For finer grain control, you can write your own input pipeline using `tf.data`. This section shows how to do just that, beginning with the file paths from the zip we downloaded earlier.

# %%
list_ds = tf.data.Dataset.list_files(str(data_dir/'*/*'), shuffle=False)
list_ds = list_ds.shuffle(image_count, reshuffle_each_iteration=False)


# %%
for f in list_ds.take(5):
  print(f.numpy())

# %% [markdown]
# The tree structure of the files can be used to compile a `class_names` list.

# %%
class_names = np.array(sorted([item.name for item in data_dir.glob('*') if item.name != "LICENSE.txt"]))
print(class_names)

# %% [markdown]
# Split the dataset into train and validation:

# %%
val_size = int(image_count * 0.2)
train_ds = list_ds.skip(val_size)
val_ds = list_ds.take(val_size)

# %% [markdown]
# You can see the length of each dataset as follows:

# %%
print(tf.data.experimental.cardinality(train_ds).numpy())
print(tf.data.experimental.cardinality(val_ds).numpy())

# %% [markdown]
# Write a short function that converts a file path to an `(img, label)` pair:

# %%
def get_label(file_path):
  # convert the path to a list of path components
  parts = tf.strings.split(file_path, os.path.sep)
  # The second to last is the class-directory
  one_hot = parts[-2] == class_names
  # Integer encode the label
  return tf.argmax(one_hot)


# %%
def decode_img(img):
  # convert the compressed string to a 3D uint8 tensor
  img = tf.image.decode_jpeg(img, channels=3)
  # resize the image to the desired size
  return tf.image.resize(img, [img_height, img_width])


# %%
def process_path(file_path):
  label = get_label(file_path)
  # load the raw data from the file as a string
  img = tf.io.read_file(file_path)
  img = decode_img(img)
  return img, label

# %% [markdown]
# Use `Dataset.map` to create a dataset of `image, label` pairs:

# %%
# Set `num_parallel_calls` so multiple images are loaded/processed in parallel.
train_ds = train_ds.map(process_path, num_parallel_calls=AUTOTUNE)
val_ds = val_ds.map(process_path, num_parallel_calls=AUTOTUNE)


# %%
for image, label in train_ds.take(1):
  print("Image shape: ", image.numpy().shape)
  print("Label: ", label.numpy())

# %% [markdown]
# ### Configure dataset for performance
# %% [markdown]
# To train a model with this dataset you will want the data:
# 
# * To be well shuffled.
# * To be batched.
# * Batches to be available as soon as possible.
# 
# These features can be added using the `tf.data` API. For more details, see the [Input Pipeline Performance](../../guide/performance/datasets) guide.

# %%
def configure_for_performance(ds):
  ds = ds.cache()
  ds = ds.shuffle(buffer_size=1000)
  ds = ds.batch(batch_size)
  ds = ds.prefetch(buffer_size=AUTOTUNE)
  return ds

train_ds = configure_for_performance(train_ds)
val_ds = configure_for_performance(val_ds)

# %% [markdown]
# ### Visualize the data
# 
# You can visualize this dataset similarly to the one you created previously.

# %%
image_batch, label_batch = next(iter(train_ds))

plt.figure(figsize=(10, 10))
for i in range(9):
  ax = plt.subplot(3, 3, i + 1)
  plt.imshow(image_batch[i].numpy().astype("uint8"))
  label = label_batch[i]
  plt.title(class_names[label])
  plt.axis("off")

# %% [markdown]
# ### Continue training the model
# 
# You have now manually built a similar `tf.data.Dataset` to the one created by the `keras.preprocessing` above. You can continue training the model with it. As before, we will train for just a few epochs to keep the running time short.

# %%
model.fit(
  train_ds,
  validation_data=val_ds,
  epochs=3
)

# %% [markdown]
# ## Using TensorFlow Datasets
# 
# So far, this tutorial has focused on loading data off disk. You can also find a dataset to use by exploring the large [catalog](https://www.tensorflow.org/datasets/catalog/overview) of easy-to-download datasets at [TensorFlow Datasets](https://www.tensorflow.org/datasets). As you have previously loaded the Flowers dataset off disk, let's see how to import it with TensorFlow Datasets. 
# %% [markdown]
# Download the flowers [dataset](https://www.tensorflow.org/datasets/catalog/tf_flowers) using TensorFlow Datasets.

# %%
(train_ds, val_ds, test_ds), metadata = tfds.load(
    'tf_flowers',
    split=['train[:80%]', 'train[80%:90%]', 'train[90%:]'],
    with_info=True,
    as_supervised=True,
)

# %% [markdown]
# The flowers dataset has five classes.

# %%
num_classes = metadata.features['label'].num_classes
print(num_classes)

# %% [markdown]
#  Retrieve an image from the dataset.

# %%
get_label_name = metadata.features['label'].int2str

image, label = next(iter(train_ds))
_ = plt.imshow(image)
_ = plt.title(get_label_name(label))

# %% [markdown]
# As before, remember to batch, shuffle, and configure each dataset for performance.

# %%
train_ds = configure_for_performance(train_ds)
val_ds = configure_for_performance(val_ds)
test_ds = configure_for_performance(test_ds)

# %% [markdown]
# You can find a complete example of working with the flowers dataset and TensorFlow Datasets by visiting the [Data augmentation](https://www.tensorflow.org/tutorials/images/data_augmentation) tutorial.
# %% [markdown]
# ## Next steps
# 
# This tutorial showed two ways of loading images off disk. First, you learned how to load and preprocess an image dataset using Keras preprocessing layers and utilities. Next, you learned how to write an input pipeline from scratch using tf.data. Finally, you learned how to download a dataset from TensorFlow Datasets. As a next step, you can learn how to add data augmentation by visiting this [tutorial](https://www.tensorflow.org/tutorials/images/data_augmentation). To learn more about tf.data, you can visit this [guide](https://www.tensorflow.org/guide/data).

