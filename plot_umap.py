import sklearn.datasets
import pandas as pd
import numpy as np
import umap
import umap.plot
pendigits = sklearn.datasets.load_digits()
mnist = sklearn.datasets.fetch_openml('mnist_784')
fmnist = sklearn.datasets.fetch_openml('Fashion-MNIST')
mapper = umap.UMAP().fit(pendigits.data)
umap.plot.points(mapper)
