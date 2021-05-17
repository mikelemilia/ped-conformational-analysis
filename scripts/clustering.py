from sklearn_extra import cluster
from sklearn.metrics import silhouette_score
import pandas
import numpy as np
import matplotlib.pyplot as plt
from functions import *

folder = 'data/features_files/PED00020e001/'
features_files = extract_filenames(folder, "csv")
points = []
for features_file in features_files:

    path = folder + features_file + '.csv'
    df = pandas.read_csv(path, index_col=0, header=None, error_bad_lines=False)
    print(df)

    features = []
    for i in range(df.shape[0]):
        features.append(df.iloc[i][1])

    points.append(features)


# # -------------------------------------------------------------------------------- esempio clustering
#
# # random points
# num_points = 10000
# coords_x = np.array(np.random.randint(-50, 50, num_points))
# coords_y = np.array(np.random.randint(-50, 50, num_points))
# coords = np.column_stack((coords_y, coords_x))
#
# clustering and silhouette
# K_set = range(2, 3)
# silhouettes = []
#
# for k in K_set:
#     kMed = cluster.KMedoids(n_clusters=k, metric="euclidean")
#     labels = kMed.fit_predict(coords)
#
#     s = silhouette_score(coords, labels)
#     silhouettes.append(s)
#
#     # for i in u_labels:
#     #     plt.scatter(coords[labels == i, 0], coords[labels == i, 1], label=i)
#     # plt.legend()
#     # plt.show()
#
# print(silhouettes)
# s_max = max(silhouettes)
# k_opt = K_set[silhouettes.index(s_max)]
#
# print(k_opt)
# print(s_max)
