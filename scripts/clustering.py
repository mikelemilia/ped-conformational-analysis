import matplotlib.pyplot as plt
from sklearn_extra import cluster
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA, KernelPCA
from sklearn.cluster import AffinityPropagation
import networkx as nx
import pandas
from functions import *

folder = 'data/features_files/'
features_files = extract_filenames(folder, 'csv')

for features_file in features_files:

    path = folder + features_file + '.csv'
    print(path)
    df = pandas.read_csv(path, index_col=None, header=None)

    p = np.full((df.shape[0], df.shape[1]), None)
    for row in range(0, df.shape[0]):
        p[row, :] = np.array(df.iloc[row])

    points = p[:, 2:]

    pca = PCA(n_components=20)
    pca_points = pca.fit_transform(points)
    # print(np.sum(pca.explained_variance_ratio_))

    # kpca = KernelPCA(n_components=20, kernel='poly')
    # kpca_points = kpca.fit_transform(points)

    K_set = range(3, 9)
    silhouettes = []

    for k in K_set:
        kMed = cluster.KMedoids(n_clusters=k, metric='euclidean', init='k-medoids++', method='pam', max_iter=500)
        labels = kMed.fit_predict(pca_points)

        s = silhouette_score(pca_points, labels)
        silhouettes.append(s)

        # u_labels = np.unique(labels)
        # for i in u_labels:
        #     plt.scatter(points[labels == i, 0], points[labels == i, 1], label=i)
        # plt.legend()
        # plt.show()

    # clustering = AffinityPropagation(random_state=2)
    # clusters = clustering.fit_predict(pca_points)

    s_max = max(silhouettes)
    k_opt = K_set[silhouettes.index(s_max)]

    print(k_opt)
    print(s_max)

    kMed = cluster.KMedoids(n_clusters=k_opt, metric='euclidean', init='k-medoids++')
    labels = kMed.fit_predict(pca_points)
    medoids = kMed.medoid_indices_

    g = nx.Graph()
    for medoid in medoids:
        g.add_node(medoid)
    for i in range(len(medoids)):
        for j in range(i+1, len(medoids)):
            g.add_edge(medoids[i], medoids[j], weight=np.sum(points[medoids[i]]-points[medoids[j]]))

    nx.draw(g, with_labels=True)
    plt.show()

    f = extract_features_csv(p[1,:])
    print(f)
