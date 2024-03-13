from precode import *
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
import matplotlib.colors as mcolors
import math

data = np.load('AllSamples.npy')

initial_centers = {}
for k in range(2,11):
    centers = initial_S1(k)
    initial_centers[k] = centers
print(initial_centers)


def get_centroids(ics):
    '''
    This function takes in a dictionary of randomly chosen initial cluster centers (one for each k).
    It generates the remaining initial centroids by sequentially assigning the data point with the
    greatest average distance to the existing centroids as the next centroid assigned until all K
    centroids have been selected.
    :param ics: a dictionary
    :return: inital_cs:
    '''
    initial_cs = ics.copy()
    #for k in initial_cs.keys():
    for k in range(2,11):
        initial_cs[k] = [initial_cs[k]]
        while len(initial_cs[k]) < k:
            point_dists = {}
            idx = 0
            for point in data:  # data[idx,:]
                point_dists[idx] = []
                skip = False
                for centroid in initial_cs[k]:
                    d = math.dist(point, centroid)
                    if d == 0:
                        skip = True
                    point_dists[idx].append(d)
                    # print(centroid, point, d)
                if skip:
                    avg_d = 0
                else:
                    avg_d = sum(point_dists[idx]) / len(point_dists[idx])
                #print(point_dists[idx], avg_d)
                point_dists[idx] = avg_d
                idx += 1
            maxd_idx = (max(point_dists, key=point_dists.get))
            initial_cs[k] = np.vstack((initial_cs[k], data[maxd_idx, :]))

    return initial_cs

def assign_clusters(centers, data, K):
    '''
    This function calculates the distance from each point to the centroid provided. Points are assigned
    to a cluster based on the nearest centroid.
    :param centers: a dictionary with the number of clusters k (key) and the coordinates of the centroids of each cluster (values).
    Values are an array of lists containing x, y point coordinates.
    data: array of lists containing x, y point coordinates.
    K: The number of clusters (int)
    :return:
    '''
    #group points into clusters
    centroid_list = centers[K]
    clusters = {}
    for centroid in centroid_list:
        clusters[str(centroid)] = []
    for point in data:
        x1 = point[0]
        y1 = point[1]
        distances = []
        for centroid in centroid_list:
            x2 = centroid[0]
            y2 = centroid[1]
            d = ((x2-x1)**2 + (y2-y1)**2)**.5 #math.dist(point,centroid)
            distances.append({ "centroid": centroid, "distance": d})
        closest = min(distances, key=lambda x: x['distance'])
        key = str(closest["centroid"])
        clusters[key].append({'x': x1, 'y': y1, 'distance': closest["distance"]})

    return clusters

def get_new_centers(clusters):
    '''
    This function computes the centroid of a given cluster.
    :param clusters: a dictionary containing the existing cluster ids (keys) and the points that belong to each cluster (values)
    :return: a dictionary with the number of clusters k (key) and an array containing the centroids of each cluster (values)
      '''
    new_centroids = {}
    k = len(clusters.values())
    new_centroids[k] = []
    for cluster in clusters.values():
        xs = []
        ys = []
        for member in cluster:
            xs.append(member['x'])
            ys.append(member['y'])
        new_x = sum(xs) / len(xs)
        new_y = sum(ys) / len(ys)
        coords = [new_x, new_y]
        new_centroids[k].append(coords)
    new_centroids[k] = np.array(new_centroids[k])

    return new_centroids

def objective_function(clusters):
    '''
    Calculates the sum of the squares of the distances of each data point to its closest vector.
    :param clusters: a dictionary that has the points associated with each cluster, as well as the
    distance from each point to the centroid of that cluster
    :return: a list of the error for each k (k, error)
    '''
    k = len(clusters.values()) # get the number of clusters (for plotting later)
    js = []
    j = 0
    for cluster in clusters.values():
        dist_sq = []
        for member in cluster:
            dist_sq.append(member["distance"] ** 2)
        j = j + sum(dist_sq)
    #print(k, j)
    js.append([k, j])

    return js

def plot_loss(x, y, subtitle):
    """
    plots the output of the objective function for each k
    :param x: a list of the k values
    :param y: a list of the loss associated with each k
    """

    plt.plot(x, y, color = 'r')
    plt.scatter(x, y, color = 'r')
    for (xi, yi) in zip(x, y):
        plt.text(xi+.20, yi+10, int(yi), va='bottom', ha='center', color = 'black')
    plt.title('Loss vs K')
    plt.xlabel("Number of Clusters (K)")
    plt.ylabel("Calculated Loss")
    plt.suptitle(subtitle)
    plt.show()

def plot_clusters(all_rounds, title):
    '''
    Plots the points and centroid for every k value. Points belonging to a cluster have the same color.
    :param all_rounds: the cluster data for every k iteration
    title: The plot title
    '''

    fig, axs = plt.subplots(3, 3)
    i = 0
    j = 0
    for kval, clusters in all_rounds.items():
        cluster_colors = {}
        cycol = cycle(mcolors.TABLEAU_COLORS)
        for key in clusters.keys():
            cluster_colors[key] = next(cycol)
        for key, value in clusters.items():
            cluster_color = cluster_colors[key]
            print(key, cluster_color)
            xs = []
            ys = []
            for member in value:
                xs.append(member['x'])
                ys.append(member['y'])
            axs[i, j].scatter(xs, ys, color=cluster_color)
            key = key.replace("[", "")
            key = key.replace("]", "")
            axs[i, j].scatter(float(key.split()[0]), float(key.split()[1]), color='black', marker = 'X', sizes = [50])
        axs[i, j].set_title("K = " + str(kval))
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])

        j += 1
        if j >= 3 and i == 0:
            i = 1
            j = 0
        elif j >= 3 and i == 1:
            i = 2
            j = 0

    fig.suptitle(title, fontsize=14)
    plt.show()

def rand_start():
    '''
    Cluster using randomly chosen initial starting points. 
    :return: dictionaries containing the relevant outputs. 'all_centroids' has the coordinates of all of
    the centroids. 'all_error' has the combined error for all clusters for each k. 'all_clusters' has
    all of the centroids and points for every cluster for each k. This can be used for visualizing the
    cluster associations.
    '''
    all_centroids = {}
    all_error = {}
    all_clusters = {}

    for k in range(2, 11):
        initial_results = assign_clusters(initial_centers, data, k)
        results = initial_results
        count = 0

        while count < 5:
            count += 1
            new_centroids = get_new_centers(results)
            new_clusters = assign_clusters(new_centroids, data, k)
            final_centroids = get_new_centers(new_clusters)
            results = assign_clusters(final_centroids, data, k)
            print(k, count)
            print(final_centroids[k][0])

        all_centroids[k] = (final_centroids[k].tolist())
        all_clusters[k] = results
        all_error[k] = (objective_function(results))

    return all_centroids, all_error, all_clusters

def max_d():
    '''
    Perform clustering using the maximum distance approach to selecting the initial centroids.
    :return: dictionaries containing the relevant outputs. 'all_centroids' has the coordinates of all of
    the centroids. 'all_error' has the combined error for all clusters for each k. 'all_clusters' has
    all of the centroids and points for every cluster for each k. This can be used for visualizing the
    cluster associations.
    '''
    all_centroids = {}
    all_error = {}
    all_clusters = {}
    new_initial_centers = get_centroids(initial_centers)

    for k in range(2, 11):
        initial_results = assign_clusters(new_initial_centers, data, k)
        results = initial_results
        count = 0

        while count < 20:
            count += 1
            new_centroids = get_new_centers(results)
            new_clusters = assign_clusters(new_centroids, data, k)
            final_centroids = get_new_centers(new_clusters)
            results = assign_clusters(final_centroids, data, k)
            print(count)
            print(final_centroids)

        all_centroids[k] = (final_centroids[k].tolist())
        all_clusters[k] = results
        all_error[k] = (objective_function(results))

    return all_centroids, all_error, all_clusters


#Plot the output from the max distance clustering 
cents, errs, clusts = max_d()

knum = []
kloss = []
for val in errs.values():
    knum.append(val[0][0])
    kloss.append(val[0][1])

title = 'Max Distance Centers'

#plot loss
plot_loss(knum, kloss, title)

#plot_clusters
plot_clusters(clusts, title)

