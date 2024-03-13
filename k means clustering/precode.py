import numpy as np
import random

data = np.load('AllSamples.npy')
id = "9575"

def initial_point_idx(k,N):
	return np.random.RandomState(seed=(id+k)).permutation(N)[:k]

def init_point(data, idx):
    return data[idx,:]

def initial_S1(k):
    # print("Strategy 1: k and initial points")
    i = int(id)%150 
    random.seed(i+500)
    init_idx = initial_point_idx(i,k,data.shape[0])
    init_s1 = init_point(data, init_idx)
    return init_s1

def initial_point_idx2(k, N):
    random.seed(int(id)+k)
    return random.randint(0,N-1)

def initial_S2(idx,k):
    # print("Strategy 2: k and initial points")
    i = int(idx)%150
    random.seed(i+800)
    init_idx2 = initial_point_idx2(i, k,data.shape[0])
    init_s2 = data[init_idx2,:]
    return init_s2
