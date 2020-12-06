import numpy as np
def matmul(*matrices):
    # for multiplying matrices together, note matrices are specified in order
    result = None
    for M in matrices:
        if len(M.shape) == 1: M = np.array([[M[i]] for i in range(len(M))])
        n, p = M.shape
        if type(result) != np.ndarray: result = M
        else:
            m, n = result.shape
            temp = np.array([[0] * p for i in range(m)])
            for i in range(m):
                for j in range(p): temp[i,j] = sum([result[i,k] * M[k,j] for k in range(n)])
            result = temp
    return result
def rotX(theta, M):
    # matrices from https://en.wikipedia.org/wiki/Rotation_matrix
    rotMat = np.array([[1, 0, 0],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]])
    return matmul(rotMat, M)
def rotY(theta, M):
    # matrices from https://en.wikipedia.org/wiki/Rotation_matrix
    rotMat = np.array([[np.cos(theta), 0, np.sin(theta)],
                        [0, 1, 0],
                        [-np.sin(theta), 0, np.cos(theta)]])
    return matmul(rotMat, M)
def rotZ(theta, M):
    # matrices from https://en.wikipedia.org/wiki/Rotation_matrix
    rotMat = np.array([[np.cos(theta), -np.sin(theta), 0],
                       [np.sin(theta), np.cos(theta), 0],
                       [0, 0, 1]])
    return matmul(rotMat, M)