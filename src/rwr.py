"""General functions to run RWR algorithm
Code adapted from pyrwr: https://github.com/jinhongjung/pyrwr.git"""

from scipy.sparse import spdiags
import numpy as np
from numpy.linalg import norm


def iterate(A, q, c=0.15, epsilon=1e-9,
            max_iters=100, norm_type=1):
    """
    Perform power iteration for RWR, PPR, or PageRank

    inputs
        A : csr_matrix
            input matrix (for RWR and it variants, it should be row-normalized)
        q : ndarray
            query vector
        c : float
            restart probability
        epsilon : float
            error tolerance for power iteration
        max_iters : int
            maximum number of iterations for power iteration
        handles_deadend : bool
            if true, it will handle the deadend issue in power iteration
            otherwise, it won't, i.e., no guarantee for sum of RWR scores
            to be 1 in directed graphs
        norm_type : int
            type of norm used in measuring residual at each iteration
    outputs
        x : ndarray
            result vector
    """
    x = q
    old_x = q
    residuals = np.zeros(max_iters)

    for i in range(max_iters):
        x = (1 - c) * (A.dot(old_x)) + (c * q)

        residuals[i] = norm(x - old_x, norm_type)

        if residuals[i] <= epsilon:
            break

        old_x = x


    return x, residuals[0:i + 1]


class RWR:
    normalized = False
    def __init__(self, A):
        self.A = A
        self.m, self.n = self.A.shape
        self.node_ids = np.arange(0, self.n)
        self.normalize()

    def normalize(self):
        '''
        Perform row-normalization of the adjacency matrix
        '''
        if self.normalized is False:
            nA = self.row_normalize(self.A)
            self.nAT = nA.T
            self.normalized = True

    def row_normalize(self, A):
        '''
        Perform row-normalization of the given matrix

        inputs
            A : csr_matrix
                (n x n) input matrix where n is # of nodes
        outputs
            nA : crs_matrix
                 (n x n) row-normalized matrix
        '''
        n = A.shape[0]

        # do row-wise sum where d is out-degree for each node
        d = A.sum(axis=1)
        d = np.asarray(d).flatten()

        # handle 0 entries in d
        d = np.maximum(d, np.ones(n))
        invd = 1.0 / d

        invD = spdiags(invd, 0, n, n)

        # compute row normalized adjacency matrix by nA = invD * A
        nA = invD.dot(A)

        return nA

    def compute(self, seed_vector, c, epsilon, max_iters):
        '''
        Compute the RWR score vector w.r.t. the seed node

        inputs
            seed : int
                seed (query) node id
            c : float
                restart probability
            epsilon : float
                error tolerance for power iteration
            max_iters : int
                maximum number of iterations for power iteration
            handles_deadend : bool
                if true, it will handle the deadend issue in power iteration
                otherwise, it won't, i.e., no guarantee for sum of RWR scores
                to be 1 in directed graphs
        outputs
            r : ndarray
                RWR score vector
        '''
        self.normalize()
        q = seed_vector
        r, residuals = iterate(self.nAT, q, c, epsilon, max_iters)

        return r

    def check_seeds(self, seeds):
        for seed in seeds:
            if seed < 0 or seed >= self.n:
                raise ValueError('Out of range of seed node id')


