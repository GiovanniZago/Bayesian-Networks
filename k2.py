import itertools

import numpy as np
from scipy.special import gammaln


def construct_unique(df):
    uniques = []
    for i in range(df.shape[1]):
        uniques.append(np.unique(df[:, i]))
    return np.array(uniques)


def g(i, parents, mat, uniques):
    r = len(uniques[i])
    if not parents:
        N = np.unique(mat[:, i], return_counts=True)[1].reshape(1, -1)
    else:
        comb = np.array(list(itertools.product(*uniques[parents], uniques[i])))
        mat_par = mat[:, parents + [i]].reshape(len(mat), 1, -1)
        foo = np.all(np.equal(mat_par, comb.reshape(1, *comb.shape)), axis=2)
        N = np.sum(foo, axis=0).reshape(-1, 2)

    logprod = np.sum(gammaln(N + 1), axis=1)
    the_sum = np.sum(N, axis=1)
    log_out = np.sum((gammaln(r - 1 + 1) - gammaln(the_sum + r - 1 + 1)) + logprod)

    return log_out


def k2(nodes, max_parents, D):
    graph = np.zeros(shape=(len(nodes), len(nodes)))
    parents = [[] for i in range(len(nodes))]
    uniques = construct_unique(D)
    for i in nodes:
        if i == 0:
            continue

        pold = g(i, parents[i], D, uniques)

        while len(parents[i]) < max_parents:
            indexes = set(nodes[:i]) - set(parents[i])

            if not indexes:
                break

            values = [None] * len(indexes)
            for j, ind in enumerate(indexes):
                cand_parents = parents[i] + [ind]
                values[j] = (g(i, list(cand_parents), D, uniques), j)
            values.sort(key=lambda x: x[0], reverse=True)

            if values[0][0] > pold:
                pold = values[0][0]
                parents[i].append(values[0][1])
                graph[values[0][1], nodes[i]] = 1

            else:
                break
    return parents

def main():
    D = np.array([
        [1, 1, 0, 1, 0, 0, 1, 0, 1, 0],
        [0, 1, 0, 1, 0, 1, 1, 0, 1, 0],
        [0, 1, 1, 1, 0, 1, 1, 0, 1, 0]
    ]).T + 1

    
    parents = k2(nodes=np.arange(0, D.shape[1]), max_parents=2, D=D)

    print(parents)

if __name__ == "__main__":
    main()
