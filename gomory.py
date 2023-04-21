import numpy as np


def gomory(c,A,b):
    pass


if __name__ == "__main__":

    n, m = map(int, input().split())
    b = np.array(list(map(int, input().split())))
    c = np.array(list(map(int, input().split())))

    A = np.empty((m, n))
    for i in range(m):
        row = np.array(list(map(int, input().split())))
        A[i, :] = row

    print(gomory(c,A,b))
