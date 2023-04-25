import numpy as np

tol = 1e-7

def tableau_simplex(T,n,m,basis):
    while True:
        j = np.argmax(T[-1,:-1])
        if T[-1,j] <= tol:
            # optimal solution found
            x = np.zeros(n)
            x[basis < n] = T[basis < n, -1]
            return T,basis, x
        # step 2: select leaving variable
        i = np.argmin(T[:-1,j] / T[:-1,-1]) #Doubt about correctness of this
        if T[i,j] <= 0:
            # unbounded solution
            return None
        # step 3: perform pivot
        basis[i] = j
        T[i,:] /= T[i,j]
        for k in range(m+1):
            if k != i:
                T[k,:] -= T[k,j] * T[i,:]

def gomory(c, A, b):
    m, n_old = A.shape
    n = n_old + m
    # Converting to standard form Ax = b
    A = np.hstack((A, np.eye(A.shape[0])))
    c = np.hstack((c, np.zeros(A.shape[0])))
    # Minimising c
    c = -c

    # making b positive
    for i in range(m):
        if(b[i] < 0):
            b[i] = -b[i]
            A[i,:] = -A[i,:]
    
    # Phase1
    # Ax + Iy = b
    c_phase1 = np.concatenate((np.zeros(n), np.ones(m)))
    A_phase1 = np.hstack((A, np.eye(m)))
    basis_indices = np.arange(n, n+m)
    
    # Constructing initial tableau
    T = np.zeros((m+1, n+m+1))
    # T[1:, 1:] = A
    # T[1:,0] = b
    # T[0, 1:] = c_phase1 - np.ones(m).T @ A
    # T[0,0] = -np.ones(m) @ b
    T[:-1, :-1] = A
    T[0:-1:,-1] = b
    T[-1, :-1] = c_phase1 - np.ones(m).T @ A
    T[-1,-1] = -np.ones(m) @ b

    # Solving the initial tableau
    T, basis_indices, x = tableau_simplex(T, n, m, basis_indices)

    

    



if __name__ == "__main__":

    n, m = map(int, input().split())
    b = np.array(list(map(int, input().split())))
    c = np.array(list(map(int, input().split())))

    A = np.empty((m, n))
    for i in range(m):
        row = np.array(list(map(int, input().split())))
        A[i, :] = row

    print(gomory(c,A,b))
