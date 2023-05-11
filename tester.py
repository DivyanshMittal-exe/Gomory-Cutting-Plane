import unittest
import numpy as np
from ilp import gomory_cab
from scipy.optimize import linprog


class TestLinearProgramming(unittest.TestCase):

    def test_linear_programming(self):
        
        while(True):
            n = np.random.randint(1, 10)
            m = np.random.randint(1, 10)

            A = np.random.randint(-1000, 1001, size=(m, n))
            
            # I do not know what this line does
            b = np.random.randint(-1000, 1001, size=m)
            # b = np.random.randint(np.amax(A, axis=1), 1001, size=m)

            
            # b = np.random.randint(-1000, 1001, size=m)
            c = np.random.randint(-1000, 1001, size=n)


            int_constraint = [1]*n
            lin_prog_solve = linprog(-c, A_ub=A, b_ub=b,method="highs",integrality=int_constraint)
            if lin_prog_solve.x is None:
                continue
            
            print(lin_prog_solve)
            # return
            our_out = gomory_cab(c, A, b)
            
            print(lin_prog_solve)
            

            int_constraint = [1] * n
            
            print("Correct Output - ", lin_prog_solve.x, c@lin_prog_solve.x)
            print("Our Output - ", our_out, c@our_out)
            # if lin_prog_solve.x != our_out


            if lin_prog_solve.x is not None:
                break
        
        print(our_out)
        
        # lin_prog_solve = [int(i) for i in lin_prog_solve.x]
        # print(lin_prog_solve)
        np.testing.assert_allclose(lin_prog_solve.x, our_out, rtol=1e-9)

        # self.assertTrue(False)
        # self.assertTrue(np.array_equal(lin_prog_solve.x, our_out))


if __name__ == "__main__":
    unittest.main()