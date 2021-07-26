import numpy as np
from scipy.linalg import toeplitz
from pyeda.inter import *
import constraint
import matplotlib.pyplot as plt

class Solver():
    def __init__(self, H, W, N) -> None:
        partialH = toeplitz([1, 1]+[0]*(H-2))
        partialW = toeplitz([1, 1]+[0]*(W-2))
        self.A = np.kron(partialH, partialW)
        self.H = H
        self.W = W
        self.N = N

    def getMove(self, uncovered, partialInfos):
        ravel_unconvered = np.ravel(uncovered)
        #D = np.diag(ravel_unconvered)
        #Dnot = np.diag(1-ravel_unconvered)
        #P = D @ self.A @ Dnot
        P = np.copy(self.A)
        P[:, np.nonzero(ravel_unconvered)[0]] = 0
        P[np.nonzero(1-ravel_unconvered)[0], :] = 0


        plt.matshow(np.linalg.pinv(P))
        plt.colorbar()
        print(np.linalg.pinv(P) @ np.ravel(partialInfos))
        plt.show()
        exit()
        ravel_partialInfos = np.ravel(partialInfos)

        problem = constraint.Problem()

        for rowidx, row in enumerate(P):
            nzrow = np.nonzero(row)[0]
            if len(nzrow) == 0:
                continue
            Xs = ['x'+str(i) for i in nzrow]
            for x in Xs:
                if x not in problem._variables:
                    problem.addVariable(x, [0, 1])
            problem.addConstraint(
                constraint.ExactSumConstraint(ravel_partialInfos[rowidx]), Xs)

        farCells = [i for i in range(self.W*self.H)
                    if (ravel_unconvered[i] == 0) and ('x'+str(i) not in problem._variables)]

        if len(farCells) == 0:
            problem.addConstraint(constraint.ExactSumConstraint(self.N))
        else:
            problem.addConstraint(constraint.MaxSumConstraint(self.N))


        numberOfSolutions = 0
        tallySolutions = {v: 0 for v in problem._variables.keys()}
        solutions = []
        for sol in problem.getSolutionIter():
            solutions.append(sol)
            numberOfSolutions += 1
            for v in sol.keys():
                tallySolutions[v] += sol[v]
            if numberOfSolutions > 1000:
                #TODO: potrebbe convenire colpire sul perimetro
                v = np.random.choice(farCells)
                hopeMove = np.unravel_index(v, (self.H, self.W))
                return False, [hopeMove]
                #exit(-1)

        reverseTally = {}
        for var, times in tallySolutions.items():
            if times not in reverseTally:
                reverseTally[times] = [var]
            else:
                reverseTally[times] += [var]

        if 0 in reverseTally:
            safeMoves = np.unravel_index(
                [int(v[1:]) for v in reverseTally[0]], (self.H, self.W))
            return True, list(zip(*safeMoves))
        else:
            # Conjecture: tutte le celle che non compaiono nel sistema lineare
            # hanno la stessa probabilità di essere una bomba
            # tale probabilità è (bombe rimanenti)/(celle chiuse non sul perimetro)
            # dove bombe rimanenti sono il numero di bombe non adiacenti al perimetro

            bestChance = min(reverseTally.keys())  # best chance sul perimetro

            # numero di bombe non sul perimetro nei vari scenari
            numFarBombs = self.N-np.mean([sum(s.values()) for s in solutions])

            if bestChance*len(farCells) <= numFarBombs*numberOfSolutions:
                v = np.random.choice([int(v[1:]) for v in reverseTally[bestChance]])
            else:
                v = np.random.choice(farCells)

            hopeMove = np.unravel_index(v, (self.H, self.W))
            return False, [hopeMove]

    def print(self, uncovered, partialInfos):
        s = ''
        for r in range(self.H):
            for c in range(self.W):
                if uncovered[r, c] == 0:
                    s += chr(9608)
                else:
                    s += str(partialInfos[r, c])
                s += ' '
            s += '\n'
        return s