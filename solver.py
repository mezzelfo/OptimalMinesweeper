import numpy as np
from scipy.linalg import toeplitz
from pyeda.inter import *
import constraint
import matplotlib.pyplot as plt
import traceback

import networkx as nx
from qpsolvers import solve_qp

from scipy.special import comb
import sympy

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
        ravel_partialInfos = np.ravel(partialInfos)
        #D = np.diag(ravel_unconvered)
        #Dnot = np.diag(1-ravel_unconvered)
        #P = D @ self.A @ Dnot
        P = np.copy(self.A)
        P[np.nonzero(1-ravel_unconvered)[0], :] = 0 # Setto righe a zero tramite moltiplicazione D*Pnot perchè
                                                    # le righe corrispondenti sono accoppiate nel vettore delle informazioni
                                                    # con le celle non ancora scoperte, così da non vincolare il sistema lineare
                                                    # a delle infomazioni non ancora disponibili
                                                    # MATRICE DI ADIACENZA: cancello gli archi uscenti dalle celle non ancora scoperte
        P[:, np.nonzero(ravel_unconvered)[0]] = 0   # Setto colonne a zero tramite moltiplicazione P*D perchè
                                                    # se che le incognite corrispondenti non danno contributo al sistema
                                                    # infatti so già che non sono bombe perchè sono già scoperte
                                                    # MATRICE DI ADIACENZA: cancello gli archi entranti in celle già scoperte
                                                    # => tutti gli archi partono da celle scoperte e vanno in celle non scoperte
                                                    # ho ottenuto la matrice di adiacenza di un grafo bibartito diretto, probabilmente non connesso
                                                    # la presenza di un arco (necessariamente da cella scoperta a cella non scoperta) indica
                                                    # la dipendenza alla cumulata della bomba      
        #print(comb(P.sum(axis = 1), ravel_partialInfos).prod())
        P2 = P[np.ix_(np.sum(P,1) > 0, np.sum(P,0) > 0)]
        m,n = P2.shape
        b = ravel_partialInfos[np.sum(P,1) > 0]
        
        xp,_,rank,__ = np.linalg.lstsq(P2,b,rcond=None)

        if m >= n and rank == n:
            xp = xp.round().astype(np.int32)
            try:
                assert np.allclose(P2 @ xp, b)
                assert np.all((xp == 0) | (xp == 1))
            except:
                print(xp)
                print(P2 @ xp)
                print(b)
                raise RuntimeError()
            print(xp == 0)
            return
        else:
            symP2 = sympy.Matrix(P2)
            nullspacebasis = symP2.nullspace()
            nullspacebasis = np.asarray(nullspacebasis).astype(np.int32)
            nullspacebasis = nullspacebasis.squeeze().T
            try:
                assert np.allclose(P2 @ nullspacebasis,0)
            except:
                print(P2)
                print(nullspacebasis)
                raise RuntimeError()
            print(P2.shape, np.all(nullspacebasis == 0,1).shape)
            return

        # s = np.sum(P2,1)
        # 

        # print('P2 =',P2,';')
        # print(f'b = {b}\';')
        # print(f's = {s}\';')
        # print(P2.shape)
        # print('NEXT')

        # info_rows = s == b
        # keep_rows = s != b
        # variables_constr_to_1 = np.unique(np.nonzero(P2[info_rows,:])[1])
        # b = b - np.sum(P2[:,variables_constr_to_1],1)
        
        # P2 = P2[keep_rows,:]
        # b = b[keep_rows]

        # #print(variables_constr_to_1)

        # # determined_vars_to_1 = np.unique(np.nonzero(P2[s == b,:])[1])
        # # mask = np.zeros(np.shape(P2)[1], dtype=np.int8)
        # # mask[determined_vars_to_1] = 1

        # # print(determined_vars_to_1)
        # # print(mask)
        # print('P2 =',P2,';')
        # print(f'b = {b}\';')
        # print('s =',np.sum(P2,1),';')
        #print(P2.shape)
        #print(partialInfos)
        #print()
        #print(uncovered)
        raise RuntimeError()
        #assert np.all(P == P.T)

        problem = constraint.Problem()

        for rowidx, row in enumerate(P):
            if np.any(row) == 0:
                continue
            Xs = ['x'+str(i) for i in np.nonzero(row)[0]]
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
                P2 = P[np.ix_(np.sum(P,1) > 0, np.sum(P,0) > 0)]    # tolgo le righe o colonne vuote
                                                            # MATRICE DI ADIACENZA: l'indice di riga indica celle scoperte, quello di colonna celle coperte (vicine a info)
                                                            # Probabilmente è TUM
                assert np.all(ravel_unconvered[np.sum(P,1) > 0])
                print(P2)
                print(ravel_partialInfos[np.sum(P,1) > 0])
                G = nx.from_numpy_matrix(P, create_using=nx.DiGraph)
                nx.draw_networkx(G, pos = {i:np.unravel_index(i,(self.H,self.W), order='F') for i in range(self.H * self.W)})
                plt.matshow(uncovered)
                plt.matshow(partialInfos)
                plt.matshow(P)
                plt.show()
                exit()
                #TODO: potrebbe convenire colpire sul perimetro
                v = np.random.choice(farCells)
                hopeMove = np.unravel_index(v, (self.H, self.W))
                return False, [hopeMove]
                #exit(-1)
        
        

        # print(yolo)
        # print(tallySolutions)
        # print(list(problem._variables.keys()))
        # Quad Prog sembra funzionare min x*(x-1)
        orderedvars = sorted(problem._variables.keys(), key=lambda v: int(v[1:]))
        P2 = P[np.ix_(np.sum(P,1) > 0, np.sum(P,0) > 0)]    # tolgo le righe o colonne vuote
        numVars = P2.shape[1]
        coinvolto = P2.sum(axis = 0)
        assert len(coinvolto) == numVars
        try:
            bella = solve_qp(
                P = 2*np.eye(numVars),
                q = -np.ones(numVars),
                #P = 2.0*np.diag(coinvolto),
                #q = -1.0*coinvolto,
                G = None,
                h = None,
                A = P2,
                b = ravel_partialInfos[np.sum(P,1) > 0],
                lb = np.zeros(numVars),
                ub = np.ones(numVars),
            )
            # bella = cvxopt.solvers.qp(
            #     P=cvxopt.matrix(2*np.eye(numVars).astype(float)),
            #     q=cvxopt.matrix(-np.ones(numVars).astype(float)),
            #     G=cvxopt.matrix(np.vstack([np.eye(numVars),-np.eye(numVars)]).astype(float)),
            #     h=cvxopt.matrix(np.hstack([np.ones(numVars),np.zeros(numVars)]).astype(float)),
            #     A=cvxopt.matrix(P2.astype(float)),
            #     b=cvxopt.matrix(ravel_partialInfos[np.sum(P,1) > 0].astype(float)))['x']
        except Exception:
            self.debug(solutions,tallySolutions,problem,P2,ravel_partialInfos,P,uncovered,partialInfos,orderedvars)
        
        fromenumeration = np.asarray([tallySolutions[v] for v in orderedvars])
        fromqp = (np.asarray(bella).T * len(solutions)).round().astype(int).squeeze()
        if np.all(fromenumeration != 0) and np.all(fromenumeration != numberOfSolutions):
            try:
                assert np.all(fromenumeration == fromqp)
            except Exception:
                print('Errore qp differente da enumerazione')
                print(numberOfSolutions)
                print(fromqp)
                print(fromenumeration)
                self.debug(solutions,tallySolutions,problem,P2,ravel_partialInfos,P,uncovered,partialInfos,orderedvars)
        #print(np.asarray(bella).T * len(solutions))
        #exit()


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

    def debug(self,solutions,tallySolutions,problem,P2,ravel_partialInfos,P,uncovered,partialInfos,orderedvars):
        print('errore')
        print('board')
        print(self.print(uncovered, partialInfos))
        traceback.print_exc()
        print('vars',orderedvars)
        print('numero soluzioni',len(solutions))
        print('tally',np.asarray([tallySolutions[v] for v in orderedvars]))
        print('A =',P2,';')
        print('b =',ravel_partialInfos[np.sum(P,1) > 0],';')
        G = nx.from_numpy_matrix(P, create_using=nx.DiGraph)
        nx.draw_networkx(G, pos = {i:np.unravel_index(i,(self.H,self.W), order='F') for i in range(self.H * self.W)})
        plt.matshow(uncovered)
        plt.matshow(partialInfos)
        plt.matshow(P)
        plt.show()
        exit()