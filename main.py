from minesweeper import MineSweeper
from solver import Solver
import numpy as np

ms = MineSweeper('easy')
solver = Solver(H = 9, W = 9, N = 10)

def play(rStart,cStart, TOT = 1000):
    vinte = 0
    perse = 0
    for game in range(TOT):
        ms.reset()
        while ms.makeMove(rStart,cStart) == False:
            ms.reset()
        loose = False
        while (not ms.gameOver()) and (not loose):
            uncovered, partialInfos = ms.getState()
            safe, moves = solver.getMove(uncovered, partialInfos)
            for (r,c) in moves:
                if ms.uncovered[r,c] == 0:
                    res = ms.makeMove(r,c)
                    if res == False:
                        loose = True
                        assert safe == False
                        break
        if ms.gameOver():
            vinte += 1
        else:
            perse += 1

    print(vinte,perse,vinte/(vinte+perse),perse/(vinte+perse))
    return vinte

# res = np.zeros((9,9))
# for r in range(9):
#     for c in range(9):
#         res[r,c] = play(r,c,TOT=10)

play(0,0,TOT = 1000)

