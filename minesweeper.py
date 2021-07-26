import numpy as np
from scipy.signal import convolve2d


class MineSweeper():
    def __init__(self, gameLevel):
        if gameLevel == 'easy':
            self.N = 10
            self.H = self.W = 9
        else:
            raise RuntimeError('Modality not yet implemented')
        self.kernel = np.ones((3, 3), dtype=np.uint8)
        self.kernel[1, 1] = 0
        self.reset()

    def reset(self):
        self.bombs, self.infos = self.generateBoard()
        self.uncovered = np.zeros_like(self.infos)

    def generateBoard(self):
        bombs_map = np.zeros((self.H, self.W), dtype=np.uint8)
        while bombs_map.sum() < self.N:
            h = np.random.randint(0, self.H)
            w = np.random.randint(0, self.W)
            if bombs_map[h, w] == 0:
                bombs_map[h, w] = 1
        info_map = convolve2d(bombs_map, self.kernel, 'same')
        return bombs_map, info_map

    def makeMove(self, r, c):
        if self.uncovered[r, c] == 1:
            raise RuntimeError('Cell already uncovered')
        if self.bombs[r, c] == 1:
            return False
        self.uncover(r, c)

    def uncover(self, r, c):
        if r < 0 or r >= self.W or c < 0 or c >= self.H:
            return
        if self.uncovered[r, c] == 1:
            return
        self.uncovered[r, c] = 1
        if self.infos[r, c] == 0:
            for i in [-1, 0, +1]:
                for j in [-1, 0, +1]:
                    self.uncover(r+i, c+j)

    def __str__(self):
        s = ''
        for r in range(self.H):
            for c in range(self.W):
                if self.uncovered[r, c] == 0:
                    s += chr(9608)
                else:
                    s += str(self.infos[r, c])
                s += ' '
            s += '\n'
        return s

    def getState(self):
        return self.uncovered, self.infos*self.uncovered

    def gameOver(self):
        return np.all(self.uncovered+self.bombs == np.ones_like(self.bombs))
