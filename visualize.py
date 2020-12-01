from cmu_112_graphics import *
from simpleMolecule import SimpleMolecule
import numpy as np

class MyApp(App):
    # how to use *argv: https://www.geeksforgeeks.org/args-kwargs-python/
    @staticmethod
    def matmul(*matrices):
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
    @staticmethod
    def rotX(theta, M):
        rotMat = np.array([[1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta), np.cos(theta)]])
        return MyApp.matmul(rotMat, M)
    @staticmethod
    def rotY(theta, M):
        rotMat = np.array([[np.cos(theta), 0, np.sin(theta)],
                            [0, 1, 0],
                            [-np.sin(theta), 0, np.cos(theta)]])
        return MyApp.matmul(rotMat, M)
    def appStarted(self):
        self.e1 = np.array([[1],[0],[0]])
        self.e2 = np.array([[0],[1],[0]])
        self.lastX, self.lastY = None, None
        self.offset = np.array([[self.width / 2],[self.height/2],[0]])
        self.vector = np.array([[-100,100],[0,0],[0,0]])
        self.thetaX = 0
        self.thetaY = 0
        self.rotator = np.eye(3)
        self.margin = 20
        self.sliderX = self.width / 2, self.margin
        self.sliderY = self.width - self.margin, self.width / 2
        self.sliderR = 10
        self.molRadius = 5
        self.xDragging = False
        self.yDragging = False
        self.projection = self.vector
    def mousePressed(self, event):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        if (sliderXPosX - self.sliderR <= event.x <= sliderXPosX + self.sliderR) and (sliderXPosY - self.sliderR <= event.y <= sliderXPosY + self.sliderR): self.xDragging = True
        elif (sliderYPosX - self.sliderR <= event.x <= sliderYPosX + self.sliderR) and (sliderYPosY - self.sliderR <= event.y <= sliderYPosY + self.sliderR): self.yDragging = True
    def mouseDragged(self, event):
        sliderXLastX, sliderXLastY = self.sliderX 
        sliderYLastX, sliderYLastY = self.sliderY
        if self.xDragging:
            if 2 * self.margin <= event.x <= self.width - 2 * self.margin: self.sliderX = event.x, sliderXLastY
            elif event.x < 2 * self.margin: self.sliderX = 2 * self.margin, sliderXLastY
            elif event.x > self.width - 2 * self.margin: self.sliderX = self.width - 2 * self.margin, sliderXLastY
        elif self.yDragging:
            if 2 * self.margin <= event.y <= self.height - 2 * self.margin: self.sliderY = sliderYLastX, event.y
            elif 2 * self.margin > event.y: self.sliderY = sliderYLastX, 2 * self.margin
            elif self.height - 2 * self.margin < event.y: self.sliderY = sliderYLastX, self.height - 2 * self.margin
        totalXSliderLength = (self.width - 2 * self.margin) - 2 * self.margin
        totalYSliderLength = (self.height - 2 * self.margin) - 2 * self.margin
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        self.thetaX = ((sliderYPosY - 2 * self.margin) / totalYSliderLength) * 2 * np.pi - np.pi
        self.thetaY = ((sliderXPosX - 2 * self.margin) / totalXSliderLength) * 2 * np.pi - np.pi
        self.rotateObject()
    def mouseReleased(self, event):
        self.xDragging, self.yDragging = False, False
    def rotateObject(self):
        self.projection = MyApp.rotX(self.thetaX, MyApp.rotY(self.thetaY, self.vector))
    def redrawAll(self, canvas):
        vector = self.projection + self.offset
        vec1, vec2 = vector[:,0], vector[:,1]
        self.drawMolecule(canvas, 'ccccc')
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        r = self.sliderR
        self.drawSlider(canvas, sliderXPosX, sliderXPosY, sliderYPosX, sliderYPosY, r)
    def drawMolecule(self, canvas, smiles):
        a = SimpleMolecule(smiles)
        bondVectors = a.bondVectors
        atomVectors = a.atomVectors
        atomVectors = SimpleMolecule.getListOfAtomVectors(atomVectors)
        bondAndAtomVectors = bondVectors + atomVectors
        bondAndAtomVectors = SimpleMolecule.normalizeBondAndAtomVectors(bondAndAtomVectors)
        bondAndAtomVectors = self.rotateBondsAndAtoms(bondAndAtomVectors)
        sortedBondAndAtomVectors = self.sortVectors(bondAndAtomVectors)
        scale = 300
        for vector in sortedBondAndAtomVectors:
            if type(vector) == list:
                vec1, vec2 = vector
                scaleFactor1 = (1 + vec1[2][0] / scale)
                scaleFactor2 = (1 + vec2[2][0] / scale)
                vec1 = vec1 * scaleFactor1
                vec2 = vec2 * scaleFactor2
                vec1 = vec1 + self.offset
                vec2 = vec2 + self.offset
                canvas.create_line(vec1[0][0], vec1[1][0], vec2[0][0], vec2[1][0], fill='black')
            else:
                atom = vector[0]
                vec = vector[1]
                scaleFactor = (1 + vec[2][0] / scale)
                radius = self.molRadius * scaleFactor
                vec = vec * scaleFactor
                if atom == 'c': color = 'red'
                elif atom == 'h': color = 'blue'
                else: color = 'green'
                #vector = MyApp.rotX(self.thetaX, MyApp.rotY(self.thetaY, vector))
                vector = vec + self.offset
                canvas.create_oval((vector[0][0] - radius), (vector[1][0] - radius), (vector[0][0] + radius), (vector[1][0] + radius), fill=color)    
    def rotateBondsAndAtoms(self, bondAndAtomVectors):
        thetaX = self.thetaX
        thetaY = self.thetaY
        result = []
        for vector in bondAndAtomVectors:
            if type(vector) == np.ndarray:
                vec1, vec2 = vector[:, 0], vector[:, 1]
                vec1 = MyApp.rotX(thetaX, MyApp.rotY(thetaY, vec1))
                vec2 = MyApp.rotX(thetaX, MyApp.rotY(thetaY, vec2))
                result.append([vec1, vec2])
            else:
                atom = vector[0]
                vec = vector[1]
                vec = MyApp.rotX(thetaX, MyApp.rotY(thetaY, vec))
                result.append((atom, vec))
        return result
    def sortVectors(self, L): # sorts vectos based on z-value (proximity to observer)
    # from https://www.cs.cmu.edu/~112/notes/notes-recursion-part1.html
        if (len(L) < 2):
            return L
        else:
            mid = len(L)//2
            left = self.sortVectors(L[:mid])
            right = self.sortVectors(L[mid:])
            return self.merge(left, right)
    def merge(self, A, B):
    # heavily inspired by https://www.cs.cmu.edu/~112/notes/notes-recursion-part1.html
        C = [ ]
        i = j = 0
        while ((i < len(A)) or (j < len(B))):
            try:
                if type(A[i]) == list: # i.e. is it a bondVector
                    valInA = min(A[i][0][2][0], A[i][1][2][0])
                else: # i.e. is it a tuple of an atomVector
                    valInA = A[i][1][2][0]
                if type(B[j]) == np.ndarray: # similarly as for A[i], check for B[j]
                    valInB = min(B[j][0][2][0], B[j][1][2][0])
                else:
                    valInB = B[j][1][2][0]
            except: pass
            if ((j == len(B)) or ((i < len(A)) and (valInA <= valInB))):
                C.append(A[i])
                i += 1
            else:
                C.append(B[j])
                j += 1
        return C
    def drawSlider(self, canvas, sliderXPosX, sliderXPosY, sliderYPosX, sliderYPosY, r):
        canvas.create_rectangle(2 * self.margin, self.margin, self.width - 2 * self.margin, self.margin, fill='black')
        canvas.create_rectangle(self.width - self.margin, 2 * self.margin, self.width - self.margin, self.height - 2 * self.margin, fill='black')
        canvas.create_oval(sliderYPosX - r, sliderYPosY - r, sliderYPosX + r, sliderYPosY + r, fill='yellow')
        canvas.create_oval(sliderXPosX - r, sliderXPosY - r, sliderXPosX + r, sliderXPosY + r, fill='yellow')
MyApp(width=400, height=400)