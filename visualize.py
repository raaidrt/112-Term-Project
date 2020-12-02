from cmu_112_graphics import *
from simpleMolecule import SimpleMolecule
import numpy as np
class MyApp(App):
    # how to use *argv: https://www.geeksforgeeks.org/args-kwargs-python/
    @staticmethod
    def matmul(*matrices):
        # for multiplying matrices together, note args are specified in order
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
        # matrices from https://en.wikipedia.org/wiki/Rotation_matrix
        rotMat = np.array([[1, 0, 0],
                            [0, np.cos(theta), -np.sin(theta)],
                            [0, np.sin(theta), np.cos(theta)]])
        return MyApp.matmul(rotMat, M)
    @staticmethod
    def rotY(theta, M):
        # matrices from https://en.wikipedia.org/wiki/Rotation_matrix
        rotMat = np.array([[np.cos(theta), 0, np.sin(theta)],
                            [0, 1, 0],
                            [-np.sin(theta), 0, np.cos(theta)]])
        return MyApp.matmul(rotMat, M)
    def appStarted(self):
        self.initializeWindowDims(400,400)
        self.intiializeSliderStates()
        self.initializeVectorStuff()
        self.moleculeInView = None
        self.maxRows, self.maxCols = 5, 5
        self.maxCellsInShell = 9
        self.cellInShellMargin = 3
        self.molecules = {}
        self.commands = []
        self.currentCommand = ">   "
        self.cellHiLited = False
    def initializeWindowDims(self, molViewWidth, molViewHeight):
        self.molViewWidth, self.molViewHeight = molViewWidth, molViewHeight
        self.shellWidth, self.shellHeight = self.width - self.molViewHeight, self.height
        self.molOptionsWidth, self.molOptionsHeight = self.molViewWidth, self.height - self.molViewHeight
    def initializeVectorStuff(self):
        self.molRadius = 5
        self.offset = np.array([[self.molViewWidth / 2],[self.molViewHeight/2],[0]])
    def intiializeSliderStates(self):
        self.thetaX = 0
        self.thetaY = 0
        self.margin = 20
        self.sliderX = self.molViewWidth / 2, self.margin
        self.sliderY = self.molViewWidth - self.margin, self.molViewWidth / 2
        self.sliderR = 10
        self.xDragging = False
        self.yDragging = False
    def mousePressed(self, event):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        self.shellModel(event)
        if (sliderXPosX - self.sliderR <= event.x <= sliderXPosX + self.sliderR) and (sliderXPosY - self.sliderR <= event.y <= sliderXPosY + self.sliderR): self.xDragging = True
        elif (sliderYPosX - self.sliderR <= event.x <= sliderYPosX + self.sliderR) and (sliderYPosY - self.sliderR <= event.y <= sliderYPosY + self.sliderR): self.yDragging = True
        self.moleculeModel(event)
    def mouseDragged(self, event):
        self.sliderModel(event)
    def sliderModel(self, event):
        sliderXLastX, sliderXLastY = self.sliderX 
        sliderYLastX, sliderYLastY = self.sliderY
        if self.xDragging:
            if 2 * self.margin <= event.x <= self.molViewWidth - 2 * self.margin: self.sliderX = event.x, sliderXLastY
            elif event.x < 2 * self.margin: self.sliderX = 2 * self.margin, sliderXLastY
            elif event.x > self.molViewWidth - 2 * self.margin: self.sliderX = self.molViewWidth - 2 * self.margin, sliderXLastY
        elif self.yDragging:
            if 2 * self.margin <= event.y <= self.molViewHeight - 2 * self.margin: self.sliderY = sliderYLastX, event.y
            elif 2 * self.margin > event.y: self.sliderY = sliderYLastX, 2 * self.margin
            elif self.molViewHeight - 2 * self.margin < event.y: self.sliderY = sliderYLastX, self.molViewHeight - 2 * self.margin
        totalXSliderLength = (self.molViewWidth - 2 * self.margin) - 2 * self.margin
        totalYSliderLength = (self.molViewHeight - 2 * self.margin) - 2 * self.margin
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        self.thetaX = ((sliderYPosY - 2 * self.margin) / totalYSliderLength) * 2 * np.pi - np.pi
        self.thetaY = ((sliderXPosX - 2 * self.margin) / totalXSliderLength) * 2 * np.pi - np.pi
    def moleculeModel(self, event):
        count = 0
        for molecule in self.molecules:
            if count >= self.maxCols * self.maxRows: pass
            x0, y0, x1, y1 = self.getCellBoundsOfMolCell(count)
            if x0 <= event.x <= x1 and y0 <= event.y <= y1:
                self.moleculeInView = molecule
            count += 1
    def getCellBoundsOfMolCell(self, count):
        row = count // self.maxRows
        col = count % self.maxCols
        cellWidth, cellHeight = (self.molOptionsWidth - 2 * self.margin) / self.maxCols, (self.molOptionsHeight - 2 * self.margin) / self.maxRows
        x0, y0 = self.margin + col * cellWidth, self.margin + row * cellHeight
        x1, y1 = x0 + cellWidth, y0 + cellHeight
        offset = self.height - self.molOptionsHeight
        return x0, y0 + offset, x1, y1 + offset
    def shellModel(self, event):
        cellHeight = cellHeight = (self.shellHeight - 2 * self.margin) / self.maxCellsInShell
        lastIndex = len(self.commands)
        x0, y0, x1, y1 = self.getShellMargins(lastIndex, cellHeight)
        self.cellHiLited = (x0 <= event.x <= x1 and y0 <= event.y <= y1)
    def keyPressed(self, event):
        if self.cellHiLited:
            if event.key == 'Enter':
                self.evaluate(self.currentCommand)
                self.commands.append(self.currentCommand)
                if len(self.commands) >= self.maxCellsInShell - 1: self.commands.pop(0)
                self.currentCommand = ">   "
            elif event.key == 'Delete':
                self.currentCommand = self.currentCommand[:-1]
            else:
                self.currentCommand += event.key
    def evaluate(self, s):
        L = s.split('->')
        self.molecules[L[0][4:]] = SimpleMolecule(L[1])
    def mouseReleased(self, event):
        self.xDragging, self.yDragging = False, False
    def redrawAll(self, canvas):
        self.drawSlider(canvas)
        self.drawMolOptionsView(canvas)
        self.drawMolView(canvas)
        self.drawShell(canvas)
    def drawMolView(self, canvas):
        for molecule in self.molecules:
            if molecule == self.moleculeInView:
                self.drawMolecule(canvas, self.molecules[molecule].smiles)
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
    def sortVectors(self, L): # sorts vectors based on z-value (proximity to observer)
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
    def drawMolOptionsView(self, canvas):
        count = 0
        for molecule in self.molecules:
            x0, y0, x1, y1 = self.getCellBoundsOfMolCell(count)
            color = "white"
            if molecule == self.moleculeInView: color = "lightGreen"
            canvas.create_rectangle(x0, y0, x1, y1, fill=color)
            canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=molecule)
            count += 1
            if count >= self.maxRows * self.maxCols: break
    def drawSlider(self, canvas):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        r = self.sliderR
        canvas.create_rectangle(2 * self.margin, self.margin, self.molViewWidth - 2 * self.margin, self.margin, fill='black')
        canvas.create_rectangle(self.molViewWidth - self.margin, 2 * self.margin, self.molViewWidth - self.margin, self.molViewHeight - 2 * self.margin, fill='black')
        canvas.create_oval(sliderYPosX - r, sliderYPosY - r, sliderYPosX + r, sliderYPosY + r, fill='yellow')
        canvas.create_oval(sliderXPosX - r, sliderXPosY - r, sliderXPosX + r, sliderXPosY + r, fill='yellow')
    def drawShell(self, canvas):
        cellHeight = (self.shellHeight - 2 * self.margin) / self.maxCellsInShell
        count = 0
        for command in self.commands:
            x0, y0, x1, y1 = self.getShellMargins(count, cellHeight)
            canvas.create_rectangle(x0, y0, x1, y1)
            canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=command, font=f'Montserrat {int(cellHeight - 10)}')
            count += 1
        lastIndex = len(self.commands)
        x0, y0, x1, y1 = self.getShellMargins(lastIndex, cellHeight)
        if self.cellHiLited: color = 'lightGreen'
        else: color = 'white'
        canvas.create_rectangle(x0, y0, x1, y1, fill=color)
        currentCommand = self.currentCommand
        canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=self.currentCommand, font=f'Montserrat {int(cellHeight - 10)}')
    def getShellMargins(self, count, cellHeight):
        x1 = self.width - self.margin - self.cellInShellMargin
        x0 = self.molViewWidth + self.margin + self.cellInShellMargin
        y0 = self.margin + self.cellInShellMargin + count * cellHeight
        y1 = y0 + cellHeight + self.cellInShellMargin
        return x0, y0, x1, y1

MyApp(width=1200, height=600)