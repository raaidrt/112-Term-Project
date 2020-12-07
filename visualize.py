from cmu_112_graphics import *
from simpleMolecule import SimpleMolecule
import numpy as np
from matrixOperations import *
class SimpleMoleculeViewAndShell(App):    
    # appStarted to initialization all values and states
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
        self.buttonMargins = 2
        self.splashX, self.splashY = None, None
        self.showSplash = False
        self.splashMolecule = None
        self.splashDims = 40
        self.zScale = 300
        self.scalar = 1
        self.error = None
        self.cellInError = False
        self.errorSplashCoords = (None, None)
        self.showErrorSplash = False
    # helper method for appStarted
    def initializeWindowDims(self, molViewWidth, molViewHeight):
        self.margin = 20
        self.molViewWidth, self.molViewHeight = molViewWidth, molViewHeight
        self.shellWidth, self.shellHeight = self.width - self.molViewHeight, self.height - self.margin
        self.molOptionsWidth, self.molOptionsHeight = self.molViewWidth, self.height - self.molViewHeight
    # helper method for appStarted
    def initializeVectorStuff(self):
        self.molRadius = 5
        self.offset = np.array([[self.molViewWidth / 2],[self.molViewHeight/2],[0]])
    # helper method for appStarted
    def intiializeSliderStates(self):
        self.thetaX = 0
        self.thetaY = 0
        self.sliderX = self.molViewWidth / 2, self.margin
        self.sliderY = self.molViewWidth - self.margin, self.molViewWidth / 2
        self.sliderZoom = self.molViewWidth / 2, self.molViewHeight - self.margin
        self.sliderR = 10
        self.xDragging = False
        self.yDragging = False
        self.zoomDragging = False
    # for when the mouse is pressed on the canvas
    def mousePressed(self, event):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        sliderZoomPosX, sliderZoomPosY = self.sliderZoom
        self.shellModel(event)
        if (sliderXPosX - self.sliderR <= event.x <= sliderXPosX + self.sliderR) and (sliderXPosY - self.sliderR <= event.y <= sliderXPosY + self.sliderR): self.xDragging = True
        elif (sliderYPosX - self.sliderR <= event.x <= sliderYPosX + self.sliderR) and (sliderYPosY - self.sliderR <= event.y <= sliderYPosY + self.sliderR): self.yDragging = True
        elif (sliderZoomPosX - self.sliderR <= event.x <= sliderZoomPosX + self.sliderR) and (sliderZoomPosY - self.sliderR <= event.y <= sliderZoomPosY + self.sliderR): self.zoomDragging = True
        self.moleculeModel(event)
        self.clearShellButtonModel(event)
        self.clearMemoryButtonModel(event)
    # helper for interacting with the clearMemoryButton
    def clearMemoryButtonModel(self, event):
        clearButtonHeight = 3 * self.margin - 2 * self.buttonMargins
        clearButtonWidth = self.shellWidth / 2 - 1.5 * self.buttonMargins
        x1 = self.width - self.buttonMargins
        x0 = x1 - clearButtonWidth
        y0 = self.height - 3 * self.margin + self.buttonMargins
        y1 = y0 + clearButtonHeight
        if x0 <= event.x <= x1 and y0 <= event.y <= y1: 
            self.molecules = {}
            self.moleculeInView = None
    # for when the mouse hovers over a particular atom or cell in the shell
    def mouseMoved(self, event):
        if self.cellInError and self.coordsOnCurrentCell(event.x, event.y):
            self.errorSplashCoords = event.x, event.y
        else: self.errorSplashCoords = None, None
        if self.moleculeInView == None: return
        atomVectors = self.molecules[self.moleculeInView].atomVectors
        atomVectors = SimpleMolecule.getListOfAtomVectors(atomVectors)
        atomVectors = SimpleMolecule.normalizeBondAndAtomVectors(atomVectors)
        atomVectors = self.rotateBondsAndAtoms(atomVectors)
        sortedAtomVectors = self.sortVectors(atomVectors)[::-1]
        scale = self.zScale
        for vector in sortedAtomVectors:
            vector = vector
            atom = vector[0]
            vec = vector[1]
            scaleFactor = (1 + vec[2][0] / scale) * self.scalar
            r = self.molRadius * scaleFactor
            vec = vec * scaleFactor
            vec = vec + self.offset
            x, y = vec[0][0], vec[1][0]
            if x - r <= event.x <= x + r and y - r <= event.y <= y + r:
                self.splashX, self.splashY = event.x, event.y
                self.showSplash = True
                self.splashMolecule = atom
                return
        self.splashX, self.splashY = None, None
        self.showSplash = False
        self.splashMolecule = None
    # helper function to check if coords are on top of current cell (bottom-most cell in the shell)
    def coordsOnCurrentCell(self, x, y):
        cellHeight = (self.shellHeight - 2 * self.margin) / self.maxCellsInShell
        lastIndex = len(self.commands)
        x0, y0, x1, y1 = self.getShellMargins(lastIndex, cellHeight)
        return x0 <= x <= x1 and y0 <= y <= y1
    # for when the mouse is pressed then dragged over the canvas -> used for the sliders
    def mouseDragged(self, event):
        self.sliderModel(event)
    # helper for mouseDragged
    def sliderModel(self, event):
        sliderXLastX, sliderXLastY = self.sliderX 
        sliderYLastX, sliderYLastY = self.sliderY
        sliderZoomLastX, sliderZoomLastY = self.sliderZoom
        if self.xDragging:
            if 2 * self.margin <= event.x <= self.molViewWidth - 2 * self.margin: self.sliderX = event.x, sliderXLastY
            elif event.x < 2 * self.margin: self.sliderX = 2 * self.margin, sliderXLastY
            elif event.x > self.molViewWidth - 2 * self.margin: self.sliderX = self.molViewWidth - 2 * self.margin, sliderXLastY
        elif self.yDragging:
            if 2 * self.margin <= event.y <= self.molViewHeight - 2 * self.margin: self.sliderY = sliderYLastX, event.y
            elif 2 * self.margin > event.y: self.sliderY = sliderYLastX, 2 * self.margin
            elif self.molViewHeight - 2 * self.margin < event.y: self.sliderY = sliderYLastX, self.molViewHeight - 2 * self.margin
        elif self.zoomDragging:
            if 2 * self.margin <= event.x <= self.molViewWidth - 2 * self.margin: self.sliderZoom = event.x, sliderZoomLastY
            elif event.x < 2 * self.margin: self.sliderZoom = 2 * self.margin, sliderZoomLastY
            elif event.x > self.molViewWidth - 2 * self.margin: self.sliderZoom = self.molViewWidth - 2 * self.margin, sliderZoomLastY
        totalXSliderLength = (self.molViewWidth - 2 * self.margin) - 2 * self.margin
        totalYSliderLength = (self.molViewHeight - 2 * self.margin) - 2 * self.margin
        totalZoomSliderLength = totalXSliderLength
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        sliderZoomPosX, sliderZoomPosY = self.sliderZoom
        self.thetaX = ((sliderYPosY - 2 * self.margin) / totalYSliderLength) * 2 * np.pi - np.pi
        self.thetaY = ((sliderXPosX - 2 * self.margin) / totalXSliderLength) * 2 * np.pi - np.pi
        self.scalar = ((sliderZoomPosX - 2 * self.margin) / totalZoomSliderLength + 0.5) ** 2
    # for manipulating the molecule view model
    def moleculeModel(self, event):
        count = 0
        for molecule in self.molecules:
            if count >= self.maxCols * self.maxRows: pass
            x0, y0, x1, y1 = self.getCellBoundsOfMolCell(count)
            if x0 <= event.x <= x1 and y0 <= event.y <= y1:
                self.moleculeInView = molecule
            count += 1
    # for manipulating the shell model upon interacting with the clearShell button
    def clearShellButtonModel(self, event):
        clearButtonHeight = 3 * self.margin - 2 * self.buttonMargins
        clearButtonWidth = self.shellWidth / 2 - 1.5 * self.buttonMargins
        x0, y0 = self.buttonMargins + self.molViewWidth, self.height - 3 * self.margin + self.buttonMargins
        x1, y1 = x0 + clearButtonWidth, y0 + clearButtonHeight
        if x0 <= event.x <= x1 and y0 <= event.y <= y1: self.commands = []
    # helper to get getBounds of molecule boxes below the molView
    def getCellBoundsOfMolCell(self, count):
        row = count // self.maxRows
        col = count % self.maxCols
        cellWidth, cellHeight = (self.molOptionsWidth - 2 * self.margin) / self.maxCols, (self.molOptionsHeight - 2 * self.margin) / self.maxRows
        x0, y0 = self.margin + col * cellWidth, self.margin + row * cellHeight
        x1, y1 = x0 + cellWidth, y0 + cellHeight
        offset = self.height - self.molOptionsHeight
        return x0, y0 + offset, x1, y1 + offset
    # for manipulating the shell model
    def shellModel(self, event):
        cellHeight = cellHeight = (self.shellHeight - 2 * self.margin) / self.maxCellsInShell
        lastIndex = len(self.commands)
        x0, y0, x1, y1 = self.getShellMargins(lastIndex, cellHeight)
        self.cellHiLited = (x0 <= event.x <= x1 and y0 <= event.y <= y1)
    # for interacting with the shell by pressing keys
    def keyPressed(self, event):
        if event.key in ['Up','Down','Left','Right','Space','Tab']: return
        if self.cellHiLited:
            if event.key == 'Enter':
                try: self.evaluate(self.currentCommand)
                except: 
                    self.evaluateError(self.currentCommand)
                    self.cellInError = True
                    return
                self.commands.append(self.currentCommand)
                if len(self.commands) >= self.maxCellsInShell - 1: self.commands.pop(0)
                self.currentCommand = ">   "
            elif event.key == 'Delete':
                if len(self.currentCommand) > 4:
                    self.currentCommand = self.currentCommand[:-1]
                self.cellInError = False
            else:
                self.currentCommand += event.key
    # evaluating strings of commands entered in the shell
    def evaluate(self, s):
        testStr = s[4:]
        for sign in ['==','>=','<=','>','<','!=']:
            if sign in testStr and '->' not in testStr:
                op1 = testStr.split(sign)[0]
                op2 = testStr.split(sign)[1]
                if sign == '==': 
                    assert(self.molecules[op1] == self.molecules[op2])
                elif sign == '>': 
                    assert(self.molecules[op1] > self.molecules[op2])
                elif sign == '>=':
                    assert(self.molecules[op1] >= self.molecules[op2])
                elif sign == '<':
                    assert(self.molecules[op1] < self.molecules[op2])
                elif sign == '<=': 
                    assert(self.molecules[op1] <= self.molecules[op2])
                elif sign == '!=': 
                    assert(self.molecules[op1] != self.molecules[op2])
                return
        if '->' in s:
            L = s.split('->')
            varName = L[0][4:]
            if '[' in L[1]:
                variable = L[1].split('[')[0]
                args = [_ for _ in L[1].split('[')[1].split(',')]
                args = [L[1].split('[')[1].split(',')[0], L[1].split('[')[1].split(',')[1][:-1]]
                args = [int(_) for _ in args]
                if len(args) == 1: start = end = args[0]
                else: start, end = args
                molecule = self.molecules[variable].molecule
                returnValue = SimpleMolecule(SimpleMolecule.iloc(molecule, start, end))
            elif '+' in L[1]:
                operands = L[1].split('+')
                op1 = self.molecules[operands[0]].smiles
                op2 = self.molecules[operands[1]].smiles
                returnValue = SimpleMolecule(op1 + op2)
            else:
                returnValue = SimpleMolecule(L[1])
            self.molecules[varName] = returnValue
        else: raise Exception()
    # evaluating command strings that contain some error in the shell
    def evaluateError(self, s):
        testStr = s[4:]
        for sign in ['==','>=','<=','>','<','!=']:
            if sign in testStr and '->' not in testStr:
                self.error = "Error: assertion error while\nusing relational operators\ni.e. assertion was False or syntax error. \nPress Backspace to try again!"
                return
        if '->' in s:
            L = s.split('->')
            varName = L[0][4]
            if '[' in L[1]:
                self.error = "Error: when using indices\nuse only pre-set variable names and\nindex numbers separated by commas\nPress Backspace and try agian!"
            elif '+' in L[1]:
                self.error = "Error: the '+' operator must\ntake operands that are previously\nstored variable names\nPress Backspace and try again!"
            else:
                self.error = "Error: when assigning molecules,\nplease only use 'c' and 'h' atoms,\nor put atoms in parentheses\nto make it a branch.\nPress Backspace and try again!"
        else:
            self.error = "Error: must use assignment operator\n'->' somewhere in command. \nPress Backspace and try again!"
    # for when the mouse has been released
    def mouseReleased(self, event):
        self.xDragging, self.yDragging, self.zoomDragging = False, False, False
    # canonical redrawAll function
    def redrawAll(self, canvas):
        self.drawMolView(canvas)
        self.drawSlider(canvas)
        self.drawMolOptionsView(canvas)
        self.drawShell(canvas)
        self.drawClearShellButton(canvas)
        self.drawSplash(canvas)
        self.drawClearMemoryButton(canvas)
    # helper to draw clearMemoryButton
    def drawClearMemoryButton(self, canvas):
        clearButtonHeight = 3 * self.margin - 2 * self.buttonMargins
        clearButtonWidth = self.shellWidth / 2 - 1.5 * self.buttonMargins
        x1 = self.width - self.buttonMargins
        x0 = x1 - clearButtonWidth
        y0 = self.height - 3 * self.margin + self.buttonMargins
        y1 = y0 + clearButtonHeight
        canvas.create_rectangle(x0, y0, x1, y1, fill="lightGreen")
        canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, font=f"Arial {int(clearButtonHeight - 2 * self.buttonMargins)}", text="Clear Memory")
    # helper to draw the molView
    def drawMolView(self, canvas):
        for molecule in self.molecules:
            if molecule == self.moleculeInView:
                self.drawMolecule(canvas, self.molecules[molecule].smiles)
    # helper to draw molecule in the molView
    def drawMolecule(self, canvas, smiles):
        a = SimpleMolecule(smiles)
        bondVectors = a.bondVectors
        atomVectors = a.atomVectors
        atomVectors = SimpleMolecule.getListOfAtomVectors(atomVectors)
        bondAndAtomVectors = bondVectors + atomVectors
        bondAndAtomVectors = SimpleMolecule.normalizeBondAndAtomVectors(bondAndAtomVectors)
        bondAndAtomVectors = self.rotateBondsAndAtoms(bondAndAtomVectors)
        sortedBondAndAtomVectors = self.sortVectors(bondAndAtomVectors)
        scale = self.zScale
        for vector in sortedBondAndAtomVectors:
            if type(vector) == list:
                vec1, vec2 = vector
                scaleFactor1 = (1 + vec1[2][0] / scale) * self.scalar
                scaleFactor2 = (1 + vec2[2][0] / scale) * self.scalar
                vec1 = vec1 * scaleFactor1
                vec2 = vec2 * scaleFactor2
                vec1 = vec1 + self.offset
                vec2 = vec2 + self.offset
                canvas.create_line(vec1[0][0], vec1[1][0], vec2[0][0], vec2[1][0], fill='black')
            else:
                atom = vector[0]
                vec = vector[1]
                scaleFactor = (1 + vec[2][0] / scale) * self.scalar
                radius = self.molRadius * scaleFactor
                vec = vec * scaleFactor
                if atom == 'c': color = 'red'
                elif atom == 'h': color = 'blue'
                else: color = 'green'
                #vector = rotX(self.thetaX, rotY(self.thetaY, vector))
                vector = vec + self.offset
                canvas.create_oval((vector[0][0] - radius), (vector[1][0] - radius), (vector[0][0] + radius), (vector[1][0] + radius), fill=color)    
    # helper to rotate the bonds and atoms by the angles specified
    def rotateBondsAndAtoms(self, bondAndAtomVectors):
        thetaX = self.thetaX
        thetaY = self.thetaY
        result = []
        for vector in bondAndAtomVectors:
            if type(vector) == np.ndarray:
                vec1, vec2 = vector[:, 0], vector[:, 1]
                vec1 = rotX(thetaX, rotY(thetaY, vec1))
                vec2 = rotX(thetaX, rotY(thetaY, vec2))
                result.append([vec1, vec2])
            else:
                atom = vector[0]
                vec = vector[1]
                vec = rotX(thetaX, rotY(thetaY, vec))
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
    def merge(self, A, B): # mergeSort modified from CMU 112 website
    # heavily inspired by https://www.cs.cmu.edu/~112/notes/notes-recursion-part1.html
        C = [ ]
        i = j = 0
        while ((i < len(A)) or (j < len(B))):
            try: # <- my contribution starts here
                if type(A[i]) == list: # i.e. is it a bondVector
                    valInA = min(A[i][0][2][0], A[i][1][2][0])
                else: # i.e. is it a tuple of an atomVector
                    valInA = A[i][1][2][0]
                if type(B[j]) == np.ndarray: # similarly as for A[i], check for B[j]
                    valInB = min(B[j][0][2][0], B[j][1][2][0])
                else:
                    valInB = B[j][1][2][0]
            except: pass
            if ((j == len(B)) or ((i < len(A)) and (valInA <= valInB))): # <- my contribution ends here
                C.append(A[i])
                i += 1
            else:
                C.append(B[j])
                j += 1
        return C
    # draws molOptionsView (below molView)
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
    # helper to draw the sliders
    def drawSlider(self, canvas):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        sliderZoomPosX, sliderZoomPosY = self.sliderZoom
        r = self.sliderR
        canvas.create_rectangle(2 * self.margin, self.margin, self.molViewWidth - 2 * self.margin, self.margin, fill='black')
        canvas.create_rectangle(self.molViewWidth - self.margin, 2 * self.margin, self.molViewWidth - self.margin, self.molViewHeight - 2 * self.margin, fill='black')
        canvas.create_rectangle(2 * self.margin, self.molViewHeight - self.margin, self.molViewWidth - 2 * self.margin, self.molViewHeight - self.margin, fill='black')
        canvas.create_oval(sliderYPosX - r, sliderYPosY - r, sliderYPosX + r, sliderYPosY + r, fill='yellow')
        canvas.create_oval(sliderXPosX - r, sliderXPosY - r, sliderXPosX + r, sliderXPosY + r, fill='yellow')
        canvas.create_oval(sliderZoomPosX - r, sliderZoomPosY - r, sliderZoomPosX + r, sliderZoomPosY + r, fill='cyan')
    # helper to draw the shell and all the cells in it
    def drawShell(self, canvas):
        cellHeight = (self.shellHeight - 2 * self.margin) / self.maxCellsInShell
        count = 0
        insideCellMargin = 5
        for command in self.commands:
            x0, y0, x1, y1 = self.getShellMargins(count, cellHeight)
            canvas.create_rectangle(x0, y0, x1, y1)
            canvas.create_text(x0 + insideCellMargin, y0 + insideCellMargin, text=command, font=f'Arial {int(cellHeight - 2 * insideCellMargin)}', anchor='nw')
            count += 1
        lastIndex = len(self.commands)
        x0, y0, x1, y1 = self.getShellMargins(lastIndex, cellHeight)
        if self.cellHiLited: color = 'lightGreen'
        else: color = 'white'
        if self.cellInError: color = 'red'
        canvas.create_rectangle(x0, y0, x1, y1, fill=color)
        currentCommand = self.currentCommand
        canvas.create_text(x0 + insideCellMargin, y0 + insideCellMargin, text=self.currentCommand, font=f'Arial {int(cellHeight - 2 * insideCellMargin)}', anchor='nw')
        if self.cellInError and self.errorSplashCoords != (None, None):
            x, y = self.errorSplashCoords
            x0, y0 = x - 250, y
            x1, y1 = x, y + 100
            canvas.create_rectangle(x0, y0, x1, y1, fill="orange")
            canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=self.error)
    # helper to get shellBounds
    def getShellMargins(self, count, cellHeight):
        x1 = self.width - self.margin - self.cellInShellMargin
        x0 = self.molViewWidth + self.margin + self.cellInShellMargin
        y0 = self.margin + self.cellInShellMargin + count * cellHeight
        y1 = y0 + cellHeight + self.cellInShellMargin
        return x0, y0, x1, y1
    # helper to draw the clearShellButton
    def drawClearShellButton(self, canvas):
        clearButtonHeight = 3 * self.margin - 2 * self.buttonMargins
        clearButtonWidth = self.shellWidth / 2 - 1.5 * self.buttonMargins
        x0, y0 = self.buttonMargins + self.molViewWidth, self.height - 3 * self.margin + self.buttonMargins
        x1, y1 = x0 + clearButtonWidth, y0 + clearButtonHeight
        canvas.create_rectangle(x0, y0, x1, y1, fill="red")
        canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, font=f"Arial {int(clearButtonHeight - 2 * self.buttonMargins)}", text="Clear Shell")
    # helper that draws the splash screens (when hovering on top of the atoms)
    def drawSplash(self, canvas):
        if self.showSplash:
            x0, y0 = self.splashX, self.splashY
            x1, y1 = x0 + 3 * self.splashDims, y0 + self.splashDims
            canvas.create_rectangle(x0, y0, x1, y1, fill='lightGreen')
            canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=f"{self.splashMolecule.upper()}")