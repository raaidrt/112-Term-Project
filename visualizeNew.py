from cmu_112_graphics import *
import numpy as np
class MyApp(App):
    # how to use *argv:
    # https://www.geeksforgeeks.org/args-kwargs-python/
    @staticmethod
    def matmul(*argv):
        result = None
        for M in argv:
            n, p = M.shape
            if type(result) != np.ndarray:
                result = M
            else:
                m, n = result.shape
                temp = np.array([[0] * p for i in range(m)])
                for i in range(m):
                    for j in range(p):
                        temp[i,j] = sum([result[i,k] * M[k,j] for k in range(n)])
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
        self.xDragging = False
        self.yDragging = False
        self.projection = self.vector
    def mousePressed(self, event):
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        if (sliderXPosX - self.sliderR <= event.x <= sliderXPosX + self.sliderR) and (sliderXPosY - self.sliderR <= event.y <= sliderXPosY + self.sliderR):
            self.xDragging = True
        elif (sliderYPosX - self.sliderR <= event.x <= sliderYPosX + self.sliderR) and (sliderYPosY - self.sliderR <= event.y <= sliderYPosY + self.sliderR):
            self.yDragging = True
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
        self.rotateObject()
    def mouseReleased(self, event):
        self.xDragging, self.yDragging = False, False
    def rotateObject(self):
        totalXSliderLength = (self.width - 2 * self.margin) - 2 * self.margin
        totalYSliderLength = (self.height - 2 * self.margin) - 2 * self.margin
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        thetaX = ((sliderYPosY - 2 * self.margin) / totalYSliderLength) * 2 * np.pi - np.pi
        thetaY = ((sliderXPosX - 2 * self.margin) / totalXSliderLength) * 2 * np.pi - np.pi
        self.projection = MyApp.rotX(thetaX, MyApp.rotY(thetaY, self.vector))
    def redrawAll(self, canvas):
        vector = self.projection + self.offset
        vec1, vec2 = vector[:,0], vector[:,1]
        print(vec1, vec2)
        canvas.create_line(vec1[0], vec1[1], vec2[0], vec2[1])
        r1 = 50 * (100 + vec1[2]) / 500
        r2 = 50 * (100 + vec2[2]) / 500
        if vec1[2] < vec2[2]:
            canvas.create_oval(vec1[0] - r1, vec1[1] - r1, vec1[0] + r1, vec1[1] + r1, fill="red", width=2 * (100 + vec1[2]) / 1000)
            canvas.create_oval(vec2[0] - r2, vec2[1] - r2, vec2[0] + r2, vec2[1] + r2, fill="green", width=2 * (100 + vec1[2]) / 1000)
        else:
            canvas.create_oval(vec2[0] - r2, vec2[1] - r2, vec2[0] + r2, vec2[1] + r2, fill="green", width=2 * (100 + vec1[2]) / 1000)
            canvas.create_oval(vec1[0] - r1, vec1[1] - r1, vec1[0] + r1, vec1[1] + r1, fill="red", width=2 * (100 + vec1[2]) / 1000)
        sliderXPosX, sliderXPosY = self.sliderX
        sliderYPosX, sliderYPosY = self.sliderY
        r = self.sliderR
        self.drawSlider(canvas, sliderXPosX, sliderXPosY, sliderYPosX, sliderYPosY, r)
    def drawSlider(self, canvas, sliderXPosX, sliderXPosY, sliderYPosX, sliderYPosY, r):
        canvas.create_rectangle(2 * self.margin, self.margin, self.width - 2 * self.margin, self.margin, fill='black')
        canvas.create_rectangle(self.width - self.margin, 2 * self.margin, self.width - self.margin, self.height - 2 * self.margin, fill='black')
        canvas.create_oval(sliderYPosX - r, sliderYPosY - r, sliderYPosX + r, sliderYPosY + r, fill='yellow')
        canvas.create_oval(sliderXPosX - r, sliderXPosY - r, sliderXPosX + r, sliderXPosY + r, fill='yellow')
MyApp(width=400, height=400)