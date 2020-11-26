from cmu_112_graphics import *
import numpy as np
def appStarted(app):
    app.l1 = np.array([[app.width // 2 - 50, app.height // 2 - 50], [app.width // 2 + 50, app.height // 2 + 50]])
    app.plane = np.array([[1,0,0],[0,1,0],[0,0,0]])
    app.mouseDragging = False
    app.lastX, app.lastY = None, None
def mousePressed(app, event):
    app.mouseDragging = not app.mouseDragging
def mouseDragged(app, event):
    if app.mouseDragging:
        event.x, event.y = x, y
        if (app.lastX, app.lastY) == (None, None):
            app.lastX, app.lastY = x, y
            return
        else:
            dx, dy = x - app.lastX, y - app.lastY
            thetaX, thetaY =  np.arctan(dx), np.arctan(dy)
            # rotation matrices from https://en.wikipedia.org/wiki/Rotation_matrix
            rotationX = np.array([1,0,0],[0,np.cos(thetaX),-np.sin(thetaX)],[0,np.sin(thetaX),np.cos(thetaX)])
            rotationY = np.array([np.cos(thetaY),0,np.sin(thetaY)],[0,1,0],[-np.sin(thetaY),0,np.cos(thetaY)],)
            app.plane = matmult([rotationX, rotationY, app.plane])
    else:
        return
def matmult(L, result = None):
    if L == []:
        return result
    else:
        M1, M2 = result, L[0] # M1 is m x n, M2 is n x p
        n, p = M2.shape
        if M1 == None: M1 = np.eye(n)
        m, n = M1.shape
        shape = (m, p)
        result = np.array([[0] * shape[1] for i in range(shape[0])])
        for i in range(len(M1)):
            for j in range(len(M2.T)):
                result[i,j] = sum([M1[i,k] * M2[k,j] for k in range(n)])
        return matmult(L[1:], result)
def redrawAll(app, canvas):
    lx1, ly1, lx2, ly2 = app.l1[0][0], app.l1[0][1], app.l1[1][0], app.l1[1][1]
    vec = np.array([[lx1, lx2],[ly1, ly2]])
    A = app.plane
    proj = matmult([A, np.linalg.inv(matmult(A.T, A)), A.T, vec])
    v1, v2 = proj[:,0], proj[:,1]
    canvas.create_line(v1[0,0], v1[1,0], v2[0,0], v2[0,1])
runApp(width=800, height=400)