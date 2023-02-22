import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)

class nevillesMethod:
    def __init__(self, xValues, yValues, value):
        self.xValues = np.array(xValues)
        self.yValues = np.array(yValues)
        self.value = value

    def nevillesAlgorithm(self):
        matrix = np.zeros((3, 3))

        for k, row in enumerate(matrix):
             row[0] = self.yValues[k]

        numPoints = len(self.xValues)
        for i in range(1,numPoints):
            for j in range(1, numPoints):
                first_multiplication = (self.value - self.xValues[i]) * matrix[i][j-1]
                second_multiplication = (self.value -self.xValues[i-1]) * matrix[i-1][j-1]

                denominator = self.xValues[i] - self.xValues[i-1]

                coefficient = (second_multiplication - first_multiplication) / denominator
                matrix[i][j] = coefficient
        return matrix[1][1]

class NewtonsDivideDifferences:
    def __init__(self, xValues, yValues, value):
        self.size = len(xValues)
        self.xValues = np.array(xValues)
        self.yValues = np.array(yValues)
        self.value = value
        self.size = len(self.xValues)
        self.matrix = np.zeros((self.size, self.size))
        self.polyCoeff = list(())

        for index, row in enumerate(self.matrix):
             row[0] = self.yValues[index]

        for i in range(1, self.size):
            for j in range(1, i + 1):
                numerator = self.matrix[i][j-1] - self.matrix[i-1][j-1]
                denominator = self.xValues[i] - self.xValues[i - j]
                operation = numerator / denominator
                self.matrix[i][j] = '{0:.7g}'.format(operation)

    def getResults(self):
        reoccuring_x_span = 1
        reoccuring_px_result = self.matrix[0][0]
        polyCoeff = []
        temp = []
       
        for index in range(1, self.size):
            polynomial_coefficient = self.matrix[index][index]
            reoccuring_x_span *= (self.value - self.xValues[index - 1])
            mult_operation = polynomial_coefficient * reoccuring_x_span
            reoccuring_px_result += mult_operation
            self.polyCoeff.append(polynomial_coefficient)    
        return reoccuring_px_result
        
class HermiteInterpolation:
    def __init__(self, xVals, yVals, fPrime):
        self.numPoints = len(xVals)
        self.xVals = xVals
        self.yVals = yVals
        self.fPrime = fPrime
        self.matrix = np.zeros((2*self.numPoints, 2*self.numPoints))

        j = 0
        for x in range(0, 2 * self.numPoints, 2):
            self.matrix[x][0] = self.xVals[j]
            self.matrix[x + 1][0] = self.xVals[j]
            j += 1
        j = 0
        for x in range(0, 2 * self.numPoints, 2):
            self.matrix[x][1] = self.yVals[j]
            self.matrix[x + 1][1] = self.yVals[j]
            j += 1
        j = 0
        for x in range(1, 2 * self.numPoints, 2):
            self.matrix[x][2] = self.fPrime[j]
            if x == 2* self.numPoints - 1:
                break
            self.matrix[x + 1][2] = self.fPrime[j]
            j += 1
        self.matrix[4][2] = -1.18

    def divideDiff(self):
        size = len(self.matrix)

        for i in range(2, size):
            for j in range(2, i + 2):
                if j >= len(self.matrix[i]) or self.matrix[i][j] != 0:
                    continue
                left = self.matrix[i][j - 1]
                diag = self.matrix[i - 1][j - 1]
                numerator = left - diag
                denominator = self.matrix[i][0] - self.matrix[i - 2][0]
                operation = numerator / denominator
                self.matrix[i][j] = operation

class CubicSpline:
    def __init__(self, xVals, yVals):
        self.xVals = np.array(xVals)
        self.yVals = np.array(yVals)
        self.deltaX = np.diff(self.xVals)
        self.deltaY = np.diff(self.yVals)
        self.numPoints = len(self.xVals)
        self.h = [self.xVals[i + 1] - self.xVals[i] for i in range(self.numPoints - 1)]
        self.A = np.zeros((self.numPoints,self.numPoints))
        self.A[0,0] = 1
        self.A[-1,-1] = 1
        self.b = np.zeros((self.numPoints,1))

        for i in range(1, self.numPoints-1):
            self.A[i, i-1] = self.deltaX[i-1]
            self.A[i, i+1] = self.deltaX[i]
            self.A[i,i] = 2*(self.deltaX[i-1]+self.deltaX[i])
            self.b[i,0] = 3*(self.deltaY[i] / self.deltaX[i] - self.deltaY[i-1] / self.deltaX[i-1])

    def solveMatrix(self):
        tolerance = 10e-100
        xDiff = tolerance + 1
        n = self.A.shape[0]
        x = np.zeros(len(self.A))
        xPrev = np.zeros(len(self.A))
        while (xDiff > tolerance):
            for i in range(0, n):
                s = 0
                for j in range(0, n):
                    if i != j:
                        s += self.A[i,j] * xPrev[j] 

                x[i] = (self.b[i] - s) / self.A[i,i]

            xDiff = (np.sum((x - xPrev)**2))**0.5 
            xPrev = x.copy()
        return x
np.set_printoptions(precision=7, suppress=True, linewidth=100)

q1 = nevillesMethod([3.6,3.8,3.9], [1.675, 1.436, 1.318], 3.7)
print("{:.7f}".format(q1.nevillesAlgorithm()))
print("\n") 
q2 = NewtonsDivideDifferences([7.2, 7.4, 7.5, 7.6], [23.5492, 25.3913, 26.8224, 27.4589], 7.3)
q2Ans = q2.getResults()
print(q2.polyCoeff)
print("\n")
print(q2Ans)
print("\n")
q4 = HermiteInterpolation([3.6, 3.8, 3.9], [1.675, 1.436, 1.318], [-1.195, -1.188, -1.182])
q4.divideDiff()
print(q4.matrix)
print("\n")
q5 = CubicSpline([2, 5, 8, 10], [3, 5, 7, 9])
print(q5.A)
print("\n")
print(np.transpose(q5.b))
print("\n")
print(q5.solveMatrix())