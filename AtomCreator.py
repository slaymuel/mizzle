from pyquaternion import Quaternion
import numpy as np
import random

#Returns rotational matrix which rotates around x-axis
def xRotate(vector, angle):
    rotMatrix = [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]
    dotProd = np.dot(rotMatrix, vector)
    return (dotProd)

#Returns rotational matrix that randomly rotates around given vector
def random_rotate(vector, angle):
    rotMatrix = np.array([[np.cos(angle)+vector[0]**2*(1-np.cos(angle)), vector[0]*vector[1]*(1 - np.cos(angle))-vector[2]*np.sin(angle), vector[0]*vector[2]*(1-np.cos(angle))+vector[1]*np.sin(angle)],
                [vector[1]*vector[0]*(1 - np.cos(angle))+vector[2]*np.sin(angle), np.cos(angle)+vector[1]**2*(1-np.cos(angle)), vector[1]*vector[2]*(1-np.cos(angle))-vector[0]*np.sin(angle)],
                [vector[2]*vector[0]*(1 - np.cos(angle))-vector[1]*np.sin(angle), vector[2]*vector[1]*(1 - np.cos(angle))+vector[0]*np.sin(angle), np.cos(angle)+vector[2]**2*(1-np.cos(angle))]])
    return rotMatrix

#Calculate rotation matrix by rotating z-axis (0, 0, 1) to align with MO-vector
def align(vec1, vec2):
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    crossProd = np.cross(vec1, vec2)    #sine
    dotProd = np.dot(vec1, vec2)    #cosine

    #Edge cases: If vectors are parallell since dot product of normalized vectors is the cosine of the angle.
    if(dotProd < 1.01 and dotProd > 0.99):
        rotMatrix = I

    elif(dotProd < -0.99 and dotProd > -1.01):
        #print("dotprod is -1")
        #Find orthonormal vector to both (cross product) and rotate pi
        #mag = np.sqrt(crossProd[0]**2 + crossProd[1]**2 + crossProd[2]**2)
        #ortVec = crossProd/mag

        #Will always be the same? Are there edge cases.....?
        #rotMatrix = -np.array([[-1+2*ortVec[0]**2, 2*ortVec[0]*ortVec[1], 2*ortVec[0]*ortVec[2]],
        #                   [2*ortVec[0]*ortVec[1], -1+2*ortVec[1]**2, 2*ortVec[1]*ortVec[2]],
        #                   [2*ortVec[0]*ortVec[2], 2*ortVec[1]*ortVec[2], -1+2*ortVec[2]**2]])
        rotMatrix = -I      #np.array([[1, 0, 0],[0, -1, 0],[0, 0, -1]])

    #Need to take into account when vectors are orthagonal i.e when dot product is 0?
    else:
        #skew-symmetric: transpose equals negative. Used to represent cross product as matrix multiplication
        skewSym = np.array([[0, -crossProd[2], crossProd[1]],[crossProd[2], 0, -crossProd[0]],[-crossProd[1], crossProd[0], 0]])
        prod = np.matmul(skewSym, skewSym)*(1/(1+dotProd))
        rotMatrix = np.add(np.add(I, skewSym), prod)

    return rotMatrix

def add_hydroxyl(coords, vectors, theta):
    atoms = np.empty([0, 3], dtype=float)
    elements = []
    i=0
    #Loop over all coordinates where oxygen is to be placed
    while i < len(coords):
        O = np.array([0, 0, 0])
        H = np.array([np.sin(theta/2), 0, np.cos(theta/2)])

        #No need to rotate O since it lies on the x-axis
        angle = np.arccos(np.cos(np.radians(115))/np.cos(np.radians(104.5/2)))
        H = xRotate(H, angle)

        #Align z axis to the directional vector
        rotMatrix = align([0, 0, 1], vectors[i])
        O = np.dot(rotMatrix, O)
        H = np.dot(rotMatrix, H)



        #Random rotation along directional vector
        #randRotMatrix = self.randomRot(vectors[i], random.uniform(0.1, 2*np.pi))
        #randRotMatrix = self.randomRot(vectors[i], 1.2)

        #axis around which to rotate the molecule
        axis = vectors[i]
        randAngle = random.uniform(0.1, 2*np.pi)
        #O = Quaternion(axis=axis,angle=theta).rotate(O)
        H = Quaternion(axis=axis,angle=randAngle).rotate(H)

        #O = randRotMatrix.dot(O)
        #H1 = randRotMatrix.dot(H1)
        #H1 = randRotMatrix.dot(H2)

        #Translate to correct coordinates
        transVector = [coords[i][0] - O[0], coords[i][1] - O[1], coords[i][2] - O[2]]
        O = np.array([O[0] + transVector[0], O[1] + transVector[1], O[2] + transVector[2]])
        H = np.array([H[0] + transVector[0], H[1] + transVector[1], H[2] + transVector[2]])

        i += 1


        #Save atoms to be added to pdb file
        atoms = np.vstack((atoms, O))
        elements.extend('O')
        atoms = np.vstack((atoms, H))
        elements.extend('H')
    return atoms, elements

def add_water(coords, vectors, theta):
    atoms = np.empty([0, 3], dtype=float)
    elements = []
    i=0
    #Loop over all coordinates where oxygen is to be placed
    while i < len(coords):
        O = np.array([0, 0, 0])
        H1 = np.array([np.sin(theta/2), 0, np.cos(theta/2)])
        H2 = np.array([-np.sin(theta/2), 0, np.cos(theta/2)])

        #No need to rotate O since it lies on the x-axis
        angle = np.arccos(np.cos(115)/np.cos(104.5/2))
        H1 = xRotate(H1, angle)
        H2 = xRotate(H2, angle)

        #Align z axis to the directional vector
        rotMatrix = align([0, 0, 1], vectors[i])
        O = np.dot(rotMatrix, O)
        H1 = np.dot(rotMatrix, H1)
        H2 = np.dot(rotMatrix, H2)



        #Random rotation along directional vector
        #randRotMatrix = self.randomRot(vectors[i], random.uniform(0.1, 2*np.pi))
        #randRotMatrix = self.randomRot(vectors[i], 1.2)

        axis = vectors[i]
        randAngle = random.uniform(0.1, 2*np.pi)
        #O = Quaternion(axis=axis,angle=theta).rotate(O)
        H1 = Quaternion(axis=axis,angle=randAngle).rotate(H1)
        H2 = Quaternion(axis=axis,angle=randAngle).rotate(H2)

        #O = randRotMatrix.dot(O)
        #H1 = randRotMatrix.dot(H1)
        #H1 = randRotMatrix.dot(H2)

        #Translate to correct coordinates
        transVector = [coords[i][0] - O[0], coords[i][1] - O[1], coords[i][2] - O[2]]
        O = np.array([O[0] + transVector[0], O[1] + transVector[1], O[2] + transVector[2]])
        H1 = np.array([H1[0] + transVector[0], H1[1] + transVector[1], H1[2] + transVector[2]])
        H2 = np.array([H2[0] + transVector[0], H2[1] + transVector[1], H2[2] + transVector[2]])

        i += 1


        #Save atoms to be added to pdb file
        atoms = np.vstack((atoms, O))
        elements.extend('O')
        atoms = np.vstack((atoms, H1))
        elements.extend('H')
        atoms = np.vstack((atoms, H2))
        elements.extend('H')

    numOfOverlaps = 0
    i = 0
    j = 0
    while i < len(atoms):
        j = i + 1
        while j < len(atoms):
            if(overlap(atoms[i], atoms[j])):
                numOfOverlaps += 1
            j += 1
        i += 1
    numOfOverlaps = numOfOverlaps-3*len(atoms)/3
    print(str(numOfOverlaps) + " overlapping atoms")

    return atoms, elements


def overlap(point1, point2):
    vector = np.array([point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]])
    distance = np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    if(distance < 2):
        return True

    else:
        return False