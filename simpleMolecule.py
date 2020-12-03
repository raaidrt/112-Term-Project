import numpy as np
class SimpleMolecule(object):
    @staticmethod
    def simpleMolecule(s, firstPass = True):
        if s == '': return {}
        elif s == 'h': return {s:[]}
        else:
            if firstPass: restOfMolecule = SimpleMolecule.getRestOfMolecule(s,bonded=False)
            else: restOfMolecule = SimpleMolecule.getRestOfMolecule(s, bonded=True)
            return {s[0]:[(1, SimpleMolecule.simpleMolecule(atom, firstPass=False)) for atom in restOfMolecule]}
    @staticmethod
    def getTrailingMolecule(s):
        s = s[1:]
        while len(s) > 0:
            parensCount = 0
            if s[0] == '(': parensCount += 1
            elif s[0] == 'c' or s[0] == 'n' or s[0] == 's': return [s]
            else: s = s[1:]
            while parensCount != 0:
                s = s[1:]
                if s[0] == '(': parensCount += 1
                elif s[0] == ')': parensCount -= 1
        return []
    @staticmethod
    def getBranch(s):
        result = []
        s = s[1:]
        while len(s) > 0:
            parensCount = 0
            temp = ''
            if s[0] == '(': parensCount += 1
            elif s[0] == 'c' or s[0] == 'n' or s[0] == 's': return result
            while parensCount != 0:
                temp += s[0]
                s = s[1:]
                if s[0] == '(': parensCount += 1
                elif s[0] == ')': parensCount -= 1
            temp = temp[1:]
            if temp != '': result += [temp]
            s = s[1:]
        return result
    @staticmethod
    def getNonBranch(s):
        result = []
        while len(s) > 0:
            parensCount = 0
            if s[0] == '(': parensCount += 1
            elif s[0] == 'c' or s[0] == 'n' or s[0] == 's': return result
            else: result += [s[0]]
            while parensCount != 0:
                s = s[1:]
                if s[0] == '(': parensCount += 1
                elif s[0] == ')': parensCount -= 1
            s = s[1:]
        return result
    @staticmethod
    def getAtomVectors(s):
        molecule = SimpleMolecule.simpleMolecule(s)
        chainLength = SimpleMolecule.getChainLength(molecule)
        if s == '': return
        atom = s[0]
        if atom != 'c': raise Exception('Only carbon chains are supported')
        vectorAtChainPosition = np.array([0,0,0])
        result = {atom:[vectorAtChainPosition]}
        restOfMolecule = SimpleMolecule.getRestOfMolecule(s)
        firstAtom = restOfMolecule.pop(0)
        firstAtomListOfVectors = result.get(firstAtom, [])
        tetrahedralAngle = (109.5 / 180) * np.pi
        firstAtomVector = np.array([50 * np.cos((np.pi / 2) + tetrahedralAngle / 2), 50 * np.sin((np.pi / 2) + tetrahedralAngle / 2), 0])
        firstAtomListOfVectors += [firstAtomVector]
        result[firstAtom] = firstAtomListOfVectors
        for chainPosition in range(chainLength):
            secondAtom = restOfMolecule[0]
            thirdAtom = restOfMolecule[1]
            fourthAtom = restOfMolecule[2]
            if chainPosition % 2 == 0:
                secondAtomVector = np.array([0, -50 * np.cos(tetrahedralAngle / 2), np.sin(tetrahedralAngle / 2)])
                thirdAtomVector = np.array([0, -50 * np.cos(tetrahedralAngle / 2), -np.sin(tetrahedralAngle / 2)])
                fourthAtomVector  = np.array([50 * np.cos((np.pi / 2) - (tetrahedralAngle / 2)), 50 * np.sin((np.pi / 2) - (tetrahedralAngle / 2)), 0])
            else:
                secondAtomVector = np.array([0, 50 * np.cos(tetrahedralAngle / 2), np.sin(tetrahedralAngle / 2)])
                thirdAtomVector = np.array([0, 50 * np.cos(tetrahedralAngle / 2), -np.sin(tetrahedralAngle / 2)])
                fourthAtomVector = np.array([50 * np.cos((3 * np.pi / 2) + (tetrahedralAngle / 2)), 50 * np.sin((3 * np.pi / 2) + (tetrahedralAngle / 2)), 0])
            result[secondAtom] = result.get(secondAtom, []) + [vectorAtChainPosition + secondAtomVector]
            result[thirdAtom] = result.get(thirdAtom, []) + [vectorAtChainPosition + thirdAtomVector]
            fourthAtom = fourthAtom[0]
            result[fourthAtom] = result.get(fourthAtom, []) + [vectorAtChainPosition + fourthAtomVector]
            vectorAtChainPosition = vectorAtChainPosition + fourthAtomVector
            atom = restOfMolecule[-1][0]
            restOfMolecule = SimpleMolecule.getRestOfMolecule(restOfMolecule[-1], bonded=True)
        return result
    @staticmethod
    def getAtomAndBondVectors(s):
        molecule = SimpleMolecule.simpleMolecule(s)
        chainLength = SimpleMolecule.getChainLength(molecule)
        if s == '': return
        atom = s[0]
        if atom != 'c': raise Exception('Only carbon chains are supported')
        vectorAtChainPosition = np.array([0,0,0])
        atomVectors = {atom:[vectorAtChainPosition]}
        bondVectors = []
        restOfMolecule = SimpleMolecule.getRestOfMolecule(s)
        firstAtom = restOfMolecule.pop(0)
        firstAtomListOfVectors = atomVectors.get(firstAtom, [])
        tetrahedralAngle = (109.5 / 180) * np.pi
        firstAtomVector = np.array([50 * np.cos((np.pi / 2) + tetrahedralAngle / 2), 50 * np.sin((np.pi / 2) + tetrahedralAngle / 2), 0])
        radiusVector = np.array([5 * np.cos((np.pi / 2) + tetrahedralAngle / 2), 5 * np.sin((np.pi / 2) + tetrahedralAngle / 2), 0])
        firstAtomBondVector = firstAtomVector - radiusVector
        vectorAtChainPositionBondVector = vectorAtChainPosition + radiusVector
        firstAtomListOfVectors += [firstAtomVector]
        atomVectors[firstAtom] = firstAtomListOfVectors
        bondVectors.append(np.array([vectorAtChainPositionBondVector, firstAtomBondVector]).T)
        for chainPosition in range(chainLength):
            secondAtom = restOfMolecule[0]
            thirdAtom = restOfMolecule[1]
            fourthAtom = restOfMolecule[2]
            if chainPosition % 2 == 0:
                secondAtomVector = np.array([0, -50 * np.cos(tetrahedralAngle / 2), 50 * np.sin(tetrahedralAngle / 2)])
                secondRadiusVector = np.array([0, -5 * np.cos(tetrahedralAngle / 2), 5 * np.sin(tetrahedralAngle / 2)])
                thirdAtomVector  = np.array([0, -50 * np.cos(tetrahedralAngle / 2), -50 * np.cos(tetrahedralAngle / 2)])
                thirdRadiusVector  = np.array([0, -5 * np.cos(tetrahedralAngle / 2), -5 * np.cos(tetrahedralAngle / 2)])
                fourthAtomVector  = np.array([50 * np.cos((np.pi / 2) - (tetrahedralAngle / 2)), 50 * np.sin((np.pi / 2) - (tetrahedralAngle / 2)), 0])
                fourthRadiusVector = np.array([5 * np.cos((np.pi / 2) - (tetrahedralAngle / 2)), 5 * np.sin((np.pi / 2) - (tetrahedralAngle / 2)), 0])
            else:
                secondAtomVector = np.array([0, 50 * np.cos(tetrahedralAngle / 2), 50 * np.cos(tetrahedralAngle / 2)])
                secondRadiusVector = np.array([0, 5 * np.cos(tetrahedralAngle / 2), 5 * np.cos(tetrahedralAngle / 2)])
                thirdAtomVector = np.array([0, 50 * np.cos(tetrahedralAngle / 2), -50 * np.cos(tetrahedralAngle / 2)])
                thirdRadiusVector  = np.array([0, 5 * np.cos(tetrahedralAngle / 2), -5 * np.cos(tetrahedralAngle / 2)])
                fourthAtomVector = np.array([50 * np.cos((3 * np.pi / 2) + (tetrahedralAngle / 2)), 50 * np.sin((3 * np.pi / 2) + (tetrahedralAngle / 2)), 0])
                fourthRadiusVector = np.array([5 * np.cos((3 * np.pi / 2) + (tetrahedralAngle / 2)), 5 * np.sin((3 * np.pi / 2) + (tetrahedralAngle / 2)), 0])
            atomVectors[secondAtom] = atomVectors.get(secondAtom, []) + [vectorAtChainPosition + secondAtomVector]
            atomVectors[thirdAtom] = atomVectors.get(thirdAtom, []) + [vectorAtChainPosition + thirdAtomVector]
            fourthAtom = fourthAtom[0]
            atomVectors[fourthAtom] = atomVectors.get(fourthAtom, []) + [vectorAtChainPosition + fourthAtomVector]
            bondVectors.append(np.array([secondAtomVector - secondRadiusVector + vectorAtChainPosition, vectorAtChainPosition + secondRadiusVector]).T)
            bondVectors.append(np.array([thirdAtomVector - thirdRadiusVector + vectorAtChainPosition, vectorAtChainPosition + thirdRadiusVector]).T)
            bondVectors.append(np.array([fourthAtomVector - fourthRadiusVector + vectorAtChainPosition, vectorAtChainPosition + fourthRadiusVector]).T)
            vectorAtChainPosition = vectorAtChainPosition + fourthAtomVector
            atom = restOfMolecule[-1][0]
            restOfMolecule = SimpleMolecule.getRestOfMolecule(restOfMolecule[-1], bonded=True)
        return (atomVectors, bondVectors)
    @staticmethod
    def getChainLength(molecule):
        if molecule == {}: return 0
        for atom in molecule:
            if molecule[atom] == [ ]: return 0
            else: return 1 + SimpleMolecule.getChainLength(molecule[atom][-1][1])
    @staticmethod
    def getRestOfMolecule(s, bonded = False):
        atom = s[0]
        branch = SimpleMolecule.getBranch(s)
        nonBranch = SimpleMolecule.getNonBranch(s)
        domains = 0
        if atom == 'c': 
            domains = 4
        elif atom == 'n': 
            domains = 3
        trailingMolecule = SimpleMolecule.getTrailingMolecule(s)
        hydrogens = domains - len(nonBranch) - len(branch) - len(trailingMolecule)
        if bonded: hydrogens -= 1
        restOfMolecule = ['h' for i in range(hydrogens)] + branch + nonBranch + trailingMolecule
        return restOfMolecule
    @staticmethod
    def getListOfAtomVectors(d):
        result = [ ]
        for atom in d:
            for vector in d[atom]:
                result.append((atom, vector))
        return result
    @staticmethod
    def normalizeBondAndAtomVectors(vectors):  
        count, average = 0, np.array([0,0,0])
        for vector in vectors:
            if not type(vector) == np.ndarray:
                average = average + vector[1]
                count += 1
        average = average / count
        result = [ ]
        for vector in vectors:
            if type(vector) == np.ndarray:
                result.append((vector.T - average).T)
            else:
                result.append((vector[0], vector[1] - average))
        return result
    def __init__(self, s):
        self.smiles = s
        self.molecule = SimpleMolecule.simpleMolecule(s)
        atomAndBondVectors = SimpleMolecule.getAtomAndBondVectors(s)
        self.atomVectors = atomAndBondVectors[0]
        self.bondVectors = atomAndBondVectors[1]
    def __repr__(self):
        return self.smiles
    def __getitem__(self, start, stop, step):
        pass