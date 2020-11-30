import numpy as np
from visualize import *
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
                secondAtomVector = np.array([0, -50, np.sin(tetrahedralAngle / 2)])
                thirdAtomVector  = np.array([0, -50, -np.sin(tetrahedralAngle / 2)])
                fourthAtomVector  = np.array([50 * np.cos((np.pi / 2) - (tetrahedralAngle / 2)), 50 * np.sin((np.pi / 2) - (tetrahedralAngle / 2)), 0])
            else:
                secondAtomVector = np.array([0, 50, np.sin(tetrahedralAngle / 2)])
                thirdAtomVector = np.array([0, 50, -np.sin(tetrahedralAngle / 2)])
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
    def getBondVectors(s):
        pass
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
    def __init__(self, s):
        self.smiles = s
        self.molecule = SimpleMolecule.simpleMolecule(s)
        self.bondVectors = SimpleMolecule.getBondVectors(s)
        self.atomVectors = SimpleMolecule.getAtomVectors(s)
    def __repr__(self):
        return (str(self.molecule))