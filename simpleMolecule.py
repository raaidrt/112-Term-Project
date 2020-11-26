class SimpleMolecule(object):
    @staticmethod
    def simpleMolecule(s, d={}):
        if s == '':
            return d
        else:
            s1, s2 = s[0], s[1:]
            if s1 == 'c':
                if d == {}:
                    return {(0,'C*'):[(1,'H'), \
                    (1,SimpleMolecule.simpleMolecule(s2)), (1,'H')]}
    @staticmethod
    def moleculeToStr(d, s = ''):
        if d == {}:
            return s
        else:
            for bond, atom in d:
                if s == '': s += atom[0]
                else:
                    rest = d[(bond, atom)]
    def __init__(self, s):
        self.molecule = SimpleMolecule.simpleMolecule(s)
    def __repr__(self):
        return (str(self.molecule))