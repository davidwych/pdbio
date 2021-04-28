# pdbio.py
# written by David Wych (2020) -- UCI and LANL
# Convenience script for reading a writing PDB files
# and manipulating atom and crystallographic info
# minimilist and unpolished -- use at you peril

class Atom():
    ''' Object used to store data from each PDB line starting with
        ATOM or HETATM
    '''
    def __init__(self, line=None):
        if line == None:
            pass

        else:
            self.index    = int(line[6:11])
            self.atomtype = line[12:16].strip()
            self.resname  = line[17:20].strip()
            self.chainid  = line[21]
            self.resid    = int(line[22:27])
            self.x        = float(line[30:38])
            self.y        = float(line[38:46])
            self.z        = float(line[46:54])
            self.occ      = float(line[54:60])
            self.bfac     = float(line[60:66])
            self.element  = line[76:78].strip()
            if len(self.element) == 0:
                if self.atomtype[0] in ['1','2','3']:
                    self.element = self.atomtype[1]
                else:
                    self.element = self.atomtype[0]

class CrystInfo():
    ''' Object used to store data from CRYST1 record
    '''
    def __init__(self, line=None):
        if line == None:
            pass
        else:
            self.a      = float(line[6:15])
            self.b      = float(line[15:24])
            self.c      = float(line[24:33])
            self.alpha  = float(line[33:40])
            self.beta   = float(line[40:47])
            self.gamma  = float(line[47:54])
            self.sgroup = line[55:66]
            self.zval   = line[66:70]

class PDBFile():
    ''' Object used to store all the data from a PDB file
    '''
    def __init__(self, ifilename=None, crystinfo=None, contents=None):
        if ifilename == None:
            self.filename = None
            self.crystinfo = crystinfo
            self.contents = contents
            if "TER\n" not in self.contents:
                self.contents.append("TER\n")
        else:
            self.filename = ifilename
            with open(self.filename, 'r') as f:
                lines = f.readlines()
                self.contents = []
                for line in lines:
                    if line[:6] == "CRYST1":
                        self.crystinfo = CrystInfo(line)
                    if line[:6] in ["ATOM  ", "HETATM"]:
                        self.contents.append(Atom(line))
                    if line[:3] == "TER":
                        self.contents.append("TER\n")
                if "TER\n" not in self.contents:
                    self.contents.append("TER\n")

        if len(self.contents) != 0:
            self.atoms = [el for el in self.contents if not isinstance(el, str)]
            if len(self.atoms) == 0:
                self.n_residues = 0
                self.n_atoms    = 0
            else:
                _res = self.atoms[0].resid; self.n_residues = 1
                for atom in self.atoms:
                    if atom.resid != _res:
                        _res = atom.resid
                        self.n_residues += 1

                self.n_atoms = len(self.atoms)

        else:
            self.n_residues = 0
            self.n_atoms    = 0

    def __repr__ (self):
        if self.filename == None:
            return "{} atoms; {} residues".format(len(self.atoms), self.n_residues)
        else:
            return "{}: {} atoms; {} residues".format(self.filename, len(self.atoms), self.n_residues)

    def combine(self, other):
        _contents = self.contents + other.contents
        self.contents = _contents
        _n_residues = self.n_residues + other.n_residues
        self.n_residues = _n_residues
        return self

    def n_atoms(self):
        return len(self.atoms)

    def renumber_residues(self):
        _res = self.atoms[0].resname; _r = 0; _atom = self.atoms[0].atomtype
        for atom in self.atoms:
            if (atom.resname == _res) and (atom.atomtype != _atom):
                atom.resid = _r
            elif (atom.resname == _res) and (atom.atomtype == _atom):
                _r += 1
                atom.resid = _r
            elif (atom.resname != _res):
                _r += 1
                _atom = atom.atomtype
                atom.resid = _r
                _res = atom.resname

    def make_model(self, num):
        self.contents = ["MODEL     {}\n".format(num)] + self.contents + ["ENDMDL\n"]
        return self

    def write(self, ofilename):
        with open(ofilename, 'w') as of:
            of.write("CRYST1{:>9}{:>9}{:>9}{:>7}{:>7}{:>7} {:>11}{:>4}\n".format("{:.3f}".format(self.crystinfo.a),
                                                                                 "{:.3f}".format(self.crystinfo.b),
                                                                                 "{:.3f}".format(self.crystinfo.c),
                                                                                 "{:.2f}".format(self.crystinfo.alpha),
                                                                                 "{:.2f}".format(self.crystinfo.beta),
                                                                                 "{:.2f}".format(self.crystinfo.gamma),
                                                                                                 self.crystinfo.sgroup,
                                                                                                 self.crystinfo.zval))
            for el in self.contents:
                if isinstance(el, str):
                    of.write(el)

                else:
                    if len(el.atomtype) <= 3:
                        el.atomtype = " " + el.atomtype

                    of.write("ATOM  {:>5} {:<4}{:>4} {}{:>4}    {:>8}{:>8}{:>8}{:>6}{:>6}          {:>2}  \n".format(el.index,
                                                                                                                     el.atomtype,
                                                                                                                     el.resname,
                                                                                                                     el.chainid,
                                                                                                                     el.resid,
                                                                                                                     "{:.3f}".format(el.x),
                                                                                                                     "{:.3f}".format(el.y),
                                                                                                                     "{:.3f}".format(el.z),
                                                                                                                     "{:.2f}".format(el.occ),
                                                                                                                     "{:.2f}".format(el.bfac),
                                                                                                                     el.element))
