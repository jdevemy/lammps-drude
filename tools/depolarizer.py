#!/usr/bin/env python
# depolarizer.py - remove Drude oscillators to LAMMPS data file.

"""
Add Drude oscillators to LAMMPS data file.

The format of file containing specification of Drude oscillators is:

# type  dm/u  dq/e  k/(kJ/molA2)  alpha/A3  thole
C3H     1.0   0.0   4184.0        2.051     2.6
...

where:
dm is the mass to put on the Drude particle (taken from its core),
dq is the charge to put on the Drude particle (taken from its core),
k is the garmonic force constant of the bond between Drude and core,
alpha is the polarizability,
thole is a Thole damping parameter.

A Drude particle is created for each atom in the LAMMPS data file
that corresponds to an atom type given in the Drude file.
Since LAMMPS uses numbers for atom types in the data file, a comment
after each line in the Masses section has to be introduced to allow
identification of the atom types within the force field database:

Masses

   1   12.011  # C3H
   2   12.011  # CTO
...


This script will add new atom types, new bond types, new atoms and
new bonds to the data file. It will also add a new seciton called
"Drude Types" containing a flag 0 (non-polarizable atom), 1 (Drude core)
or 2 (Drude particle).

"""

import argparse


# keywords of header and main sections (from data.py in Pizza.py)

hkeywords = ["atoms", "ellipsoids", "lines", "triangles", "bodies",
             "bonds", "angles", "dihedrals", "impropers",
             "atom types", "bond types", "angle types", "dihedral types",
             "improper types", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]

skeywords = [["Masses", "atom types"],
             ["Drude Types", "atom types"],
             ["Pair Coeffs", "atom types"],
             ["Bond Coeffs", "bond types"], ["Angle Coeffs", "angle types"],
             ["Dihedral Coeffs", "dihedral types"],
             ["Improper Coeffs", "improper types"],
             ["BondBond Coeffs", "angle types"],
             ["BondAngle Coeffs", "angle types"],
             ["MiddleBondTorsion Coeffs", "dihedral types"],
             ["EndBondTorsion Coeffs", "dihedral types"],
             ["AngleTorsion Coeffs", "dihedral types"],
             ["AngleAngleTorsion Coeffs", "dihedral types"],
             ["BondBond13 Coeffs", "dihedral types"],
             ["AngleAngle Coeffs", "improper types"],
             ["Atoms", "atoms"], ["Ellipsoids", "ellipsoids"],
             ["Lines", "lines"], ["Triangles", "triangles"],
             ["Bodies", "bodies"],
             ["Bonds", "bonds"],
             ["Angles", "angles"], ["Dihedrals", "dihedrals"],
             ["Impropers", "impropers"], ["Velocities", "atoms"],
             ["Molecules", "atoms"]]


def massline(att):
    return "{0:4d} {1:8.3f}  # {2}\n".format(att['id'], att['m'], att['type'])

def bdtline(bdt):
    return "{0:4d} {1}\n".format(bdt['id'], bdt['note'])

def atomline(at):
    return "{0:7d} {1:7d} {2:4d} {3:8.4f} {4:13.6e} {5:13.6e} {6:13.6e} "\
           "{7}\n".format(at['n'], at['mol'], at['id'], at['q'],
                          at['x'], at['y'], at['z'], at['note'])

def bondline(bd):
    return "{0:7d} {1:4d} {2:7d} {3:7d}  # {4}\n".format(bd['n'], bd['id'],
                                            bd['i'], bd['j'], bd['note'])


# --------------------------------------


class Data(object):

    def __init__(self, datafile):
        '''read LAMMPS data file (from data.py in Pizza.py)'''

        # for extract method
        self.atomtypes = []
        self.bondtypes = []
        self.atoms = []
        self.bonds = []

        self.nselect = 1

        f = open(datafile, "r")

        self.title = f.readline()
        self.names = {}

        headers = {}
        while 1:
            line = f.readline().strip()
            if '#' in line:
                line = line[:line.index('#')].strip()
            if len(line) == 0:
                continue
            found = 0
            for keyword in hkeywords:
                if keyword in line:
                    found = 1
                    words = line.split()
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                      keyword == "zlo zhi":
                        headers[keyword] = (float(words[0]), float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = \
                          (float(words[0]), float(words[1]), float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
            if not found:
                break

        sections = {}
        while 1:
            if len(line) > 0:
                found = 0
                for pair in skeywords:
                    keyword, length = pair[0], pair[1]
                    if keyword == line:
                        found = 1
                        if length not in headers:
                            raise RuntimeError("data section {} has no matching"\
                                               " header value".format(line))
                        f.readline()
                        list_data = []
                        for _ in range(headers[length]):
                            list_data.append(f.readline())
                        sections[keyword] = list_data
                if not found:
                    raise RuntimeError("invalid section {} in data"\
                                       " file".format(line))
            #f.readline()
            line = f.readline()
            if not line:
                break
            if '#' in line:
                line = line[:line.index('#')]
            line = line.strip()

        f.close()
        self.headers = headers
        self.sections = sections

    def write(self, filename):
        '''write out a LAMMPS data file (from data.py in Pizza.py)'''

        with open(filename, "w") as f:
            f.write(self.title + '\n')
            for keyword in hkeywords:
                if keyword in self.headers:
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                       keyword == "zlo zhi":
                        pair = [ str(p) for p in self.headers[keyword] ]
                        f.write(pair[0] + ' ' + pair[1] + ' ' + keyword + '\n')
                    elif keyword == "xy xz yz":
                        triple = [ str(t) for t in self.headers[keyword] ]
                        f.write(triple[0] + ' ' + triple[1] + ' ' + triple[2] +
                                ' ' + keyword + '\n')
                    else:
                        f.write(str(self.headers[keyword]) + ' ' +
                                keyword + '\n')
            for pair in skeywords:
                keyword = pair[0]
                if keyword == "Drude Types":
                    continue
                if keyword in self.sections:
                    f.write("\n{}\n\n".format(keyword))
                    for line in self.sections[keyword]:
                        f.write(line)

    def extract(self):
        """extract atom IDs, bond types and atoms from data"""

        # extract atom IDs
        for line in self.sections['Masses']:
            tok = line.split()
            if len(tok) < 4:
                print("warning: missing type for atom ID " + tok[0])
                continue
            atomtype = {}
            atomtype['id'] = int(tok[0])
            atomtype['m'] = float(tok[1])
            atomtype['type'] = tok[3]
            self.atomtypes.append(atomtype)

        # extract atom drude types
        for itype, line in enumerate(self.sections['Drude Types']):
            tok = line.split()
            if len(tok) < 2:
                print("warning: missing Drude Type for atom type " + tok[0])
                continue
            if self.atomtypes[itype]['id'] != int(tok[0]):
                print("warning: Drude Types or Masses are in wrong order")
            self.atomtypes[itype]['dflag'] = int(tok[1])

        # extract bond type data
        for line in self.sections['Bond Coeffs']:
            tok = line.split()
            bondtype = {}
            bondtype['id'] = int(tok[0])
            bondtype['note'] = ' '.join(tok[1:])
            self.bondtypes.append(bondtype)

        # extract atom registers
        for line in self.sections['Atoms']:
            tok = line.split()
            atom = {}
            atom['n'] = int(tok[0])
            atom['mol'] = int(tok[1])
            atom['id'] = int(tok[2])
            atom['q'] = float(tok[3])
            atom['x'] = float(tok[4])
            atom['y'] = float(tok[5])
            atom['z'] = float(tok[6])
            note = ''
            for i in range(7, len(tok)):
                note += tok[i] + ' '
            atom['note'] = note.strip()
            self.atoms.append(atom)

        # extract bond data
        for line in self.sections['Bonds']:
            tok = line.split()
            bond = {}
            bond['n'] = int(tok[0])
            bond['id'] = int(tok[1])
            bond['i'] = int(tok[2])
            bond['j'] = int(tok[3])
            bond['note'] = ' '.join(tok[4:])
            self.bonds.append(bond)

    def depolarize(self):
        """remove Drude particles"""

        self.extract()

        atom_tp_map = {}
        bond_tp_map = {}
        atom_map = {}
        bond_map = {}
        q = {}
        atom_tp = {}
        m = {}

        for att in self.atomtypes:
            if att['dflag'] != 2:
                atom_tp_map[att['id']] = len(atom_tp_map) + 1
            m[att['id']] = att['m']
        for atom in self.atoms:
            if atom['id'] in atom_tp_map:
                atom_map[atom['n']] = len(atom_map) + 1
            q[atom['n']] = atom['q']
            atom_tp[atom['n']] = atom['id']
        for bond in self.bonds:
            if bond['i'] in atom_map and bond['j'] in atom_map:
                bond_map[bond['n']] = len(bond_map) + 1
                if bond['id'] not in bond_tp_map:
                    bond_tp_map[bond['id']] = len(bond_tp_map) + 1
            else:
                if bond['i'] in atom_map:
                    q[bond['i']] += q[bond['j']]
                    if atom_tp[bond['j']] in m:
                        m[atom_tp[bond['i']]] += m.pop(atom_tp[bond['j']])
                else:
                    q[bond['j']] += q[bond['i']]
                    if atom_tp[bond['i']] in m:
                        m[atom_tp[bond['j']]] += m.pop(atom_tp[bond['i']])

        size = len(self.atomtypes)
        for iatom_tp in reversed(xrange(size)):
            att = self.atomtypes[iatom_tp]
            if att['id'] not in atom_tp_map:
                del self.atomtypes[iatom_tp]
            else:
                att['m'] = m[att['id']]
                att['id'] = atom_tp_map[att['id']]

        size = len(self.bondtypes)
        for ibond_tp in reversed(xrange(size)):
            bdt = self.bondtypes[ibond_tp]
            if bdt['id'] not in bond_tp_map:
                del self.bondtypes[ibond_tp]
            else:
                bdt['id'] = bond_tp_map[bdt['id']]

        size = len(self.atoms)
        for iatom in reversed(xrange(size)):
            atom = self.atoms[iatom]
            if atom['n'] not in atom_map:
                del self.atoms[iatom]
            else:
                atom['q'] = q[atom['n']]
                atom['n'] = atom_map[atom['n']]
                atom['id'] = atom_tp_map[atom['id']]

        size = len(self.bonds)
        for ibond in reversed(xrange(size)):
            bond = self.bonds[ibond]
            if bond['n'] not in bond_map:
                del self.bonds[ibond]
            else:
                bond['n'] = bond_map[bond['n']]
                bond['id'] = bond_tp_map[bond['id']]
                bond['i'] = atom_map[bond['i']]
                bond['j'] = atom_map[bond['j']]

        self.sections['Atoms'] = []
        for atom in self.atoms:
            self.sections['Atoms'].append(atomline(atom))
        self.sections['Masses'] = []
        for att in self.atomtypes:
            self.sections['Masses'].append(massline(att))
        self.sections['Bonds'] = []
        for bond in self.bonds:
            self.sections['Bonds'].append(bondline(bond))
        self.sections['Bond Coeffs'] = []
        for bdt in self.bondtypes:
            self.sections['Bond Coeffs'].append(bdtline(bdt))

# --------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description = 'Remove Drude dipole polarization from LAMMPS data file')
    parser.add_argument('infile', help = 'input LAMMPS data file')
    parser.add_argument('outfile', help = 'output LAMMPS data file')
    args = parser.parse_args()

    data = Data(args.infile)
    data.depolarize()
    data.write(args.outfile)

if __name__ == '__main__':
    main()
