#Version 0.1.1
#27/03/2025
#Author: Tom Wood

import sys
import re

class Structure:
    """Holds the information contained within the cif file in a dictionary format"""
    def __init__(self, props):
        self.props = props
    
    def __str__(self):
        s = ''
        cell_params = ["a", "b", "c", "alpha", "beta", "gamma"]
        if "_chemical_formula_structural" in self.props:
            s += ''.join(f"{self.props['_chemical_formula_structural']}")
            s += '\n'
        else:
            s += "Unknown phase\n"
        if "_space_group_name_H-M_alt" in self.props:
            s += f"Space group = {''.join(self.props['_space_group_name_H-M_alt'].split())}\n"
        elif "_symmetry_space_group_name_H-M" in self.props:
            s += f"Space group = {''.join(self.props['_symmetry_space_group_name_H-M'].split())}\n"
        if "_space_group_IT_number" in self.props:
            s += f"Space group number = {self.props['_space_group_IT_number']}\n"
        elif "_symmetry_Int_Tables_number" in self.props:
            s += f"Space group number = {self.props['_symmetry_Int_Tables_number']}\n"
        s += "Cell parameters:\n"
        for cp in cell_params:
            s += f"\t{cp} = "
            if f"_cell_length_{cp}" in self.props:
                cpstr = f"_cell_length_{cp}"
                s += f"{self.props[cpstr]}\n"
            elif f"_cell_angle_{cp}" in self.props:
                castr = f"_cell_angle_{cp}"
                s += f"{self.props[castr]}\n"
            else:
                s += "unknown\n"
        for k in self.props.keys():
            if 'loop' in k:
                if '_atom_site_label' in self.props[k]:
                    loop = self.props[k]
                    s += 'Atomic sites:\n'
                    for i in range(len(loop['_atom_site_label'])):
                        s += loop['_atom_site_label'][i]
                        for dim in ['x', 'y', 'z']:
                            s += f" {dim} "
                            s += loop[f'_atom_site_fract_{dim}'][i]
                        s += f" occ {loop['_atom_site_occupancy'][i]} "
                        if '_atom_site_U_iso_or_equiv' in loop:
                            s += f" beq {float(loop['_atom_site_U_iso_or_equiv'][i].split('(')[0]) * 8 * 3.14159**2:.3f}"
                        elif '_atom_site_B_iso_or_equiv' in loop:
                            s += f" beq {loop['_atom_site_B_iso_or_equiv'][i]}"
                        else:
                            s += ' beq .'
                        s += '\n'
        return s

    @classmethod
    def from_cif(cls, fname):
        """Read a cif file and return a Structure instance
        
        Args:
            fname (str): name of cif file
        """
        props = {}
        loop_props = False
        loop_items = False
        loop_count = -1
        lpcount = 0
        cont_prop = False
        multi_line = False
        with open(fname, 'r') as f:
            for line in f:
                #this bit is clunky, but couldn't find an elegant regex way of doing it.
                l = re.split('\'|\"', line)
                if len(l) > 1:
                    l = []
                    inquot = False
                    quot_start = None
                    s_start = 0
                    for i, char in enumerate(line):
                        if char == '\'' or char == '\"':
                            if inquot:
                                inquot = False
                                l.append(line[quot_start + 1:i])
                                s_start = i + 1
                            else:
                                inquot = True
                                quot_start = i
                        else:
                            if inquot is False:
                                if char == ' ' or i == len(line) - 1:
                                    l.append(line[s_start:i])
                                    s_start = i + 1
                else:
                    l = line.split()
                if not len(l):
                    continue
                if l[0][0] == '#':
                    continue
                if cont_prop:
                    if l[0] == ';':
                        if multi_line:
                            props.update({current_prop : new_prop})
                            cont_prop = False
                        else:
                            multi_line = True
                            new_prop = ''
                        continue
                    else:
                        if multi_line:
                            new_prop += line
                            continue
                        else:
                            props.update({current_prop : l[0]})
                            cont_prop = False
                if loop_props:
                    ln = f"loop{loop_count}"
                    if ln not in props:
                        props.update({ln : {}})
                        current_props = {}
                    if l[0][0] == '_':
                        props[ln].update({l[0] : []})
                        current_props.update({lpcount : l[0]})
                        lpcount += 1
                        continue
                    else:
                        loop_props = False
                        loop_items = True
                if loop_items:
                    if l[0] == 'loop_':
                        loop_items = False
                        loop_props = True
                        loop_count += 1
                        lpcount = 0
                        continue
                    elif l[0][0] == '_':
                        loop_items = False
                    else:
                        i = 0
                        for item in l:
                            if not item:
                                continue
                            props[f"loop{loop_count}"][current_props[i]].append(item)
                            i += 1
                        continue
                if l[0][0] == '_':
                    if len(l) == 1:
                        current_prop = l[0]
                        cont_prop = True
                    else:
                        props.update({l[0] : l[1]})
                elif l[0] == 'loop_':
                    loop_props = True
                    loop_count += 1
                    lpcount = 0
        return cls(props)
    
    def write_topas_output(self, fname_out):
        """Writes a TOPAS-style str into a separate file
        
        Args:
            fname_out (str): name of file to write
        """
        s = ''
        cell_params = ["a", "b", "c", "alpha", "beta", "gamma"]
        if "_chemical_formula_structural" in self.props:
            s += 'phase_name "'
            s += ''.join(f"{self.props['_chemical_formula_structural']}".split())
            s += '"\n'
        if "_space_group_name_H-M_alt" in self.props:
            s += f"space_group \"{''.join(self.props['_space_group_name_H-M_alt'].split())}\"\n"
        elif "_symmetry_space_group_name_H-M" in self.props:
            s += f"space_group \"{''.join(self.props['_symmetry_space_group_name_H-M'].split())}\"\n"
        else:
            if "_space_group_IT_number" in self.props:
                s += f"space_group {self.props['_space_group_IT_number']}\n"
            elif "_symmetry_Int_Tables_number" in self.props:
                s += f"space_group {self.props['_symmetry_Int_Tables_number']}\n"
        for cp in cell_params:
            s += f"{cp[:2]} "
            if f"_cell_length_{cp}" in self.props:
                cpstr = f"_cell_length_{cp}"
                s += f"{self.props[cpstr].split('(')[0]}\n"
            elif f"_cell_angle_{cp}" in self.props:
                castr = f"_cell_angle_{cp}"
                s += f"{self.props[castr].split('(')[0]}\n"
            else:
                s += ".\n"
        for k in self.props.keys():
            if 'loop' in k:
                if '_atom_site_label' in self.props[k]:
                    loop = self.props[k]
                    for i in range(len(loop['_atom_site_label'])):
                        site_lab = re.split('(\d+)', loop['_atom_site_label'][i])[0]
                        s += f"site {site_lab}"
                        for dim in ['x', 'y', 'z']:
                            s += f" {dim} "
                            s += loop[f'_atom_site_fract_{dim}'][i].split('(')[0]
                        s += f" occ {site_lab} {loop['_atom_site_occupancy'][i].split('(')[0]} "
                        if '_atom_site_U_iso_or_equiv' in loop:
                            s += f" beq {float(loop['_atom_site_U_iso_or_equiv'][i].split('(')[0]) * 8 * 3.14159**2:.3f}"
                        elif '_atom_site_B_iso_or_equiv' in loop:
                            s += f" beq {loop['_atom_site_B_iso_or_equiv'][i].split('(')[0]}"
                        else:
                            s += ' beq .'
                        s += '\n'
        with open(fname_out, 'w') as f:
            f.write(s)
        return

if __name__ == "__main__":
    fname_in = sys.argv[1]
    if len(sys.argv) < 3:
        fname_out = ''.join(fname_in.split('.')[:-1]) + '_ciftop.txt'
    else:
        fname_out = sys.argv[2]
    struc = Structure.from_cif(fname_in)
    struc.write_topas_output(fname_out)