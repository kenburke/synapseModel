from clustering import io
import pytest
import os

# @pytest.mark.parametrize("filename,names,numbers", [
#     ("276.pdb", ["HIS", "HIS", "HIS", "HIS", "ASP"], [55, 57, 201, 230, 301]),
#     ("4629.pdb", ["ASP", "THR", "ARG", "SER", "LYS", "TYR", "SER", "ASN", "ASP"], [10, 14, 41, 118, 151, 157, 176, 177, 180]),
# ])
# def test_residues(filename, names, numbers):
#     filepath = os.path.join("data", filename)
# 
#     activesite = io.read_active_site(filepath)
# 
#     assert [residue.type for residue in activesite.residues] == names
#     assert [residue.number for residue in activesite.residues] == numbers
# 
# 
# @pytest.mark.parametrize("filename,residue_number,atoms,xs,ys,zs", [
#     ("276.pdb", 0,
#         ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
#         [42.050, 42.764, 42.377, 42.961, 44.309, 45.049, 45.926, 45.019, 46.397, 45.876],
#         [26.570, 27.774, 28.890, 29.025, 27.606, 28.823, 28.732, 30.140, 29.964, 30.848],
#         [22.707, 23.102, 22.153, 21.077, 23.162, 23.706, 24.777, 23.323, 25.028, 24.167]),
# ])
# def test_atoms(filename, residue_number, atoms, xs, ys, zs):
#     filepath = os.path.join("data", filename)
# 
#     activesite = io.read_active_site(filepath)
# 
#     residue = activesite.residues[residue_number]
# 
#     assert [atom.type for atom in residue.atoms] == atoms
#     assert [atom.coords for atom in residue.atoms] == list(zip(xs, ys, zs))
