"""Protect the example/library boundary from sibling-script dependencies."""
import ast
import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / 'application/Chern-number-calculation.py'


def test_example_uses_only_standard_or_installed_libraries():
    source = EXAMPLE.read_text(encoding='utf-8')
    tree = ast.parse(source)
    allowed = set(sys.stdlib_module_names) | {'numpy', 'matplotlib', 'BerryPy', 'pygbz2d'}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert all(n.name.split('.')[0] in allowed for n in node.names)
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0 and node.module.split('.')[0] in allowed
    assert 'sys.path' not in source and 'playground' not in source


def test_experimental_import_does_not_load_plotting_or_demos():
    # A fresh process prevents modules cached by other tests from hiding an
    # accidental optional dependency or a plotting-backend side effect.
    code = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith(('matplotlib', 'BerryPy', 'playground', 'demo_', 'chern_on_gbz')):
            raise AssertionError(fullname)
sys.meta_path.insert(0, Block())
sys.path.insert(0, sys.argv[1])
from pygbz2d.experimental import build_band_mesh, refine_mesh, integrate_chern
'''
    subprocess.run([sys.executable, '-c', code, str(ROOT/'src')], check=True,
                   capture_output=True, text=True, timeout=30)


def test_copied_example_can_integrate_and_write_flux(tmp_path):
    copied = tmp_path / EXAMPLE.name
    copied.write_bytes(EXAMPLE.read_bytes())
    spec = importlib.util.spec_from_file_location('copied_chern_example', copied)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    assert example.DATA == tmp_path/'data/QWZ-benchmark'
    from pygbz2d.experimental import regular_torus_mesh
    verts, tri = regular_torus_mesh(13)
    H = lambda b1, b2: example.hamiltonian(b1, b2, 1.2, .4)
    energies = np.array([min(np.linalg.eigvals(H(*np.exp(1j*v))), key=lambda z: z.real)
                         for v in verts])
    path = tmp_path/'mesh.npz'
    np.savez(path, verts=verts, triangles=tri, E=energies, mu1=np.zeros(len(verts)), mu2=np.zeros(len(verts)))
    result = example.calculate_chern(path, H, mode='biorthogonal', save_flux=True)
    assert abs(result.chern+1) < 1e-12
    with np.load(tmp_path/'mesh_chern_flux.npz') as saved:
        np.testing.assert_array_equal(saved['flux'], saved['flux_complex'].real)
        assert saved['link_rule'].item() == 'biorthogonal-reciprocal-v1'
    assert all('playground' not in row['path'] for row in example.source_provenance())
