# compare H2O result with pyscf

from time import time_ns
import numpy as np
from pyscf import gto
from pyscf import scf

a0 = 0.529_177_210_903
hartree = 27.211386245988

def H2O(r:float, theta:float):
    theta = theta * np.pi / 180
    return gto.M(
        atom = [
            ['O', [0, 0, 0]],
            ['H', [r * np.cos(theta/2), r * np.sin(theta/2), 0]],
            ['H', [r * np.cos(theta/2), -r * np.sin(theta/2), 0]]],
        basis="sto-3g",
        verbose=0
    )

def H2O_pes(r_range, theta_range):
    Nr = len(r_range)
    Nt = len(theta_range)
    Emat = np.zeros((Nr, Nt))
    for i, r in enumerate(r_range):
        for j, theta in enumerate(theta_range):
            hf = scf.RHF(H2O(r, theta))
            hf.chkfile = None
            Emat[i,j] = hf.kernel()
    pos = np.argmin(Emat)
    i, j = pos // Nt, pos % Nt
    print(f"minimal energy is {Emat[i,j]*hartree}eV, with the bond length = {r_range[i]:.2f}Å, bond angle = {theta_range[j]:.1f}ᵒ")

if __name__ == "__main__":
    start = time_ns()
    H2O_pes(np.linspace(0.9, 1.1, 21), np.linspace(85, 115, 31))
    stop = time_ns()
    print("time used:", (stop-start)/1e9)
