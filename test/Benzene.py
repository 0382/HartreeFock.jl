# compare Benzene result with pyscf

from time import time_ns
import numpy as np
from pyscf import gto
from pyscf import scf

a0 = 0.529_177_210_903
hartree = 27.211386245988

def Benzene(dCC:float, dCH:float):
    sqrt3_half = np.sqrt(3)/2
    RH = dCC + dCH
    return gto.M(
        atom = [
            ['C', np.array([  1.,          0., 0.]) * dCC],
            ['C', np.array([ 0.5,  sqrt3_half, 0.]) * dCC],
            ['C', np.array([-0.5,  sqrt3_half, 0.]) * dCC],
            ['C', np.array([ -1.,          0., 0.]) * dCC],
            ['C', np.array([-0.5, -sqrt3_half, 0.]) * dCC],
            ['C', np.array([ 0.5, -sqrt3_half, 0.]) * dCC],
            ['H', np.array([  1.,          0., 0.]) * RH],
            ['H', np.array([ 0.5,  sqrt3_half, 0.]) * RH],
            ['H', np.array([-0.5,  sqrt3_half, 0.]) * RH],
            ['H', np.array([ -1.,          0., 0.]) * RH],
            ['H', np.array([-0.5, -sqrt3_half, 0.]) * RH],
            ['H', np.array([ 0.5, -sqrt3_half, 0.]) * RH],
        ],
        basis = "sto-3g",
        verbose = 0
    )

def Benzene_pes(dCC_range, dCH_range):
    NC = len(dCC_range)
    NH = len(dCH_range)
    Emat = np.zeros((NC, NH))
    for i in range(NC):
        for j in range(NH):
            ben = Benzene(dCC_range[i], dCH_range[j])
            hf = scf.RHF(ben)
            hf.chkfile = None
            Emat[i,j] = hf.kernel()
    pos = np.argmin(Emat)
    i, j = pos // NH, pos % NH
    print(f"minimal energy is {Emat[i,j]*hartree}eV, with dCC = {dCC_range[i]:.2f}Å, dCH = {dCH_range[j]:.2f}Å")

if __name__ == "__main__":
    start = time_ns()
    Benzene_pes(np.linspace(1.3, 1.5, 11), np.linspace(1, 1.2, 11))
    stop = time_ns()
    print("time used:", (stop-start)/1e9)
