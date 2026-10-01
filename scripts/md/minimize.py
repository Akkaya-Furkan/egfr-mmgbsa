"""P02 - energy minimization (OpenMM, HIP platform) for each system.

Two-stage: (1) restrain protein+ligand heavy atoms (relax water/ions/H's first),
(2) unrestrained minimization to convergence.
"""
import sys
import time
from pathlib import Path

from openmm import app, unit, Platform, CustomExternalForce
from openmm import LangevinMiddleIntegrator

P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
genotype = sys.argv[1]

prmtop = app.AmberPrmtopFile(str(P02 / f"system_{genotype}.prmtop"))
inpcrd = app.AmberInpcrdFile(str(P02 / f"system_{genotype}.inpcrd"))

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True)

# Stage 1: restrain protein+ligand heavy atoms (k=100 kcal/mol/A^2) so water/ions relax first
restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
restraint.addGlobalParameter("k", 100.0 * unit.kilocalorie_per_mole / unit.angstrom**2)
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")
positions = inpcrd.positions
for atom in prmtop.topology.atoms():
    if atom.residue.name not in ("WAT", "Na+", "Cl-", "HOH", "NA", "CL"):
        if atom.element is not None and atom.element.symbol != "H":
            restraint.addParticle(atom.index, positions[atom.index].value_in_unit(unit.nanometer))
system.addForce(restraint)

integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
platform = Platform.getPlatformByName("HIP")
simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

t0 = time.time()
e0 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"{genotype} stage1 initial energy: {e0}")
simulation.minimizeEnergy(maxIterations=2000)
e1 = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"{genotype} stage1 (restrained) final energy: {e1}  ({time.time()-t0:.0f}s)")

# Stage 2: remove restraint, full minimization
system.removeForce(system.getNumForces() - 1)
integrator2 = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation2 = app.Simulation(prmtop.topology, system, integrator2, platform)
state = simulation.context.getState(getPositions=True, getVelocities=False)
simulation2.context.setPositions(state.getPositions())
simulation2.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())

t1 = time.time()
simulation2.minimizeEnergy(maxIterations=5000)
e2 = simulation2.context.getState(getEnergy=True).getPotentialEnergy()
print(f"{genotype} stage2 (unrestrained) final energy: {e2}  ({time.time()-t1:.0f}s)")

final_state = simulation2.context.getState(getPositions=True, getVelocities=False)
with open(P02 / f"minimized_{genotype}.pdb", "w") as f:
    app.PDBFile.writeFile(prmtop.topology, final_state.getPositions(), f)
import pickle
with open(P02 / f"minimized_{genotype}.state.pkl", "wb") as f:
    pickle.dump({
        "positions": final_state.getPositions(asNumpy=True),
        "box": final_state.getPeriodicBoxVectors(asNumpy=True),
    }, f)
print(f"{genotype}: DONE, total {time.time()-t0:.0f}s")
