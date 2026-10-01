"""P02 - NVT -> NPT equilibration, then a short production benchmark to measure
ns/day on this host's HIP GPU (informs the final production length decision).
"""
import pickle
import sys
import time
from pathlib import Path

from openmm import app, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
genotype = sys.argv[1]

prmtop = app.AmberPrmtopFile(str(P02 / f"system_{genotype}.prmtop"))
with open(P02 / f"minimized_{genotype}.state.pkl", "rb") as f:
    state = pickle.load(f)

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True)
platform = Platform.getPlatformByName("HIP")

# --- NVT: 50 ps, heavy-atom restraint on protein+ligand, heat 0->300 K ---
integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
from openmm import CustomExternalForce
restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
restraint.addGlobalParameter("k", 50.0 * unit.kilocalorie_per_mole / unit.angstrom**2)
restraint.addPerParticleParameter("x0"); restraint.addPerParticleParameter("y0"); restraint.addPerParticleParameter("z0")
positions0 = state["positions"]
for atom in prmtop.topology.atoms():
    if atom.residue.name not in ("WAT", "Na+", "Cl-", "HOH", "NA", "CL"):
        if atom.element is not None and atom.element.symbol != "H":
            restraint.addParticle(atom.index, positions0[atom.index])
restraint_force_idx = system.addForce(restraint)

simulation = app.Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(positions0)
simulation.context.setPeriodicBoxVectors(*state["box"])
simulation.context.setVelocitiesToTemperature(50 * unit.kelvin)

t0 = time.time()
n_heat_steps = 25000  # 50 ps
for i in range(10):
    target_T = 50 + (300 - 50) * (i + 1) / 10
    integrator.setTemperature(target_T * unit.kelvin)
    simulation.step(n_heat_steps // 10)
print(f"{genotype} NVT heating done ({time.time()-t0:.0f}s)")

# --- NPT: 200 ps, barostat on, restraint relaxed to 5 kcal/mol/A^2 ---
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))
simulation.context.reinitialize(preserveState=True)
simulation.context.setParameter("k", 5.0)
t1 = time.time()
simulation.step(100000)  # 200 ps
print(f"{genotype} NPT (restrained, relaxed) done ({time.time()-t1:.0f}s)")

# --- remove restraint entirely, short free equilibration ---
system.removeForce(restraint_force_idx)
integrator2 = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
simulation2 = app.Simulation(prmtop.topology, system, integrator2, platform)
st = simulation.context.getState(getPositions=True, getVelocities=True)
simulation2.context.setPositions(st.getPositions())
simulation2.context.setVelocities(st.getVelocities())
simulation2.context.setPeriodicBoxVectors(*st.getPeriodicBoxVectors())
t2 = time.time()
simulation2.step(25000)  # 50 ps free equil
print(f"{genotype} free equilibration (unrestrained, 50 ps) done ({time.time()-t2:.0f}s)")

# save equilibrated state for production
eq_state = simulation2.context.getState(getPositions=True, getVelocities=True)
with open(P02 / f"equilibrated_{genotype}.state.pkl", "wb") as f:
    pickle.dump({
        "positions": eq_state.getPositions(asNumpy=True),
        "velocities": eq_state.getVelocities(asNumpy=True),
        "box": eq_state.getPeriodicBoxVectors(asNumpy=True),
    }, f)

# --- benchmark: 100 ps timed production ---
t3 = time.time()
n_bench_steps = 50000  # 100 ps
simulation2.step(n_bench_steps)
dt = time.time() - t3
ns_per_day = (n_bench_steps * 0.002 / 1000) / (dt / 86400)
print(f"{genotype} BENCHMARK: {n_bench_steps} steps (100 ps) in {dt:.1f}s -> {ns_per_day:.1f} ns/day")

with open(P02 / f"benchmark_{genotype}.txt", "w") as f:
    f.write(f"{ns_per_day:.2f}\n")
