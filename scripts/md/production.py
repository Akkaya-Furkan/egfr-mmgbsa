"""P02 - production MD: 3 x 10 ns replicates per system, from the equilibrated
state, re-randomizing velocities per replicate (standard practice for
independent replicates from one equilibrated structure).
"""
import pickle
import sys
import time
from pathlib import Path

from openmm import app, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
genotype = sys.argv[1]
SEEDS = [1, 2, 3]
NS_PER_REPLICATE = 10
STEPS_PER_REPLICATE = int(NS_PER_REPLICATE * 1000 / 0.002)  # 2 fs timestep
REPORT_INTERVAL_PS = 20
REPORT_INTERVAL_STEPS = int(REPORT_INTERVAL_PS / 0.002)

prmtop = app.AmberPrmtopFile(str(P02 / f"system_{genotype}.prmtop"))
with open(P02 / f"equilibrated_{genotype}.state.pkl", "rb") as f:
    eq = pickle.load(f)

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                              constraints=app.HBonds, rigidWater=True)
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))
platform = Platform.getPlatformByName("HIP")

t_start = time.time()
for seed in SEEDS:
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integrator.setRandomNumberSeed(seed)
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(eq["positions"])
    simulation.context.setPeriodicBoxVectors(*eq["box"])
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

    dcd_path = P02 / f"production_{genotype}_seed{seed}.dcd"
    simulation.reporters.append(app.DCDReporter(str(dcd_path), REPORT_INTERVAL_STEPS))
    log_path = P02 / f"production_{genotype}_seed{seed}.log"
    simulation.reporters.append(app.StateDataReporter(
        str(log_path), REPORT_INTERVAL_STEPS * 5, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, speed=True))

    t0 = time.time()
    simulation.step(STEPS_PER_REPLICATE)
    dt = time.time() - t0
    ns_per_day = (NS_PER_REPLICATE) / (dt / 86400)
    print(f"{genotype} seed{seed}: {NS_PER_REPLICATE} ns done in {dt/3600:.2f} h "
          f"({ns_per_day:.1f} ns/day), total elapsed {(time.time()-t_start)/3600:.2f} h",
          flush=True)

    # save final frame coordinates too (for a quick sanity-check structure)
    final_state = simulation.context.getState(getPositions=True)
    with open(P02 / f"production_{genotype}_seed{seed}_final.pdb", "w") as f:
        app.PDBFile.writeFile(prmtop.topology, final_state.getPositions(), f)

print(f"{genotype}: ALL 3 REPLICATES DONE, total {(time.time()-t_start)/3600:.2f} h", flush=True)
