"""P02 - WT production MD extension (continuation of the 3 x 10 ns WT replicates).
segment 1: start from the last production frame (positions+box from cpptraj rst7),
           velocities re-drawn at 300 K (production wrote no velocities/state).
segment k>1: continue exactly from the saved OpenMM state of segment k-1.
Same force-field/integrator/barostat settings as production.py.
usage: python extend_wt.py <segment> [smoke]
"""
import sys, time
from pathlib import Path
from openmm import app, unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator

P02 = Path("/mnt/HHD/portfolio-projects/P02-MMGBSA")
segment = int(sys.argv[1])
smoke = len(sys.argv) > 2 and sys.argv[2] == "smoke"
SEEDS = [1] if smoke else [1, 2, 3]
NS = 10
STEPS = 5000 if smoke else int(NS * 1000 / 0.002)
REPORT_STEPS = 10000          # 20 ps, as in production.py
prefix = "smoke_" if smoke else ""

prmtop = app.AmberPrmtopFile(str(P02 / "system_WT.prmtop"))
system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds, rigidWater=True)
system.addForce(MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))
platform = Platform.getPlatformByName("HIP")

t_start = time.time()
for seed in SEEDS:
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    integrator.setRandomNumberSeed(seed * 100 + segment)
    sim = app.Simulation(prmtop.topology, system, integrator, platform)
    if segment == 1:
        crd = app.AmberInpcrdFile(str(P02 / f"restart_WT_seed{seed}_prod_last.rst7"))
        sim.context.setPositions(crd.positions)
        sim.context.setPeriodicBoxVectors(*crd.boxVectors)
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed * 100 + segment)
    else:
        sim.loadState(str(P02 / f"ext{segment-1}_WT_seed{seed}.state.xml"))
    stem = f"{prefix}production_WT_seed{seed}_ext{segment}"
    sim.reporters.append(app.DCDReporter(str(P02 / f"{stem}.dcd"), REPORT_STEPS))
    sim.reporters.append(app.StateDataReporter(str(P02 / f"{stem}.log"), REPORT_STEPS * 5, step=True, time=True,
                         potentialEnergy=True, temperature=True, density=True, speed=True))
    t0 = time.time()
    sim.step(STEPS)
    dt = time.time() - t0
    ns = STEPS * 0.002 / 1000
    print(f"WT seed{seed} ext{segment}: {ns:.3f} ns in {dt/3600:.3f} h ({ns/(dt/86400):.1f} ns/day), "
          f"total {(time.time()-t_start)/3600:.2f} h", flush=True)
    if not smoke:
        sim.saveState(str(P02 / f"ext{segment}_WT_seed{seed}.state.xml"))
        st = sim.context.getState(getPositions=True)
        with open(P02 / f"{stem}_final.pdb", "w") as f:
            app.PDBFile.writeFile(prmtop.topology, st.getPositions(), f)
print(f"WT ext{segment}: DONE {'(smoke)' if smoke else ''}", flush=True)
