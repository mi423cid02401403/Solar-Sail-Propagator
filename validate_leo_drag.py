"""Validate the reference trajectory and run the drag-sail lifetime grid."""

import csv
from itertools import product
import math
from pathlib import Path
import subprocess
import sys

# Replace these example values with your study inputs.
START_ALTITUDES_KM = (300, 350, 400)
SAIL_ANGLES_DEG = (0, 30, 60, 90)  # Sail normal relative to the incoming flow
MASS_KG = 120.0
BUS_AREA_M2 = 0.8
SAIL_AREA_M2 = 10.0
DRAG_COEFFICIENT = 2.2
INCLINATION_DEG = 51.6
START_DATE = (2024, 1, 1)
STOP_ALTITUDE_KM = 120.0
MAX_DAYS = 90.0
SECONDS_PER_DAY = 86400.0


EXPECTED_EPOCH = "2024-01-02T00:00:00.000"
EXPECTED_POSITION_M = (-5877813.925036346, -1766341.318500905, -2863833.1260059983)
EXPECTED_VELOCITY_M_S = (3777.9679448135, -4358.8249280127, -5059.0830795985)
POSITION_TOLERANCE_M = 1.0
VELOCITY_TOLERANCE_M_S = 0.01


def read_value(output, label):
    for line in output.splitlines():
        if line.startswith(label):
            return line[len(label):].strip()
    raise ValueError(f"Missing output line: {label}")


def read_vector(output, label):
    value = read_value(output, label)
    components = tuple(
        float(part.replace(",", "").strip())
        for part in value.strip("{}").split(";")
    )
    if len(components) != 3 or not all(map(math.isfinite, components)):
        raise ValueError(f"Invalid vector for {label} {value}")
    return components


def validate_reference():
    script = Path(__file__).resolve().with_name("leo_drag.py")
    if not script.is_file():
        raise SystemExit(f"FAIL: Put validate_leo_drag.py beside leo_drag.py. Missing: {script}")

    result = subprocess.run(
        [sys.executable, str(script)],
        cwd=script.parent,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        print(result.stderr, end="", file=sys.stderr)
        raise SystemExit("FAIL: leo_drag.py did not complete.")

    try:
        epoch = read_value(result.stdout, "Final epoch (UTC):")
        position = read_vector(result.stdout, "GCRF position (m):")
        velocity = read_vector(result.stdout, "GCRF velocity (m/s):")
    except ValueError as error:
        raise SystemExit(f"FAIL: {error}") from None

    position_error = math.dist(position, EXPECTED_POSITION_M)
    velocity_error = math.dist(velocity, EXPECTED_VELOCITY_M_S)
    print(f"Final epoch (UTC): {epoch}")
    print(f"Position difference: {position_error:.6f} m")
    print(f"Velocity difference: {velocity_error:.9f} m/s")

    if epoch != EXPECTED_EPOCH:
        raise SystemExit(f"FAIL: expected final epoch {EXPECTED_EPOCH}")
    if position_error > POSITION_TOLERANCE_M:
        raise SystemExit(f"FAIL: position differs by more than {POSITION_TOLERANCE_M:g} m.")
    if velocity_error > VELOCITY_TOLERANCE_M_S:
        raise SystemExit(f"FAIL: velocity differs by more than {VELOCITY_TOLERANCE_M_S:g} m/s.")
    print("PASS: propagation matches the reference run.")


def projected_area(angle_deg):
    if angle_deg is None or angle_deg == 90:
        return BUS_AREA_M2
    return BUS_AREA_M2 + SAIL_AREA_M2 * math.cos(math.radians(angle_deg))


def run_lifetime_grid():
    import orekit

    orekit.initVM()

    from java.io import File
    from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
    from org.orekit.data import DataContext, DirectoryCrawler
    from org.orekit.forces.drag import DragForce, IsotropicDrag
    from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
    from org.orekit.forces.gravity.potential import GravityFieldFactory
    from org.orekit.frames import FramesFactory
    from org.orekit.models.earth.atmosphere import NRLMSISE00
    from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
    from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
    from org.orekit.propagation.events import AltitudeDetector
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import Constants, IERSConventions

    if any(h <= STOP_ALTITUDE_KM for h in START_ALTITUDES_KM):
        raise SystemExit("Starting altitudes must be above the stop altitude.")
    if any(not 0 <= angle <= 90 for angle in SAIL_ANGLES_DEG):
        raise SystemExit("Sail angles must be between 0 and 90 degrees.")
    if min(MASS_KG, BUS_AREA_M2, SAIL_AREA_M2, DRAG_COEFFICIENT, MAX_DAYS) <= 0:
        raise SystemExit("Mass, areas, drag coefficient and duration must be positive.")

    project_dir = Path(__file__).resolve().parent
    data_dir = project_dir / "orekit-data"
    if not data_dir.is_dir():
        raise SystemExit(f"Missing Orekit data folder: {data_dir}")

    manager = DataContext.getDefault().getDataProvidersManager()
    manager.addProvider(DirectoryCrawler(File(str(data_dir))))
    utc = TimeScalesFactory.getUTC()
    epoch = AbsoluteDate(*START_DATE, 0, 0, 0.0, utc)
    end = epoch.shiftedBy(MAX_DAYS * SECONDS_PER_DAY)
    weather = CssiSpaceWeatherData(CssiSpaceWeatherData.DEFAULT_SUPPORTED_NAMES)
    if epoch.compareTo(weather.getMinDate()) < 0 or end.compareTo(weather.getMaxDate()) >= 0:
        raise SystemExit("CSSI data do not cover the requested window. Adjust START_DATE or MAX_DAYS.")

    gcrf = FramesFactory.getGCRF()
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
    radius = Constants.WGS84_EARTH_EQUATORIAL_RADIUS
    earth = OneAxisEllipsoid(radius, Constants.WGS84_EARTH_FLATTENING, itrf)
    gravity = GravityFieldFactory.getNormalizedProvider(8, 8)
    atmosphere = NRLMSISE00(weather, CelestialBodyFactory.getSun(), earth)
    orbit_type = OrbitType.CARTESIAN
    angle_type = PositionAngleType.TRUE
    builder = DormandPrince853IntegratorBuilder(0.001, 300.0, 1.0)
    stop_altitude_m = STOP_ALTITUDE_KM * 1000.0

    output = project_dir / "sail_lifetime_grid.csv"
    results = {}
    with output.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
            "initial_altitude_km", "case", "angle_deg", "effective_area_m2",
            "status", "lifetime_days", "simulated_days", "final_altitude_km",
        ])
        for altitude_km, angle in product(START_ALTITUDES_KM, (None, *SAIL_ANGLES_DEG)):
            area = projected_area(angle)
            case = "no_sail" if angle is None else "deployed_sail"
            label = f"{altitude_km:g} km, " + ("no sail" if angle is None else f"{angle:g} deg")
            print(f"Running {label} ...", flush=True)

            # Circular initial orbit, starting on the equator.
            orbit = KeplerianOrbit(
                radius + altitude_km * 1000.0, 0.0, math.radians(INCLINATION_DEG),
                0.0, 0.0, 0.0, angle_type, gcrf, epoch, gravity.getMu(),
            )
            propagator = NumericalPropagator(builder.buildIntegrator(orbit, orbit_type, angle_type))
            propagator.setOrbitType(orbit_type)
            propagator.setInitialState(SpacecraftState(orbit).withMass(MASS_KG))
            propagator.addForceModel(HolmesFeatherstoneAttractionModel(itrf, gravity))
            propagator.addForceModel(DragForce(atmosphere, IsotropicDrag(area, DRAG_COEFFICIENT)))
            propagator.addEventDetector(AltitudeDetector(60.0, 0.001, stop_altitude_m, earth))

            final_state = propagator.propagate(end)
            final_date = final_state.getDate()
            days = final_date.durationFrom(epoch) / SECONDS_PER_DAY
            geodetic = earth.transform(final_state.getPosition(), gcrf, final_date)
            crossed = final_date.compareTo(end) < 0 or geodetic.getAltitude() <= stop_altitude_m
            status = "threshold_reached" if crossed else "beyond_horizon"
            results[altitude_km, angle] = (days, crossed)
            writer.writerow([
                altitude_km, case, angle, area, status,
                f"{days:.6f}" if crossed else "", f"{days:.6f}",
                f"{geodetic.getAltitude() / 1000.0:.3f}",
            ])
            file.flush()
            print(f"  Lifetime: {days:.3f} days" if crossed else f"  Lifetime > {MAX_DAYS:g} days")

    print(f"Results saved to {output}")
    return results


def compare_lifetimes(earlier, later):
    """Check the expected ordering with one second of numerical slack."""
    earlier_days, earlier_crossed = earlier
    later_days, later_crossed = later
    tolerance_days = 1.0 / SECONDS_PER_DAY

    if earlier_crossed and earlier_days <= later_days + tolerance_days:
        return "PASS"
    if later_crossed and earlier_days > later_days + tolerance_days:
        return "WARN"
    return "INCONCLUSIVE"


def validate_lifetime_trends(results):
    print("\nLifetime trend checks")
    for altitude in sorted(START_ALTITUDES_KM):
        angles = sorted(SAIL_ANGLES_DEG)
        for first, second in zip(angles, angles[1:]):
            status = compare_lifetimes(results[altitude, first], results[altitude, second])
            print(f"{status}: {altitude:g} km, {first:g} deg vs {second:g} deg")
        if angles:
            status = compare_lifetimes(results[altitude, angles[0]], results[altitude, None])
            print(f"{status}: {altitude:g} km, {angles[0]:g} deg vs no sail")

    altitudes = sorted(START_ALTITUDES_KM)
    for angle in (None, *sorted(SAIL_ANGLES_DEG)):
        label = "no sail" if angle is None else f"{angle:g} deg"
        for lower, higher in zip(altitudes, altitudes[1:]):
            status = compare_lifetimes(results[lower, angle], results[higher, angle])
            print(f"{status}: {label}, {lower:g} km vs {higher:g} km")


if __name__ == "__main__":
    validate_reference()
    validate_lifetime_trends(run_lifetime_grid())
