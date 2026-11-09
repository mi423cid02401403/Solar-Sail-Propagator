from math import radians
from pathlib import Path

import orekit


def main():
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
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import Constants, IERSConventions
    data_dir = Path(__file__).resolve().parent / "orekit-data"
    if not data_dir.is_dir():
        raise SystemExit(f"Extract the Orekit data folder here:\n{data_dir}")

    manager = DataContext.getDefault().getDataProvidersManager()
    manager.addProvider(DirectoryCrawler(File(str(data_dir))))

    utc = TimeScalesFactory.getUTC()
    epoch = AbsoluteDate(2024, 1, 1, 0, 0, 0.0, utc)
    inertial = FramesFactory.getGCRF()
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)

    earth_radius_m = Constants.WGS84_EARTH_EQUATORIAL_RADIUS
    earth = OneAxisEllipsoid(
        earth_radius_m, Constants.WGS84_EARTH_FLATTENING, itrf
    )
    gravity = GravityFieldFactory.getNormalizedProvider(8, 8)

    # Sample near-circular orbit, roughly 400 km above Earth.
    initial_orbit = KeplerianOrbit(
        earth_radius_m + 400_000.0,
        0.001,          # eccentricity
        radians(51.6),  # inclination
        0.0,            # argument of perigee
        0.0,            # right ascension of the ascending node
        0.0,            # true anomaly
        PositionAngleType.TRUE,
        inertial,
        epoch,
        gravity.getMu(),
    )

    # Sample spacecraft parameters.
    mass_kg = 120.0
    drag_area_m2 = 0.8
    cd = 2.2

    weather = CssiSpaceWeatherData(
        CssiSpaceWeatherData.DEFAULT_SUPPORTED_NAMES)
    atmosphere = NRLMSISE00(weather, CelestialBodyFactory.getSun(), earth)

    orbit_type = OrbitType.CARTESIAN
    integrator = DormandPrince853IntegratorBuilder(
        0.001, 300.0, 1.0
    ).buildIntegrator(initial_orbit, orbit_type, PositionAngleType.TRUE)

    propagator = NumericalPropagator(integrator)
    propagator.setOrbitType(orbit_type)
    propagator.setInitialState(
        SpacecraftState(initial_orbit).withMass(mass_kg))
    propagator.addForceModel(HolmesFeatherstoneAttractionModel(itrf, gravity))
    propagator.addForceModel(
        DragForce(atmosphere, IsotropicDrag(drag_area_m2, cd)))

    final_state = propagator.propagate(epoch.shiftedBy(86_400.0))
    pv = final_state.getPVCoordinates()

    print("Final epoch (UTC):", final_state.getDate().toString(utc))
    print("GCRF position (m):", pv.getPosition())
    print("GCRF velocity (m/s):", pv.getVelocity())


if __name__ == "__main__":
    main()
