# Drag sail lifetime validation

From the activated `orekit-env` environment, run `python validate_leo_drag.py` in the folder containing `leo_drag.py` and `orekit-data`. The validator first checks the 24-hour reference trajectory at 1 m position and 0.01 m/s velocity tolerances. It then runs the lifetime grid and saves the results to `sail_lifetime_grid.csv`.

The grid tests starting altitudes of 300, 350, and 400 km at fixed sail angles of 0°, 30°, 60°, and 90°, with a no-sail baseline at each altitude. Angles are measured between the sail normal and the atmosphere-relative flow, with 0° face-on and 90° edge-on. The model uses a 120 kg spacecraft, a 10 m² sail, a 0.8 m² bus drag area, and a drag coefficient of 2.2. Lifetime is measured from deployment to the first descending crossing of 120 km. Cases that do not cross within the 90-day window are reported as greater than 90 days.

The reference trajectory matched and all 22 trend checks passed. Face-on lifetimes were 1.953, 5.954, and 15.152 days at 300, 350, and 400 km, respectively. The edge-on cases matched their no-sail baselines because the model assigns the sail zero projected area at 90°. At 400 km, both cases exceeded the 90-day window. These results apply to the specified spacecraft parameters, fixed-angle model, and atmospheric conditions starting on 1 January 2024.
https://gitlab.orekit.org/orekit/orekit-data/-/archive/main/orekit-data-main.zip
