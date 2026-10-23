"""Check the bundled forcing's spatial and interval-coverage contracts."""

from datetime import datetime
from pathlib import Path
import unittest

import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "data/casr_3.1/SE"
VARIABLES = [("CaSR_v3_1_A_PR0_SFC", 1, "kg m-2"),
             ("CaSR_v3_1_A_RAIN_SFC", 1, "kg m-2"),
             ("CaSR_v3_1_A_SNOW_SFC", 1, "kg m-2"),
             ("CaSR_v3_1_A_TT_1_5m", 3, "degC")]


class BundledForcingTests(unittest.TestCase):
    def test_four_variables_cover_the_complete_2020_utc_year_on_the_same_grid(self):
        reference_coordinates = None
        for name, cadence, units in VARIABLES:
            with self.subTest(variable=name), nc.Dataset(INPUTS / (name + ".nc")) as dataset:
                self.assertEqual(len(dataset.dimensions["rlat"]), 13)
                self.assertEqual(len(dataset.dimensions["rlon"]), 13)
                self.assertEqual(dataset.example_gauge_id, "02KA015")
                self.assertEqual(dataset.example_buffer_grid_cells, 5)
                time = dataset["time"]
                bounds = dataset[time.bounds][:]
                begin, end = nc.date2num([datetime(2020, 1, 1), datetime(2021, 1, 1)],
                                        time.units, time.calendar)
                self.assertEqual(len(time), 8784 // cadence)
                self.assertEqual(bounds[0, 0], begin)
                self.assertEqual(bounds[-1, 1], end)
                np.testing.assert_array_equal(bounds[:, 1], time[:])
                np.testing.assert_array_equal(bounds[1:, 0], bounds[:-1, 1])
                np.testing.assert_array_equal(bounds[:, 1] - bounds[:, 0], cadence)
                variable = dataset[name]
                self.assertEqual(variable.units, units)
                self.assertEqual(variable.dimensions, ("time", "rlat", "rlon"))
                self.assertEqual(variable.dtype, np.dtype("float32"))
                self.assertFalse(np.ma.getmaskarray(variable[:]).any())
                self.assertTrue(np.isfinite(variable[:]).all())
                self.assertEqual(dataset[variable.grid_mapping].grid_mapping_name,
                                 "rotated_latitude_longitude")
                coordinates = {key: dataset[key][:] for key in ("lat", "lon", "rlat", "rlon")}
                if reference_coordinates is None:
                    reference_coordinates = coordinates
                else:
                    for key in coordinates:
                        np.testing.assert_array_equal(coordinates[key], reference_coordinates[key])


if __name__ == "__main__":
    unittest.main()
