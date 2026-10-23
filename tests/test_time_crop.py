"""Check inclusive CF time selection and lossless copying of packed forcing."""

from pathlib import Path
import sys
import tempfile
import unittest

import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from spatial_aggregation import crop_netcdf_time_range


class TimeCropTests(unittest.TestCase):
    def test_bundled_gregorian_times_and_bounds_are_cropped_inclusively(self):
        source = ROOT / "data/casr_3.1/SE/CaSR_v3_1_A_PR0_SFC.nc"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "subset.nc"
            crop_netcdf_time_range(source, target, (2020, 1, 1, 1), "2020-01-01 03:00", copy_chunk_hours=2)
            with nc.Dataset(source) as original, nc.Dataset(target) as subset:
                self.assertEqual(len(subset.dimensions["time"]), 3)
                for name, variable in original.variables.items():
                    copied = subset[name]
                    self.assertEqual(copied.dimensions, variable.dimensions)
                    self.assertEqual(copied.dtype, variable.dtype)
                    for key in variable.ncattrs():
                        np.testing.assert_equal(copied.getncattr(key), variable.getncattr(key))
                    variable.set_auto_maskandscale(False)
                    copied.set_auto_maskandscale(False)
                    selection = tuple(slice(0, 3) if dim == "time" else slice(None) for dim in variable.dimensions)
                    np.testing.assert_array_equal(copied[...], variable[selection] if selection else variable[...])

    def test_packed_values_fill_codes_scalar_and_nonleading_time_axis_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target = Path(directory) / "packed.nc", Path(directory) / "subset.nc"
            raw = np.array([[10, -32767, 30, 40], [50, 60, 70, 80]], dtype="i2")
            with nc.Dataset(source, "w") as dataset:
                dataset.createDimension("cell", 2)
                dataset.createDimension("time", 4)
                time = dataset.createVariable("time", "f8", ("time",))
                time.units = "hours since 2001-01-01 00:00:00"
                time.calendar = "noleap"
                time[:] = [1, 2, 3, 4]
                value = dataset.createVariable("forcing", "i2", ("cell", "time"), fill_value=-32767)
                value.scale_factor = 0.1
                value.add_offset = 250.0
                value.units = "K"
                value.set_auto_maskandscale(False)
                value[:] = raw
                dataset.createVariable("crs", "i4").assignValue(42)
            crop_netcdf_time_range(source, target, "2001-01-01 02:00", "2001-01-01 03:00", copy_chunk_hours=1)
            with nc.Dataset(target) as subset:
                self.assertEqual(subset["time"].calendar, "noleap")
                np.testing.assert_array_equal(subset["time"][:], [2, 3])
                self.assertEqual(subset["forcing"]._FillValue, -32767)
                self.assertEqual(subset["forcing"].scale_factor, 0.1)
                self.assertEqual(subset["forcing"].add_offset, 250.0)
                self.assertTrue(np.ma.getmaskarray(subset["forcing"][:])[0, 0])
                self.assertAlmostEqual(subset["forcing"][1, 1], 257.0)
                subset.set_auto_maskandscale(False)
                np.testing.assert_array_equal(subset["forcing"][:], raw[:, 1:3])
                self.assertEqual(subset["crs"].getValue(), 42)

    def test_invalid_ranges_do_not_create_or_overwrite_outputs(self):
        source = ROOT / "data/casr_3.1/SE/CaSR_v3_1_A_PR0_SFC.nc"
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "untouched.nc"
            target.write_bytes(b"existing output")
            for start, end in [("2019-01-01", "2019-01-02"), ("2020-01-03", "2020-01-01")]:
                with self.subTest(start=start), self.assertRaises(ValueError):
                    crop_netcdf_time_range(source, target, start, end)
                self.assertEqual(target.read_bytes(), b"existing output")
            with self.assertRaisesRegex(ValueError, "different files"):
                crop_netcdf_time_range(source, source, "2020-01-01", "2020-01-02")


if __name__ == "__main__":
    unittest.main()
