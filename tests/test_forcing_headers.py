"""Check that weight headers describe actual forcing metadata without model assumptions."""

from pathlib import Path
import sys
import tempfile
import unittest

import netCDF4 as nc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import spatial_aggregation as spatial


class ForcingHeaderTests(unittest.TestCase):
    def test_input_title_is_used_without_inventing_a_product_or_version(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "arbitrary_filename.nc"
            with nc.Dataset(path, "w") as dataset:
                dataset.title = " Another forcing\n  dataset "
                dataset.source = "A secondary description"
            lines = spatial._gridweights_header(path, "2026-10-23")
            self.assertIn("#   Forcing dataset: Another forcing dataset", lines)
            text = "\n".join(lines)
            self.assertNotIn("Model info", text)
            self.assertNotIn("OHMS", text)
            self.assertNotIn("CaSR", text)
            self.assertNotIn("A secondary description", text)
            self.assertTrue(all("\n" not in line and line.startswith("#") for line in lines))

    def test_blank_title_falls_back_to_source_and_absent_metadata_is_omitted(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forcing.nc"
            with nc.Dataset(path, "w") as dataset:
                dataset.title = " \t\n "
                dataset.source = "Provided forcing product"
            lines = spatial._gridweights_header(path, "2026-10-23")
            self.assertIn("#   Forcing dataset: Provided forcing product", lines)
            with nc.Dataset(path, "a") as dataset:
                dataset.delncattr("source")
            lines = spatial._gridweights_header(path, "2026-10-23")
            self.assertFalse(any("Forcing dataset" in line for line in lines))
            self.assertTrue(any("Created by: climate-netcdf-aggregation" in line for line in lines))
            self.assertTrue(any("Repository: https://github.com/" in line for line in lines))


if __name__ == "__main__":
    unittest.main()
