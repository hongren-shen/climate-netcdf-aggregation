"""Check that metadata filenames retain the forcing's supported product identity."""

from pathlib import Path
import sys
import tempfile
import unittest

import netCDF4 as nc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from spatial_aggregation import _get_dataset_id_and_varstr


class DatasetIdentityTests(unittest.TestCase):
    def identity(self, variables, filename="forcing.nc", dataset_id=None,
                 title="Canadian Surface Reanalysis (CaSR)"):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / filename
            with nc.Dataset(path, "w") as dataset:
                dataset.title = title
                if dataset_id is not None:
                    dataset.dataset_id = dataset_id
                dataset.createDimension("time", 1)
                for name in variables:
                    dataset.createVariable(name, "f4", ("time",))[:] = [1]
            return _get_dataset_id_and_varstr(path)

    def test_bundled_inputs_keep_exact_dataset_and_variable_labels(self):
        for variable in ("A_PR0_SFC", "A_RAIN_SFC", "A_SNOW_SFC", "A_TT_1_5m"):
            with self.subTest(variable=variable):
                path = ROOT / "data/casr_3.1/SE" / ("CaSR_v3_1_" + variable + ".nc")
                self.assertEqual(_get_dataset_id_and_varstr(path), ("CaSR_v3_1", variable))

    def test_main_variable_version_takes_priority_over_filename(self):
        self.assertEqual(self.identity(["CaSR_v3_2_A_PR0_SFC"], "CaSR_v3_1_A_PR0_SFC.nc"),
                         ("CaSR_v3_2", "A_PR0_SFC"))

    def test_filename_supplies_version_when_main_variable_has_no_version(self):
        for title in ("Canadian Surface Reanalysis (CaSR)", ""):
            with self.subTest(title=title):
                self.assertEqual(self.identity(["A_PR0_SFC"], "CaSR_v4_0_A_PR0_SFC.nc",
                                               title=title),
                                 ("CaSR_v4_0", "A_PR0_SFC"))

    def test_generic_casr_metadata_does_not_invent_a_version(self):
        self.assertEqual(self.identity(["CaSR_A_PR0_SFC"]), ("CaSR", "A_PR0_SFC"))

    def test_explicit_dataset_id_is_preserved(self):
        self.assertEqual(self.identity(["Custom_A_PR0_SFC"], "CaSR_v3_1_A_PR0_SFC.nc",
                                       dataset_id="Custom"),
                         ("Custom", "A_PR0_SFC"))

    def test_conflicting_variable_versions_do_not_choose_a_filename_version(self):
        variables = ["CaSR_v3_1_A_PR0_SFC", "CaSR_v3_2_A_SNOW_SFC"]
        self.assertEqual(self.identity(variables, "CaSR_v3_1_forcing.nc"),
                         ("CaSR", "v3_1_A_PR0_SFC_v3_2_A_SNOW_SFC"))

    def test_non_casr_inference_is_unchanged(self):
        self.assertEqual(self.identity(["Example_PR0"], title="Another prepared product"),
                         ("Example", "PR0"))


if __name__ == "__main__":
    unittest.main()
