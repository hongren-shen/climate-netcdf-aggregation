"""Check version provenance when citation metadata changes or is unavailable."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import spatial_aggregation as spatial


class SoftwareVersionTests(unittest.TestCase):
    def test_version_updates_without_reloading_the_module(self):
        with tempfile.TemporaryDirectory() as directory:
            citation = Path(directory) / "CITATION.cff"
            with patch.object(spatial, "_CITATION_PATH", citation):
                for field, expected in [
                    ("1.0.0", "v1.0.0"),
                    ('"1.2.0" # release', "v1.2.0"),
                    ("'v2.0.0-rc.1+build.2'", "v2.0.0-rc.1+build.2"),
                ]:
                    with self.subTest(field=field):
                        citation.write_text(
                            f"cff-version: 1.2.0\nversion: {field}\n"
                            "references:\n  - type: software\n    version: 9.9.9\n",
                            encoding="utf-8",
                        )
                        self.assertEqual(spatial._software_creator(),
                                         f"climate-netcdf-aggregation {expected}")

    def test_copied_script_without_citation_still_has_a_clean_creator(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "CITATION.cff"
            with patch.object(spatial, "_CITATION_PATH", missing):
                self.assertEqual(spatial._software_creator(), "climate-netcdf-aggregation")

    def test_unusable_metadata_does_not_add_an_incorrect_version(self):
        with tempfile.TemporaryDirectory() as directory:
            citation = Path(directory) / "CITATION.cff"
            with patch.object(spatial, "_CITATION_PATH", citation):
                for contents in [
                    "cff-version: 1.2.0\nreferences:\n  - version: 9.9.9\n",
                    "version: \n",
                    "version: null\n",
                    "version: [1, 2, 3]\n",
                    'version: "1.0.0\n',
                    "version: 1.0.0\nversion: 2.0.0\n",
                ]:
                    with self.subTest(contents=contents):
                        citation.write_text(contents, encoding="utf-8")
                        self.assertEqual(spatial._software_creator(), "climate-netcdf-aggregation")
                citation.write_bytes(b"\xff\xfe\xfd")
                self.assertEqual(spatial._software_creator(), "climate-netcdf-aggregation")


if __name__ == "__main__":
    unittest.main()
