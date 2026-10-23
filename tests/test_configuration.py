"""Exercise driver configuration and spawn propagation without production jobs."""

import copy
import importlib
import io
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflow"
sys.path.insert(0, str(WORKFLOW))
driver = importlib.import_module("aggregate_driver")


def worker_snapshot(keys):
    """Inspect the imported driver globals inside a newly spawned process."""
    imported_driver = importlib.import_module("aggregate_driver")
    return {key: getattr(imported_driver, key) for key in keys}


class DriverCliTests(unittest.TestCase):
    def argument(self, destination, *, negative=False):
        parser = driver.build_argument_parser()
        action = next((item for item in parser._actions if item.dest == destination), None)
        self.assertIsNotNone(action, f"No CLI argument exposes {destination}")
        candidates = [option for option in action.option_strings
                      if option.startswith("--no-") == negative and option.startswith("--")]
        self.assertTrue(candidates, f"No {'negative' if negative else 'positive'} switch for {destination}")
        return candidates[0]

    def test_empty_arguments_preserve_editable_top_level_defaults(self):
        config = driver.parse_runtime_config([])
        self.assertIsInstance(config, dict)
        self.assertTrue(config)
        for key, value in config.items():
            with self.subTest(parameter=key):
                self.assertTrue(hasattr(driver, key), key)
                self.assertEqual(value, getattr(driver, key))
        self.assertEqual(config["hru_versions"], ["v2-0-1a", "v2-0-1b"])
        self.assertFalse(config["want_temporal_aggregation"])
        self.assertFalse(config["derive_daily_temp_minmax"])

    def test_list_path_and_performance_overrides_do_not_change_module_defaults(self):
        before = driver.parse_runtime_config([])
        arguments = [self.argument("region_for_process"), "SE", "SW",
                     self.argument("extract_gauge_list"), "02KA015", "02KF012",
                     self.argument("want_aggregated_targets"), "basin", "hru",
                     self.argument("hru_versions"), "v2-0-1b", "v2-0-1a",
                     self.argument("nc_files_process"), "custom_precip.nc", "custom_temp.nc",
                     self.argument("nc_data_folder"), "/tmp/isolated-forcing-input",
                     self.argument("output_folder"), "/tmp/isolated-forcing-output",
                     self.argument("num_cores"), "3",
                     self.argument("aggregate_compress_level"), "4",
                     self.argument("time_block_hours"), "7",
                     self.argument("output_time_chunk_hours"), "12",
                     self.argument("want_temporal_aggregation"),
                     self.argument("temporal_agg_freq"), "6H"]
        config = driver.parse_runtime_config(arguments)
        self.assertEqual(config["region_for_process"], ["SE", "SW"])
        self.assertEqual(config["extract_gauge_list"], ["02KA015", "02KF012"])
        self.assertEqual(config["want_aggregated_targets"], ["basin", "hru"])
        self.assertEqual(config["hru_versions"], ["v2-0-1b", "v2-0-1a"])
        self.assertEqual(config["nc_files_process"], ["custom_precip.nc", "custom_temp.nc"])
        self.assertEqual(config["nc_data_folder"], ["/tmp/isolated-forcing-input"])
        self.assertEqual(config["output_folder"], ["/tmp/isolated-forcing-output"])
        self.assertEqual(config["num_cores"], 3)
        self.assertEqual(config["aggregate_compress_level"], 4)
        self.assertEqual(config["time_block_hours"], 7)
        self.assertEqual(config["output_time_chunk_hours"], 12)
        self.assertTrue(config["want_temporal_aggregation"])
        self.assertEqual(config["temporal_agg_freq"], "6H")
        self.assertEqual(driver.parse_runtime_config([]), before)

    def test_negative_boolean_switches_reach_false(self):
        defaults = driver.parse_runtime_config([])
        public_parameters = {action.dest for action in driver.build_argument_parser()._actions}
        boolean_parameters = [key for key, value in defaults.items()
                              if isinstance(value, bool) and key in public_parameters]
        self.assertGreaterEqual(len(boolean_parameters), 8)
        config = driver.parse_runtime_config([self.argument(key, negative=True) for key in boolean_parameters])
        for key in boolean_parameters:
            with self.subTest(parameter=key):
                self.assertIs(config[key], False)

    def test_invalid_values_fail_during_parsing(self):
        invalid = [
            [self.argument("num_cores"), "0"],
            [self.argument("num_cores"), "-1"],
            [self.argument("aggregate_compress_level"), "10"],
            [self.argument("time_block_hours"), "0"],
            [self.argument("output_time_chunk_hours"), "0"],
            [self.argument("temporal_anchor_hour"), "24"],
            [self.argument("temporal_agg_freq"), "not-a-frequency"],
            [self.argument("want_aggregated_targets"), "hru",
             self.argument("hru_versions"), "v9-9-9a"],
            [self.argument("want_aggregated_targets"), "hru",
             self.argument("hru_versions"), "not-an-hru-version"],
        ]
        for arguments in invalid:
            with self.subTest(arguments=arguments):
                diagnostic = io.StringIO()
                with redirect_stderr(diagnostic), self.assertRaises(SystemExit) as caught:
                    driver.parse_runtime_config(arguments)
                self.assertEqual(caught.exception.code, 2)
                self.assertIn("error:", diagnostic.getvalue().lower())

    def test_owl_hydro_codes_and_legacy_names_select_the_same_layout(self):
        for code, alias, suffix in [("v2-0-1a", "semi-distributed", "a"),
                                    ("v2-0-1b", "routing-only", "b")]:
            with self.subTest(version=code):
                current = driver.normalize_hru_version(code, olrrp_version="v2-0-1")
                legacy = driver.normalize_hru_version(alias, olrrp_version="v2-0-1")
                self.assertEqual(current, legacy)
                self.assertEqual(current[0], suffix)
        with self.assertRaises(ValueError):
            driver.normalize_hru_version("v2-0-1a", olrrp_version="v2-0-2")

    def test_empty_gauge_selection_and_hourly_aliases(self):
        config = driver.parse_runtime_config(["--gauges", "--frequency", "6h"])
        self.assertEqual(config["extract_gauge_list"], [])
        self.assertEqual(config["temporal_agg_freq"], "6H")

    def test_runtime_configs_do_not_share_mutable_default_lists(self):
        first = driver.parse_runtime_config([])
        second = copy.deepcopy(driver.parse_runtime_config([]))
        first["hru_versions"].append("temporary-test-value")
        first["extract_gauge_list"].clear()
        self.assertEqual(driver.parse_runtime_config([]), second)

    def run_driver(self, arguments):
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PROJ_NETWORK="OFF")
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([sys.executable, str(WORKFLOW / "aggregate_driver.py"), *arguments],
                                    cwd=directory, env=environment, capture_output=True, text=True,
                                    timeout=30)
            self.assertEqual(list(Path(directory).iterdir()), [])
        return result

    def test_help_entrypoint_exits_before_data_io(self):
        result = self.run_driver(["--help"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("usage:", result.stdout.lower())
        self.assertIn(self.argument("hru_versions"), result.stdout)
        self.assertIn(self.argument("want_temporal_aggregation", negative=True), result.stdout)
        self.assertNotIn("Running with configuration", result.stdout)

    def test_invalid_entrypoint_values_fail_before_data_io(self):
        cases = [([self.argument("num_cores"), "0"], self.argument("num_cores")),
                 ([self.argument("temporal_agg_freq"), "not-a-frequency"], "freq"),
                 ([self.argument("want_aggregated_targets"), "hru",
                   self.argument("hru_versions"), "v9-9-9a"], "version")]
        for arguments, expected_message in cases:
            with self.subTest(arguments=arguments):
                result = self.run_driver(arguments)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn("error:", result.stderr.lower())
                self.assertIn(expected_message, result.stderr.lower())
                self.assertNotIn("FileNotFoundError", result.stderr)
                self.assertNotIn("Running with configuration", result.stdout)

    def test_spawn_initializer_propagates_actual_cli_overrides(self):
        config = driver.parse_runtime_config([
            self.argument("num_cores"), "1",
            self.argument("region_for_process"), "SW",
            self.argument("hru_versions"), "v2-0-1b",
            self.argument("want_aggregated_targets"), "hru",
            self.argument("want_gridded_forcing"),
            self.argument("want_temporal_aggregation"),
            self.argument("temporal_agg_freq"), "6H",
            self.argument("aggregate_overwrite", negative=True),
            self.argument("if_zip_gauge_folder", negative=True)])
        keys = ("num_cores", "region_for_process", "hru_versions", "want_aggregated_targets",
                "want_gridded_forcing", "want_temporal_aggregation", "temporal_agg_freq",
                "aggregate_overwrite", "if_zip_gauge_folder")
        expected = {key: config[key] for key in keys}
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=driver._configure_worker, initargs=(config,)) as pool:
            actual = pool.submit(worker_snapshot, keys).result(timeout=30)
        self.assertEqual(actual, expected)
        self.assertTrue(driver.aggregate_overwrite)
        self.assertTrue(driver.if_zip_gauge_folder)


class InputReadinessTests(unittest.TestCase):
    def test_default_inventory_is_bundled_se_24h_and_paths_are_absolute(self):
        config = driver.parse_runtime_config([])
        expected = ROOT / "data/gauge_inventory/SE/24H/gauge_info.csv"
        self.assertEqual(Path(config["gauge_info_path"][0]), expected)
        self.assertTrue(expected.is_file())
        for key in ("nc_data_folder", "gauge_info_path", "gauge_routing_product_folder", "output_folder"):
            self.assertTrue(Path(config[key][0]).is_absolute(), key)
        self.assertTrue(Path(config["grid_weights_generator"]).is_absolute())
        self.assertFalse(config["want_temporal_aggregation"])

    def test_bundled_inventory_selects_examples_without_forcing_daily_aggregation(self):
        with tempfile.TemporaryDirectory() as directory:
            forcing = Path(directory) / "SE/precip.nc"
            forcing.parent.mkdir()
            forcing.touch()
            config = driver.parse_runtime_config(["--input-root", directory, "--variables", "precip.nc"])
            table = driver.load_gauge_table_and_check_forcing(config)
            examples = table.loc[table["Obs_NM"].isin(config["extract_gauge_list"])]
            self.assertEqual(set(examples["Obs_NM"]), {"02KA015"})
            self.assertEqual(set(examples["Region"]), {"SE"})
            self.assertFalse(config["want_temporal_aggregation"])

    def test_inventory_errors_identify_missing_file_columns_and_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "gauges.csv"
            config = driver.parse_runtime_config(["--gauge-table", str(table)])
            with self.assertRaisesRegex(FileNotFoundError, "Gauge table not found"):
                driver.load_gauge_table_and_check_forcing(config)
            table.write_text("Region,Station\nSE,02KA015\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing columns: Obs_NM"):
                driver.load_gauge_table_and_check_forcing(config)
            table.write_text("Region,Obs_NM\nNW,02KA015\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "No gauges match"):
                driver.load_gauge_table_and_check_forcing(config)

    def test_missing_forcing_exits_cleanly_and_preserves_existing_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "outputs/SE"
            output.mkdir(parents=True)
            sentinel = output / "previous-result.txt"
            sentinel.write_text("keep this result", encoding="utf-8")
            environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PROJ_NETWORK="OFF")
            result = subprocess.run(
                [sys.executable, str(WORKFLOW / "aggregate_driver.py"),
                 "--input-root", str(root / "missing-forcing"),
                 "--variables", "precip.nc", "--output-root", str(root / "outputs")],
                cwd=root, env=environment, capture_output=True, text=True, timeout=30,
            )
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertIn("Prepared forcing NetCDF files not found", result.stderr)
            self.assertIn(str(root / "missing-forcing/SE/precip.nc"), result.stderr)
            self.assertNotIn("Gauge table not found", result.stderr)
            self.assertNotIn("Traceback", result.stderr)
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep this result")


class DriverFailureTests(unittest.TestCase):
    def run_workflow(self, directory, arguments, *, extra_environment=None):
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PROJ_NETWORK="OFF")
        environment["PATH"] = str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
        if extra_environment:
            environment.update(extra_environment)
        return subprocess.run(
            [sys.executable, str(WORKFLOW / "aggregate_driver.py"),
             "--variables", "CaSR_v3_1_A_PR0_SFC.nc", "--workers", "1", *arguments],
            cwd=directory, env=environment, capture_output=True, text=True, timeout=60,
        )

    def test_missing_routing_or_helper_preserves_existing_outputs(self):
        cases = [("--routing-root", "missing-routing", "Required routing files not found"),
                 ("--gridweights-generator", "missing-helper.py", "Grid-weight helper not found")]
        for option, filename, message in cases:
            with self.subTest(option=option), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                output = root / "outputs/SE"
                output.mkdir(parents=True)
                sentinel = output / "previous-result.txt"
                sentinel.write_text("retain previous output", encoding="utf-8")
                result = self.run_workflow(root, [option, str(root / filename),
                                                 "--output-root", str(output.parent)])
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertIn(message, result.stderr)
                self.assertEqual(sentinel.read_text(), "retain previous output")
                self.assertFalse((output / "02KA015.zip").exists())

    def test_missing_shapefile_companion_is_caught_before_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            routing = root / "routing/SE/02KA015"
            routing.mkdir(parents=True)
            source = ROOT / "data/Routing_product_simplified_extraction/SE/02KA015"
            for stem in ("finalcat_info_v2-0-1", "final_hru_v2-0-1a", "final_hru_v2-0-1b"):
                for suffix in (".shp", ".shx", ".dbf", ".prj"):
                    shutil.copyfile(source / (stem + suffix), routing / (stem + suffix))
            missing = routing / "final_hru_v2-0-1b.prj"
            missing.unlink()
            output = root / "outputs/SE"
            output.mkdir(parents=True)
            sentinel = output / "previous-result.txt"
            sentinel.write_text("retain previous output", encoding="utf-8")
            result = self.run_workflow(root, ["--routing-root", str(root / "routing"),
                                             "--output-root", str(output.parent)])
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertIn(str(missing), result.stderr)
            self.assertEqual(sentinel.read_text(), "retain previous output")

    def test_only_requested_or_available_hru_layouts_are_required(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            routing = root / "routing/SE/02KA015"
            routing.mkdir(parents=True)
            for stem in ("final_hru_v2-0-1a", "finalcat_info_v2-0-1", "outline"):
                for suffix in (".shp", ".shx", ".dbf", ".prj"):
                    (routing / (stem + suffix)).touch()
            arguments = ["--routing-root", str(root / "routing"), "--targets", "hru",
                         "--hru-versions", "v2-0-1a"]
            for gridded in (False, True):
                with self.subTest(gridded=gridded):
                    config = driver.parse_runtime_config(arguments + (["--gridded"] if gridded else []))
                    self.assertFalse(driver.load_gauge_table_and_check_forcing(config).empty)

    def test_helper_inside_regional_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            diagnostic = io.StringIO()
            with redirect_stderr(diagnostic), self.assertRaises(SystemExit) as caught:
                driver.parse_runtime_config(["--output-root", str(root),
                                             "--gridweights-generator", str(root / "SE/helper.py")])
            self.assertEqual(caught.exception.code, 2)
            self.assertIn("grid-weight helper", diagnostic.getvalue())

    def test_skipped_aggregation_is_an_explicit_worker_failure(self):
        for target in ("basin", "subbasin", "hru"):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                with mock.patch.multiple(driver, want_gridded_forcing=False,
                                         want_aggregated_targets=[target],
                                         hru_versions=["v2-0-1a"]), \
                        mock.patch.object(driver, "aggregate_netcdf_forcing", return_value=None):
                    with self.assertRaisesRegex(RuntimeError, "Requested output was not produced"):
                        driver._process_one_gauge("02KA015", root / "input.nc", root / "routing",
                                                  root / "outputs", "precip.nc", {})

    def test_runtime_worker_failure_exits_nonzero_without_packaging(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            helper = root / "failing_helper.py"
            helper.write_text("raise RuntimeError('intentional helper failure')\n", encoding="utf-8")
            output = root / "outputs"
            result = self.run_workflow(root, ["--gridweights-generator", str(helper),
                                             "--output-root", str(output)])
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Gauge processing failed", result.stderr)
            self.assertTrue((output / "SE/02KA015").is_dir())
            self.assertFalse((output / "SE/02KA015.zip").exists())

    def test_zip_failure_exits_nonzero_and_retains_generated_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            hook = root / "test-hook"
            hook.mkdir()
            (hook / "sitecustomize.py").write_text(
                "import shutil\n"
                "def fail_zip(*args, **kwargs):\n"
                "    raise OSError('intentional ZIP failure')\n"
                "shutil.make_archive = fail_zip\n", encoding="utf-8")
            python_path = str(hook) + os.pathsep + os.environ.get("PYTHONPATH", "")
            output = root / "outputs"
            result = self.run_workflow(root, ["--output-root", str(output)],
                                       extra_environment={"PYTHONPATH": python_path})
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("ZIP packaging failed", result.stderr)
            self.assertTrue(list((output / "SE/02KA015").glob("*.nc")))
            self.assertFalse((output / "SE/02KA015.zip").exists())


if __name__ == "__main__":
    unittest.main()
