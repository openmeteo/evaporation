import datetime as dt
import os
import shutil
import sys
import tempfile
import textwrap
from io import StringIO
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import click
import numpy as np
import pandas as pd
from click.testing import CliRunner
from htimeseries import HTimeseries
from osgeo import gdal, osr

from evaporation import cli

from .test_evaporation import senegal_tzinfo

gdal.UseExceptions()


def create_geotiff_file(filename, value):
    geo_transform = (-16.25, 1.0, 0, 16.217, 0, 1.0)
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    f = gdal.GetDriverByName("GTiff").Create(filename, 2, 1, 1, gdal.GDT_Float32)
    try:
        f.SetGeoTransform(geo_transform)
        f.SetProjection(wgs84.ExportToWkt())
        f.GetRasterBand(1).WriteArray(value)
    finally:
        # Close the dataset
        f = None


class NonExistentConfigFileTestCase(TestCase):
    def setUp(self):
        runner = CliRunner(mix_stderr=False)
        self.result = runner.invoke(cli.main, ["nonexistent.conf"])

    def test_exit_status(self):
        self.assertEqual(self.result.exit_code, 1)

    def test_error_message(self):
        self.assertIn(
            "No such file or directory: 'nonexistent.conf'", self.result.stderr
        )


class WrongPointConfigurationTestCase(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.configfilename = os.path.join(self.tempdir, "evaporation.conf")
        htsfilename = os.path.join(self.tempdir, "wind_speed.hts")
        Path(htsfilename).touch()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_nonexistent_log_level_raises_error(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                    base_dir = {self.tempdir}
                    albedo = 0.23
                    nighttime_solar_radiation_ratio = 0.8
                    elevation = 8
                    step_length = 60
                    unit_converter_pressure = x / 10.0
                    unit_converter_solar_radiation = x * 3600 /1e6
                    loglevel = NONEXISTENT_LOG_LEVEL
                    """.format(
                        self=self
                    )
                )
            )
        msg = "loglevel must be one of ERROR, WARNING, INFO, DEBUG"
        with self.assertRaisesRegex(click.ClickException, msg):
            cli.App(self.configfilename).run()

    def test_missing_step_raises_error(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                    base_dir = {self.tempdir}
                    albedo = 0.23
                    nighttime_solar_radiation_ratio = 0.8
                    elevation = 8
                    """
                ).format(self=self)
            )
        with self.assertRaisesRegex(click.ClickException, "step_length"):
            cli.App(self.configfilename).run()

    def test_missing_albedo_raises_error(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                    base_dir = {self.tempdir}
                    nighttime_solar_radiation_ratio = 0.8
                    elevation = 8
                    step_length = 60
                    """
                ).format(self=self)
            )
        with self.assertRaisesRegex(click.ClickException, "albedo"):
            cli.App(self.configfilename).run()


class WrongSpatialConfigurationTestCase(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.configfilename = os.path.join(self.tempdir, "evaporation.conf")
        htsfilename = os.path.join(self.tempdir, "wind_speed.tif")
        Path(htsfilename).touch()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_missing_elevation_raises_error(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                step_length = 60
                """
                ).format(self=self)
            )
        with self.assertRaisesRegex(click.ClickException, "elevation"):
            cli.App(self.configfilename).run()

    def test_single_albedo_with_wrong_domain_float_inputs(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 2
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        with self.assertRaisesRegex(ValueError, "Albedo must be between 0.0 and 1.0"):
            cli.App(self.configfilename).run()

    def test_seasonal_albedo_configuration_with_not_enough_arguments(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = a01.tiff a02.tiff a11.tiff a12.tiff
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        msg = "Albedo must be either one item or 12 space-separated items"
        with self.assertRaisesRegex(ValueError, msg):
            cli.App(self.configfilename).run()

    def test_seasonal_albedo_with_wrong_domain_mixin_inputs(self):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 1 2 0.34 0.24 0.45 0.4
                        0.34 0.12 a00.tif 0.78 0.78 2
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        with self.assertRaisesRegex(ValueError, "Albedo must be between 0.0 and 1.0"):
            cli.App(self.configfilename).run()


class CorrectPointConfigurationTestCase(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.configfilename = os.path.join(self.tempdir, "evaporation.conf")
        htsfilename = os.path.join(self.tempdir, "wind_speed.hts")
        Path(htsfilename).touch()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @patch("evaporation.cli.ProcessAtPoint")
    def test_executes(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                    base_dir = {self.tempdir}
                    albedo = 0.23
                    nighttime_solar_radiation_ratio = 0.8
                    step_length = 60
                    unit_converter_pressure = x / 10.0
                    unit_converter_solar_radiation = x * 3600 /1e6
                    """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessAtPoint")
    def test_albedo_configuration_as_one_number(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()


class CorrectSpatialConfigurationTestCase(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.configfilename = os.path.join(self.tempdir, "evaporation.conf")
        htsfilename = os.path.join(self.tempdir, "wind_speed.tif")
        Path(htsfilename).touch()
        self._create_albedo_files()

    def _create_albedo_files(self):
        items = (
            "00",
            "01",
            "02",
            "03",
            "04",
            "05",
            "06",
            "07",
            "08",
            "09",
            "10",
            "11",
            "12",
        )
        for item in items:
            filename = os.path.join(self.tempdir, "a{}.tif".format(item))
            create_geotiff_file(filename, np.array([[0.23, 0.44]]))

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    @patch("evaporation.cli.ProcessSpatial")
    def test_albedo_configuration_as_one_grid(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = a00.tif
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_seasonal_albedo_configuration_as_12_grids(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = a01.tif a02.tif a03.tif a04.tif a05.tif a06.tif
                         a07.tif a08.tif a09.tif a10.tif a11.tif a12.tif
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_seasonal_albedo_configuration_as_mix_numbers_and_grids(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 0.23 a02.tif a03.tif a04.tif a05.tif a06.tif
                         a07.tif a08.tif a09.tif a10.tif a11.tif a12.tif
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_run_app_seasonal_albedo_with_float_sample_inputs(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 0.10 0.23 0.34 0.24 0.45 0.46
                         0.34 0.12 0.14 0.78 0.78 0.12
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_run_app_with_seasonal_albedo_with_grid_sample_inputs(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = a00.tif a00.tif a00.tif a00.tif a00.tif a00.tif
                         a00.tif a00.tif a00.tif a00.tif a00.tif a00.tif
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_run_app_with_seasonal_albedo_with_mix_sample_inputs(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = a00.tif 0.23 a00.tif a00.tif a00.tif a00.tif
                         a00.tif a00.tif 0.23 a00.tif 0.23 a00.tif
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()

    @patch("evaporation.cli.ProcessSpatial")
    def test_seasonal_albedo_configuration_as_12_numbers(self, m):
        with open(self.configfilename, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                step_length = 60
                albedo = 0.10 0.23 0.34 0.24 0.45 0.46
                         0.34 0.12 0.14 0.78 0.78 0.12
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                """
                ).format(self=self)
            )
        cli.App(self.configfilename).run()
        m.return_value.execute.assert_called_once_with()


class CorrectConfigurationWithLogFileTestCase(TestCase):
    @patch("evaporation.cli.ProcessAtPoint")
    def test_creates_log_file(self, *args):
        with tempfile.TemporaryDirectory() as dirname:
            configfilename = os.path.join(dirname, "evaporation.conf")
            logfilename = os.path.join(dirname, "vaporize.log")
            with open(configfilename, "w") as f:
                f.write(
                    textwrap.dedent(
                        """\
                        base_dir = {}
                        albedo = 0.23
                        nighttime_solar_radiation_ratio = 0.8
                        step_length = 60
                        unit_converter_pressure = x / 10.0
                        unit_converter_solar_radiation = x * 3600 /1e6
                        logfile = {}
                        """
                    ).format(dirname, logfilename)
                )
            htsfilename = os.path.join(dirname, "somefile.hts")
            Path(htsfilename).touch()
            cli.App(configfilename).run()
            self.assertTrue(os.path.exists(logfilename))


class HtsTestCase(TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.config_file = os.path.join(self.tempdir, "vaporize.conf")
        self.savedcwd = os.getcwd()

    def tearDown(self):
        os.chdir(self.savedcwd)
        shutil.rmtree(self.tempdir)

    def setup_input_file(self, step, basename, value, missing=None):
        filename = os.path.join(self.tempdir, basename + ".hts")
        timestamp = step == "hourly" and "2014-10-01 15:00" or "2014-07-06"
        with open(filename, "w") as f:
            f.write("Title={}\n".format(basename.replace("_", " ").capitalize()))
            f.write("Timezone=CVT (UTC-0100)\n")
            if missing != "location":
                f.write("Location=-16.25 16.217 4326\n")
            if missing != "altitude":
                f.write("Altitude={}\n".format(step == "hourly" and "8" or "100"))
            f.write("\n{},{},\n".format(timestamp, value))

    def setup_hourly_input_files(self, missing=None):
        self.setup_input_file("hourly", "temperature", 38, missing=missing)
        self.setup_input_file("hourly", "humidity", 52, missing=missing)
        self.setup_input_file("hourly", "wind_speed", 3.3, missing=missing)
        self.setup_input_file("hourly", "pressure", 1013, missing=missing)
        self.setup_input_file("hourly", "solar_radiation", 681, missing=missing)

    def setup_daily_input_files(self):
        self.setup_input_file("daily", "temperature_max", 21.5)
        self.setup_input_file("daily", "temperature_min", 12.3)
        self.setup_input_file("daily", "humidity_max", 84)
        self.setup_input_file("daily", "humidity_min", 63)
        self.setup_input_file("daily", "wind_speed", 2.078)
        self.setup_input_file("daily", "sunshine_duration", 9.25)

    def setup_config_file(self, step_length):
        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                    base_dir = {self.tempdir}
                    albedo = 0.23
                    step_length = {step_length}
                    """
                ).format(self=self, step_length=step_length)
            )
            if step_length == 60:
                f.write(
                    textwrap.dedent(
                        """\
                        nighttime_solar_radiation_ratio = 0.8
                        unit_converter_pressure = x / 10.0
                        unit_converter_solar_radiation = x * 3600 / 1e6
                        """
                    )
                )

    def test_hourly(self):
        self.setup_hourly_input_files()
        self.setup_config_file(60)

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(self.tempdir, "evaporation.hts")
        assert not os.path.exists(result_filename)

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file and that the file is correct
        with open(result_filename) as f:
            t = HTimeseries(f)
        expected_result = pd.DataFrame(
            data={"value": [0.63], "flags": [""]},
            columns=("value", "flags"),
            index=[dt.datetime(2014, 10, 1, 15, 0)],
        )
        expected_result.index.name = "date"
        pd.testing.assert_frame_equal(t.data, expected_result, check_less_precise=2)

    def test_daily(self):
        self.setup_daily_input_files()
        self.setup_config_file(1440)

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(self.tempdir, "evaporation.hts")
        assert not os.path.exists(result_filename)

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file and that the file is correct
        with open(result_filename) as f:
            t = HTimeseries(f)
        expected_result = pd.DataFrame(
            data={"value": [3.9], "flags": [""]},
            columns=["value", "flags"],
            index=[dt.datetime(2014, 7, 6)],
        )
        expected_result.index.name = "date"
        pd.testing.assert_frame_equal(t.data, expected_result, check_less_precise=1)

    def test_missing_location(self):
        self.setup_hourly_input_files(missing="location")
        self.setup_config_file(60)
        msg = (
            "Incorrect or unspecified or inconsistent locations in the time series "
            "files."
        )
        with self.assertRaisesRegex(click.ClickException, msg):
            cli.App(self.config_file).run()

    def test_missing_altitude(self):
        self.setup_hourly_input_files(missing="altitude")
        self.setup_config_file(60)
        msg = (
            "Incorrect or unspecified or inconsistent altitudes in the time series "
            "files."
        )
        with self.assertRaisesRegex(click.ClickException, msg):
            cli.App(self.config_file).run()


class SpatialTestCase(TestCase):
    def setup_input_file(self, variable, value, with_date=True):
        """
        Saves value, which is an np array, to a GeoTIFF file whose name is
        based on variable.
        """
        if not with_date:
            filename = variable
        elif not isinstance(self.timestamp, dt.datetime):
            # Daily stuff
            filename = variable + "-" + self.timestamp.strftime("%Y-%m-%d")
        else:
            # Hourly stuff
            filename = variable + "-" + self.timestamp.strftime("%Y-%m-%d-%H-%M%z")

        filename = os.path.join(self.tempdir, filename + ".tif")
        f = gdal.GetDriverByName("GTiff").Create(filename, 2, 1, 1, gdal.GDT_Float32)

        nodata = -1.70141183e38

        value[np.isnan(value)] = nodata
        try:
            f.SetMetadataItem("TIMESTAMP", self.timestamp.isoformat())
            f.SetGeoTransform(self.geo_transform)
            f.SetProjection(self.wgs84.ExportToWkt())
            f.GetRasterBand(1).SetNoDataValue(nodata)
            f.GetRasterBand(1).WriteArray(value)
        finally:
            # Close the dataset
            f = None

    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.config_file = os.path.join(self.tempdir, "vaporize.conf")
        self.savedcwd = os.getcwd()

        # Prepare data common to all input files
        self.geo_transform = (-16.25, 1.0, 0, 16.217, 0, 1.0)
        self.wgs84 = osr.SpatialReference()
        self.wgs84.ImportFromEPSG(4326)

        # Prepare sample albedo grid
        filename = os.path.join(self.tempdir, "a00.tif")
        create_geotiff_file(filename, np.array([[0.23, 0.44]]))

        # Save standard error (some tests change it)
        self.saved_stderr = sys.stderr

    def tearDown(self):
        os.chdir(self.savedcwd)
        shutil.rmtree(self.tempdir)

    def test_execute_notz(self):
        # Prepare input files without time zone
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0)
        self.setup_input_file("temperature-notz", np.array([[38.0, 28.0]]))
        self.setup_input_file("humidity-notz", np.array([[52.0, 42.0]]))
        self.setup_input_file("wind_speed-notz", np.array([[3.3, 2.3]]))
        self.setup_input_file("pressure-notz", np.array([[1013.0, 1013.0]]))
        self.setup_input_file("solar_radiation-notz", np.array([[681.0, 403.0]]))

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                step_length = 60
                unit_converter_pressure = x / 10.0
                unit_converter_solar_radiation = x * 3600 / 1e6
                temperature_prefix = temperature-notz
                humidity_prefix = humidity-notz
                wind_speed_prefix = wind_speed-notz
                solar_radiation_prefix = solar_radiation-notz
                pressure_prefix = pressure-notz
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d-%H-%M%z")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Execute
        with self.assertRaisesRegex(Exception, "time zone"):
            cli.App(self.config_file).run()

        # Verify the output file still doesn't exist
        self.assertFalse(os.path.exists(result_filename))

    def test_execute_daily(self):
        # Prepare input files
        self.timestamp = dt.date(2014, 7, 6)
        self.setup_input_file("temperature_max", np.array([[21.5, 28]]))
        self.setup_input_file("temperature_min", np.array([[12.3, 15]]))
        self.setup_input_file("humidity_max", np.array([[84.0, 70.0]]))
        self.setup_input_file("humidity_min", np.array([[63.0, 60.0]]))
        self.setup_input_file("wind_speed", np.array([[2.078, 2.244]]))
        self.setup_input_file("sunshine_duration", np.array([[9.25, 9.0]]))

        # Also setup an output file that has no corresponding input files
        rogue_output_file = os.path.join(self.tempdir, "evaporation-2013-01-01.tif")
        with open(rogue_output_file, "w") as f:
            f.write("irrelevant contents")

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                elevation = 100
                step_length = 1440
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Verify the rogue output file is still here
        self.assertTrue(os.path.exists(rogue_output_file))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the rogue output file is gone
        self.assertFalse(os.path.exists(rogue_output_file))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-07-06")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(), np.array([[3.9, 4.8]]), decimal=1
        )
        fp = None

    def test_execute_daily_with_radiation(self):
        """Same as test_execute_daily, except that we use solar radiation
           instead of sunshine duration."""
        # Prepare input files
        self.timestamp = dt.date(2014, 7, 6)
        self.setup_input_file("temperature_max", np.array([[21.5, 28]]))
        self.setup_input_file("temperature_min", np.array([[12.3, 15]]))
        self.setup_input_file("humidity_max", np.array([[84.0, 70.0]]))
        self.setup_input_file("humidity_min", np.array([[63.0, 60.0]]))
        self.setup_input_file("wind_speed", np.array([[2.078, 2.244]]))
        self.setup_input_file("solar_radiation", np.array([[22.07, 21.62]]))

        # Also setup an output file that has no corresponding input files
        rogue_output_file = os.path.join(self.tempdir, "evaporation-2013-01-01.tif")
        with open(rogue_output_file, "w") as f:
            f.write("irrelevant contents")

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                elevation = 100
                step_length = 1440
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Verify the rogue output file is still here
        self.assertTrue(os.path.exists(rogue_output_file))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the rogue output file is gone
        self.assertFalse(os.path.exists(rogue_output_file))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-07-06")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(), np.array([[3.9, 4.8]]), decimal=1
        )
        fp = None

    def test_execute_hourly(self):
        # Prepare input files
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0, tzinfo=senegal_tzinfo)
        self.setup_input_file("temperature", np.array([[38.0, 28.0]]))
        self.setup_input_file("humidity", np.array([[52.0, 42.0]]))
        self.setup_input_file("wind_speed", np.array([[3.3, 2.3]]))
        self.setup_input_file("pressure", np.array([[1013.0, 1013.0]]))
        self.setup_input_file("solar_radiation", np.array([[681.0, 403.0]]))

        # Also setup an output file that has no corresponding input files
        rogue_output_file = os.path.join(
            self.tempdir, "evaporation-2013-01-01-15-00-0100.tif"
        )
        with open(rogue_output_file, "w") as f:
            f.write("irrelevant contents")

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                step_length = 60
                unit_converter_pressure = x / 10.0
                unit_converter_solar_radiation = x * 3600 / 1e6
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d-%H-%M%z")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Verify the rogue output file is still here
        self.assertTrue(os.path.exists(rogue_output_file))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the rogue output file is gone
        self.assertFalse(os.path.exists(rogue_output_file))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-10-01T15:00:00-01:00")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(), np.array([[0.63, 0.36]]), decimal=2
        )
        fp = None

    def test_execute_hourly_no_pressure(self):
        """Same as test_execute_hourly, but does not have pressure an input;
           therefore, it will calculate pressure itself."""
        # Prepare input files
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0, tzinfo=senegal_tzinfo)
        self.setup_input_file("temperature", np.array([[38.0, 28.0]]))
        self.setup_input_file("humidity", np.array([[52.0, 42.0]]))
        self.setup_input_file("wind_speed", np.array([[3.3, 2.3]]))
        self.setup_input_file("solar_radiation", np.array([[681.0, 403.0]]))

        # Also setup an output file that has no corresponding input files
        rogue_output_file = os.path.join(
            self.tempdir, "evaporation-2013-01-01-15-00-0100.tif"
        )
        with open(rogue_output_file, "w") as f:
            f.write("irrelevant contents")

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                step_length = 60
                unit_converter_solar_radiation = x * 3600 / 1e6
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d-%H-%M%z")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Verify the rogue output file is still here
        self.assertTrue(os.path.exists(rogue_output_file))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the rogue output file is gone
        self.assertFalse(os.path.exists(rogue_output_file))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-10-01T15:00:00-01:00")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(), np.array([[0.63, 0.36]]), decimal=2
        )
        fp = None

    def test_execute_hourly_without_sun(self):
        # Prepare input files, without solar radiation
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0, tzinfo=senegal_tzinfo)
        self.setup_input_file("temperature", np.array([[38.0, 28.0]]))
        self.setup_input_file("humidity", np.array([[52.0, 42.0]]))
        self.setup_input_file("wind_speed", np.array([[3.3, 2.3]]))
        self.setup_input_file("pressure", np.array([[1013.0, 1013.0]]))

        # Configuration
        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = 8
                step_length = 60
                unit_converter_pressure = x / 10.0
                unit_converter_solar_radiation = x * 3600 / 1e6
                """
                ).format(self=self)
            )

        # Execute and check exception
        msg = "Neither sunshine_duration nor solar_radiation are available"
        with self.assertRaisesRegex(click.ClickException, msg):
            cli.App(self.config_file).run()

    def test_execute_with_dem(self):
        """This is essentially the same as test_execute, but uses a GeoTIFF
        with a DEM instead of a constant elevation. The numbers are the same,
        however (all DEM gridpoints have the same value)."""

        # Prepare input files
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0, tzinfo=senegal_tzinfo)
        self.setup_input_file("temperature", np.array([[38.0, 28.0]]))
        self.setup_input_file("humidity", np.array([[52.0, 42.0]]))
        self.setup_input_file("wind_speed", np.array([[3.3, 2.3]]))
        self.setup_input_file("pressure", np.array([[1013.0, 1013.0]]))
        self.setup_input_file("solar_radiation", np.array([[681.0, 403.0]]))
        self.setup_input_file("dem", np.array([[8.0, 8.0]]), with_date=False)

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = dem.tif
                step_length = 60
                unit_converter_pressure = x / 10.0
                unit_converter_solar_radiation = x * 3600 / 1e6
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d-%H-%M%z")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-10-01T15:00:00-01:00")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(), np.array([[0.63, 0.36]]), decimal=2
        )
        fp = None

    def test_execute_with_nodata(self):
        """This is essentially the same as test_execute, but the gdal rasters
        contain cells with nodata."""

        sys.stderr = StringIO()

        # Prepare input files
        self.timestamp = dt.datetime(2014, 10, 1, 15, 0, tzinfo=senegal_tzinfo)
        nan = float("nan")
        self.setup_input_file("temperature", np.array([[38.0, nan]]))
        self.setup_input_file("humidity", np.array([[52.0, nan]]))
        self.setup_input_file("wind_speed", np.array([[3.3, 2.3]]))
        self.setup_input_file("pressure", np.array([[1013.0, nan]]))
        self.setup_input_file("solar_radiation", np.ma.array([[681.0, nan]]))
        self.setup_input_file("dem", np.array([[8.0, nan]]), with_date=False)

        with open(self.config_file, "w") as f:
            f.write(
                textwrap.dedent(
                    """\
                base_dir = {self.tempdir}
                albedo = 0.23
                nighttime_solar_radiation_ratio = 0.8
                elevation = dem.tif
                step_length = 60
                unit_converter_pressure = x / 10.0
                unit_converter_solar_radiation = x * 3600 / 1e6
                """
                ).format(self=self)
            )

        # Verify the output file doesn't exist yet
        result_filename = os.path.join(
            self.tempdir,
            "evaporation-{}.tif".format(self.timestamp.strftime("%Y-%m-%d-%H-%M%z")),
        )
        self.assertFalse(os.path.exists(result_filename))

        # Execute
        cli.App(self.config_file).run()

        # Check that it has created a file
        self.assertTrue(os.path.exists(result_filename))

        # Check that the created file is correct
        fp = gdal.Open(result_filename)
        timestamp = fp.GetMetadata()["TIMESTAMP"]
        self.assertEqual(timestamp, "2014-10-01T15:00:00-01:00")
        self.assertEqual(fp.RasterXSize, 2)
        self.assertEqual(fp.RasterYSize, 1)
        self.assertEqual(fp.GetGeoTransform(), self.geo_transform)
        # We can't just compare fp.GetProjection() to self.wgs84.ExportToWkt(),
        # because sometimes there are minor differences in the formatting or in
        # the information contained in the WKT.
        self.assertTrue(fp.GetProjection().startswith('GEOGCS["WGS 84",'))
        self.assertTrue(fp.GetProjection().endswith('AUTHORITY["EPSG","4326"]]'))
        nodatavalue = fp.GetRasterBand(1).GetNoDataValue()
        np.testing.assert_almost_equal(
            fp.GetRasterBand(1).ReadAsArray(),
            np.array([[0.63, nodatavalue]]),
            decimal=2,
        )

        # sys.stderr should be empty (no RuntimeWarning should have been
        # written there)
        self.assertEqual("", sys.stderr.getvalue())

        fp = None
