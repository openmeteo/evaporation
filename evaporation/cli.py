import configparser
import datetime as dt
import logging
import os
import sys
import traceback
from glob import glob
from io import StringIO

import click
import iso8601
import numpy as np
from hspatial import NODATAVALUE
from htimeseries import HTimeseries, TzinfoFromString
from osgeo import gdal, ogr, osr

from . import __version__
from .evaporation import PenmanMonteith

gdal.UseExceptions()


class WrongValueError(configparser.Error):
    pass


class App:
    def __init__(self, configfilename):
        self.configfilename = configfilename

    def run(self):
        self.config = AppConfig(self.configfilename)
        self.config.read()
        self._setup_logger()
        self._execute_with_error_handling()

    def _execute_with_error_handling(self):
        self.logger.info("Starting evaporation, " + dt.datetime.today().isoformat())
        try:
            self._execute()
        except Exception as e:
            self.logger.error(str(e))
            self.logger.debug(traceback.format_exc())
            self.logger.info(
                "evaporation terminated with error, " + dt.datetime.today().isoformat()
            )
            raise click.ClickException(str(e))
        else:
            self.logger.info("Finished evaporation, " + dt.datetime.today().isoformat())

    def _setup_logger(self):
        self.logger = logging.getLogger("evaporation")
        self._set_logger_handler()
        self.logger.setLevel(self.config.loglevel.upper())

    def _set_logger_handler(self):
        if getattr(self.config, "logfile", None):
            self.logger.addHandler(logging.FileHandler(self.config.logfile))
        else:
            self.logger.addHandler(logging.StreamHandler())

    def _execute(self):
        if self.config.is_spatial:
            ProcessSpatial(self.config).execute()
        else:
            ProcessAtPoint(self.config).execute()


class AppConfig:
    config_file_options = {
        "logfile": {"fallback": ""},
        "loglevel": {"fallback": "warning"},
        "base_dir": {"fallback": ""},
        "step_length": {},
        "elevation": {"fallback": ""},
        "albedo": {},
        "nighttime_solar_radiation_ratio": {"fallback": 0.6},
        "unit_converter_temperature": {"fallback": "x"},
        "unit_converter_humidity": {"fallback": "x"},
        "unit_converter_wind_speed": {"fallback": "x"},
        "unit_converter_pressure": {"fallback": "x"},
        "unit_converter_solar_radiation": {"fallback": "x"},
        "temperature_min_prefix": {"fallback": "temperature_min"},
        "temperature_max_prefix": {"fallback": "temperature_max"},
        "temperature_prefix": {"fallback": "temperature"},
        "humidity_prefix": {"fallback": "humidity"},
        "humidity_min_prefix": {"fallback": "humidity_min"},
        "humidity_max_prefix": {"fallback": "humidity_max"},
        "wind_speed_prefix": {"fallback": "wind_speed"},
        "pressure_prefix": {"fallback": "pressure"},
        "solar_radiation_prefix": {"fallback": "solar_radiation"},
        "sunshine_duration_prefix": {"fallback": "sunshine_duration"},
        "evaporation_prefix": {"fallback": "evaporation"},
    }

    def __init__(self, configfilename):
        self.configfilename = configfilename

    def read(self):
        try:
            self._parse_config()
        except (OSError, configparser.Error) as e:
            sys.stderr.write(str(e))
            raise click.ClickException(str(e))

    def _parse_config(self):
        self._read_config_file()
        self._get_config_options()
        self._parse_config_options()
        self._parse_unit_converters()
        self._parse_elevation()

    def _read_config_file(self):
        self.config = configparser.ConfigParser(interpolation=None)
        try:
            self._read_config_file_assuming_it_has_section_headers()
        except configparser.MissingSectionHeaderError:
            self._read_config_file_without_sections()

    def _read_config_file_assuming_it_has_section_headers(self):
        with open(self.configfilename) as f:
            self.config.read_file(f)

    def _read_config_file_without_sections(self):
        with open(self.configfilename) as f:
            configuration = "[General]\n" + f.read()
        self.config.read_file(StringIO(configuration))

    def _get_config_options(self):
        self.options = {
            opt: self.config.get("General", opt, **kwargs)
            for opt, kwargs in self.config_file_options.items()
        }
        for key, value in self.options.items():
            setattr(self, key, value)

    def _parse_config_options(self):
        self._parse_log_level()
        self._parse_step_length()
        self._parse_unit_converters()
        self._set_is_spatial()
        self._parse_elevation()
        self._parse_albedo()
        self._parse_nighttime_solar_radiation_ratio()

    def _parse_log_level(self):
        log_levels = ("ERROR", "WARNING", "INFO", "DEBUG")
        self.loglevel = self.options["loglevel"].upper()
        if self.loglevel not in log_levels:
            raise WrongValueError("loglevel must be one of " + ", ".join(log_levels))

    def _parse_step_length(self):
        s = self.options["step_length"]
        self._check_step_length(s)
        self.step = dt.timedelta(minutes=int(s))

    def _check_step_length(self, s):
        if not s.isdecimal() or int(s) not in (60, 1440):
            raise WrongValueError(
                '"{}" is not an appropriate time step; in this version of '
                "vaporize, the step must be either 60 or 1440.".format(s)
            )

    def _parse_unit_converters(self):
        vars = {"temperature", "humidity", "wind_speed", "pressure", "solar_radiation"}
        self.unit_converters = {}
        for var in vars:
            config_item = "unit_converter_" + var
            config_value = self.options[config_item]
            lambda_defn = "lambda x: " + config_value
            self._set_unit_converter(var, config_item, config_value, lambda_defn)

    def _set_unit_converter(self, var, config_item, config_value, lambda_defn):
        try:
            self.unit_converters[var] = eval(lambda_defn)
        except Exception as e:
            raise WrongValueError(
                "{} while parsing {} ({}): {}".format(
                    e.__class__.__name__, config_item, config_value, str(e)
                )
            )

    def _parse_elevation(self):
        s = self.options["elevation"]
        if s == "":
            self.elevation = None
        else:
            self.elevation = self._get_number_or_grid(s)
        self._check_elevation()

    def _check_elevation(self):
        if self.is_spatial:
            self._check_elevation_for_spatial()
        else:
            self._check_elevation_for_point()

    def _check_elevation_for_spatial(self):
        if self.elevation is None:
            raise WrongValueError(
                "elevation needs to be specified in the configuration file"
            )
        elif np.any(self.elevation < -427) or np.any(self.elevation > 8848):
            raise WrongValueError("The elevation must be between -427 and 8848")

    def _check_elevation_for_point(self):
        if self.elevation is not None:
            raise WrongValueError(
                "elevation should only be specified in the configuration file "
                "in spatial calculation (otherwise we get it from the hts files)"
            )

    def _parse_albedo(self):
        s = self.options["albedo"].split()
        self._check_albedo_is_one_or_twelve_items(s)
        albedo = []
        for item in s:
            albedo.append(self._get_number_or_grid(item))
            self._check_albedo_domain(albedo[-1])
        self.albedo = albedo[0] if len(s) == 1 else albedo

    def _check_albedo_is_one_or_twelve_items(self, s):
        if len(s) not in (1, 12):
            raise ValueError(
                "Albedo must be either one item or 12 space-separated items"
            )

    def _check_albedo_domain(self, albedo):
        value_to_test = albedo if isinstance(albedo, float) else albedo.all()
        if value_to_test < 0.0 or value_to_test > 1.0:
            raise ValueError("Albedo must be between 0.0 and 1.0")

    def _parse_nighttime_solar_radiation_ratio(self):
        s = self.options["nighttime_solar_radiation_ratio"]
        self.nighttime_solar_radiation_ratio = self._get_number_or_grid(s)
        self._check_nighttime_solar_radiation_ratio()

    def _check_nighttime_solar_radiation_ratio(self):
        a = self.nighttime_solar_radiation_ratio
        if a < 0.4 or a > 0.8:
            raise WrongValueError(
                "The nighttime solar radiation ratio must " "be between 0.4 and 0.8"
            )

    def _get_number_or_grid(self, s):
        """Return either a number or a grid from s, depending on its contents.

        If string s holds a valid number, return it; otherwise try to open the
        geotiff file whose filename is s and read its first band into the
        returned numpy array.
        """
        try:
            return float(s)
        except ValueError:
            return self._get_grid(s)

    def _get_grid(self, s):
        input_file = gdal.Open(os.path.join(self.options["base_dir"], s))
        result = input_file.GetRasterBand(1).ReadAsArray()
        nodata = input_file.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            result = np.ma.masked_values(result, nodata, copy=False)
        return result

    def _set_is_spatial(self):
        base_dir = self.options["base_dir"]
        base_dir_has_tif_files = bool(glob(os.path.join(base_dir, "*.tif")))
        base_dir_has_hts_files = bool(glob(os.path.join(base_dir, "*.hts")))
        self._check_tif_hts_consistency(base_dir_has_tif_files, base_dir_has_hts_files)
        self.is_spatial = base_dir_has_tif_files

    def _check_tif_hts_consistency(self, has_tif, has_hts):
        if has_tif and has_hts:
            raise WrongValueError(
                "Base directory {} contains both tif files and hts files; "
                "this is not allowed.".format(self.options["base_dir"])
            )
        elif not has_tif and not has_hts:
            raise WrongValueError(
                "Base directory {} contains neither tif files nor hts files.".format(
                    self.options["base_dir"]
                )
            )


class ProcessAtPoint:
    def __init__(self, config):
        self.config = config

    def execute(self):
        self._read_input_timeseries()
        self._setup_attrs()
        self._check_all_timeseries_are_in_same_location_and_timezone()
        self._get_location_in_wgs84()
        self._prepare_penman_monteith_parameters()
        self._prepare_resulting_htimeseries_object()
        self._determine_variables_to_use_in_calculation()
        self._calculate_evaporation()
        self._save_result()

    def _setup_attrs(self):
        atimeseries = self.input_timeseries["wind_speed"]
        self.location = atimeseries.location
        self.timezone = getattr(atimeseries, "timezone", None)

    def _read_input_timeseries(self):
        vars = {
            "temperature_min",
            "temperature_max",
            "temperature",
            "humidity",
            "humidity_min",
            "humidity_max",
            "wind_speed",
            "pressure",
            "solar_radiation",
            "sunshine_duration",
        }
        self.input_timeseries = {}
        for var in vars:
            self._get_input_timeseries_for_var(var)

    def _get_input_timeseries_for_var(self, var):
        filename = os.path.join(
            self.config.base_dir, getattr(self.config, var + "_prefix") + ".hts"
        )
        if not os.path.exists(filename):
            return
        with open(filename, "r") as f:
            self.input_timeseries[var] = HTimeseries(f)

    def _check_all_timeseries_are_in_same_location_and_timezone(self):
        for i, (name, hts) in enumerate(self.input_timeseries.items()):
            if i == 0:
                reference_hts = hts
            else:
                self._compare_location_and_timezone_of(hts, reference_hts)

    def _compare_location_and_timezone_of(self, hts1, hts2):
        self._compare_locations_of(hts1, hts2)
        self._compare_altitudes_of(hts1, hts2)
        self._compare_timezones_of(hts1, hts2)

    def _compare_locations_of(self, hts1, hts2):
        loc1 = hts1.location
        loc2 = hts2.location

        abscissas_differ = self._numbers_are_wrong_or_differ(
            loc1.get("abscissa"), loc2.get("abscissa")
        )
        ordinates_differ = self._numbers_are_wrong_or_differ(
            loc1.get("ordinate"), loc2.get("ordinate")
        )
        srids_differ = loc1.get("srid") != loc2.get("srid")

        if abscissas_differ or ordinates_differ or srids_differ:
            raise ValueError(
                "Incorrect or unspecified or inconsistent locations in the time series "
                "files."
            )

    def _numbers_are_wrong_or_differ(self, num1, num2, tolerance=1e7):
        if num1 is None or num2 is None:
            return True
        return abs(num1 - num2) > tolerance

    def _compare_altitudes_of(self, hts1, hts2):
        altitude1 = hts1.location.get("altitude")
        altitude2 = hts2.location.get("altitude")
        if self._numbers_are_wrong_or_differ(altitude1, altitude2, 1e-2):
            raise ValueError(
                "Incorrect or unspecified or inconsistent altitudes in the time series "
                "files."
            )

    def _compare_timezones_of(self, hts1, hts2):
        timezone1 = getattr(hts1, "timezone", "")
        timezone2 = getattr(hts2, "timezone", "")
        if timezone1 != timezone2:
            raise ValueError(
                "Incorrect or unspecified or inconsistent time zones in the time "
                "series files."
            )

    def _get_location_in_wgs84(self):
        source_projection = osr.SpatialReference()
        source_projection.ImportFromEPSG(self.location["srid"])
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(source_projection, wgs84)
        apoint = ogr.Geometry(ogr.wkbPoint)
        apoint.AddPoint(self.location["abscissa"], self.location["ordinate"])
        apoint.Transform(transform)
        self.location["latitude"] = apoint.GetY()
        self.location["longitude"] = apoint.GetX()

    def _prepare_penman_monteith_parameters(self):
        nsrr = self.config.nighttime_solar_radiation_ratio
        self.penman_monteith = PenmanMonteith(
            albedo=self.config.albedo,
            nighttime_solar_radiation_ratio=nsrr,
            elevation=self.location["altitude"],
            latitude=self.location["latitude"],
            longitude=self.location["longitude"],
            step_length=self.config.step,
            unit_converters=self.config.unit_converters,
        )

    def _prepare_resulting_htimeseries_object(self):
        self.pet = HTimeseries()
        minutes = int(self.config.step.total_seconds() / 60)
        self.pet.time_step = str(minutes) + ",0"
        self.pet.unit = "mm"
        self.pet.timezone = self.timezone
        self.pet.variable = "Potential Evapotranspiration"
        self.pet.precision = 2 if self.config.step == dt.timedelta(hours=1) else 1
        self.pet.location = self.location

    def _determine_variables_to_use_in_calculation(self):
        if self.config.step == dt.timedelta(hours=1):
            vars = ["temperature", "humidity", "wind_speed", "solar_radiation"]
            if "pressure" in self.input_timeseries:
                vars.append("pressure")
        elif self.config.step == dt.timedelta(days=1):
            vars = (
                "temperature_max",
                "temperature_min",
                "humidity_max",
                "humidity_min",
                "wind_speed",
                (
                    "solar_radiation"
                    if "solar_radiation" in self.input_timeseries
                    else "sunshine_duration"
                ),
            )
        self.input_vars = vars

    def _calculate_evaporation(self):
        for adatetime in self.input_timeseries["wind_speed"].data.index:
            self._calculate_evaporation_for(adatetime)

    def _calculate_evaporation_for(self, adatetime):
        try:
            kwargs = {
                v: self.input_timeseries[v].data.loc[adatetime, "value"]
                for v in self.input_vars
            }
        except (IndexError, KeyError):
            return
        kwargs["adatetime"] = self._datetime64_to_aware_datetime(adatetime)
        self.pet.data.loc[adatetime, "value"] = self.penman_monteith.calculate(**kwargs)

    def _datetime64_to_aware_datetime(self, adatetime):
        result = adatetime.to_pydatetime()
        if self.timezone:
            result = result.replace(tzinfo=TzinfoFromString(self.timezone))
        return result

    def _save_result(self):
        outfilename = self.config.evaporation_prefix + ".hts"
        outpathname = os.path.join(self.config.base_dir, outfilename)
        with open(outpathname, "w") as f:
            self.pet.write(f, format=HTimeseries.FILE)


class ProcessSpatial:
    def __init__(self, config):
        self.config = config

    def execute(self):
        self._get_geographical_reference_file()
        self._get_timestamps()
        self._get_geodata()
        self._prepare_penman_monteith_parameters()
        self._calculate_evaporation_for_all_timestamps()
        self._remove_stale_evaporation_files()

    def _get_geographical_reference_file(self):
        # Arbitrarily use the first temperature file to extract location and
        # other geographical stuff. Elsewhere consistency of such data from all
        # other files with this file will be checked.
        pattern = os.path.join(
            self.config.base_dir, self.config.wind_speed_prefix + "-*.tif"
        )
        wind_speed_files = glob(pattern)
        self.geographical_reference_file = wind_speed_files[0]

    def _get_timestamps(self):
        pattern = self.config.wind_speed_prefix + "-*.tif"
        saved_cwd = os.getcwd()
        try:
            os.chdir(self.config.base_dir)
            wind_speed_files = glob(pattern)
        finally:
            os.chdir(saved_cwd)

        # Remove the prefix from the start and the .tif from the end, leaving
        # only the date.
        prefix_len = len(self.config.wind_speed_prefix)
        start = prefix_len + 1
        self.timestamps = [item[start:-4] for item in wind_speed_files]

    def _get_geodata(self):
        """
        Retrieve geographical stuff into self.latitude, self.longitude,
        self.width, self.height, self.geo_transform, self.projection. We do
        this by retrieving the data from self.geographical_reference_file.
        """
        # Read data from GeoTIFF file
        fp = gdal.Open(self.geographical_reference_file)
        self.width, self.height = fp.RasterXSize, fp.RasterYSize
        self.geo_transform = fp.GetGeoTransform()
        self.projection = osr.SpatialReference()
        self.projection.ImportFromWkt(fp.GetProjection())

        # Find (x_left, y_top), (x_right, y_bottom)
        x_left, x_step, d1, y_top, d2, y_step = self.geo_transform
        x_right = x_left + self.width * x_step
        y_bottom = y_top + self.height * y_step

        # Transform into (long_left, lat_top), (long_right, lat_bottom)
        wgs84 = osr.SpatialReference()
        wgs84.ImportFromEPSG(4326)
        transform = osr.CoordinateTransformation(self.projection, wgs84)
        top_left = ogr.Geometry(ogr.wkbPoint)
        top_left.AddPoint(x_left, y_top)
        bottom_right = ogr.Geometry(ogr.wkbPoint)
        bottom_right.AddPoint(x_right, y_bottom)
        top_left.Transform(transform)
        bottom_right.Transform(transform)
        long_left, lat_top = top_left.GetX(), top_left.GetY()
        long_right, lat_bottom = bottom_right.GetX(), bottom_right.GetY()

        # Calculate self.latitude and self.longitude
        long_step = (long_right - long_left) / self.width
        longitudes = np.arange(
            long_left + long_step / 2.0, long_left + long_step * self.width, long_step
        )
        lat_step = (lat_top - lat_bottom) / self.height
        latitudes = np.arange(
            lat_top + lat_step / 2.0, lat_top + lat_step * self.height, lat_step
        )
        self.longitude, self.latitude = np.meshgrid(longitudes, latitudes)

    def _prepare_penman_monteith_parameters(self):
        nsrr = self.config.nighttime_solar_radiation_ratio
        self.penman_monteith = PenmanMonteith(
            albedo=self.config.albedo,
            nighttime_solar_radiation_ratio=nsrr,
            elevation=self.config.elevation,
            latitude=self.latitude,
            longitude=self.longitude,
            step_length=self.config.step,
            unit_converters=self.config.unit_converters,
        )

    def _calculate_evaporation_for_all_timestamps(self):
        for timestamp in self.timestamps:
            self._calculate_evaporation_for(timestamp)

    def _calculate_evaporation_for(self, timestamp):
        self._read_input_geofiles(timestamp)
        self._verify_that_solar_radiation_or_sunshine_duration_was_present(timestamp)
        adatetime = self._get_adatetime_arg(timestamp)
        result = self.penman_monteith.calculate(adatetime=adatetime, **self.input_data)
        self._write_result_to_file(result, timestamp, adatetime)

    def _read_input_geofiles(self, timestamp):
        self.input_data = {v: None for v in self._get_variables()}
        for var in self.input_data:
            self._read_input_geofile(var, timestamp)

    def _get_variables(self):
        if self.config.step == dt.timedelta(minutes=60):
            return {
                "temperature",
                "humidity",
                "wind_speed",
                "pressure",
                "solar_radiation",
            }
        else:
            return {
                "temperature_max",
                "temperature_min",
                "humidity_max",
                "humidity_min",
                "wind_speed",
                "solar_radiation",
                "sunshine_duration",
            }

    def _read_input_geofile(self, variable, timestamp):
        fp = self._open_geofile(variable, timestamp)
        if fp is None:
            return
        self._verify_consistency_of_geodata(fp)
        self.input_data[variable] = self._read_array_from_geofile(fp)
        fp = None

    def _open_geofile(self, variable, timestamp):
        filename_prefix = getattr(self.config, variable + "_prefix")
        filename = filename_prefix + "-" + timestamp + ".tif"
        self._filename = os.path.join(self.config.base_dir, filename)
        if variable in ("solar_radiation", "sunshine_duration"):
            if not os.path.exists(self._filename):
                # Either solar_radiation or sunshine_duration may be absent; here we
                # allow both to be absent and elsewhere we will check that one was
                # present
                return
        try:
            return gdal.Open(self._filename)
        except RuntimeError:
            if variable == "pressure":
                return
            raise

    def _verify_consistency_of_geodata(self, fp):
        consistent = all(
            (
                self.width == fp.RasterXSize,
                self.height == fp.RasterYSize,
                self.geo_transform == fp.GetGeoTransform(),
                self.projection.ExportToWkt() == fp.GetProjection(),
            )
        )
        if not consistent:
            raise Exception(
                "Not all input files have the same "
                "width, height, geo_transform and projection "
                "(offending items: {} and {})".format(
                    self.geographical_reference_file, self._filename
                )
            )

    def _read_array_from_geofile(self, fp):
        array = fp.GetRasterBand(1).ReadAsArray()
        nodata = fp.GetRasterBand(1).GetNoDataValue()
        if nodata is not None:
            array = np.ma.masked_values(array, nodata, copy=False)
        return array

    def _verify_that_solar_radiation_or_sunshine_duration_was_present(self, timestamp):
        if self.input_data["solar_radiation"] is not None:
            self.input_data.pop("sunshine_duration", None)
        elif self.input_data.get("sunshine_duration", None) is not None:
            self.input_data.pop("solar_radiation", None)
        else:
            raise RuntimeError(
                "Neither sunshine_duration nor solar_radiation are available for "
                + timestamp
            )

    def _get_adatetime_arg(self, timestamp):
        adatetime = iso8601.parse_date(
            self._timestamp_from_filename(timestamp), default_timezone=None
        )
        if self.config.step == dt.timedelta(minutes=1440):
            adatetime = adatetime.date()
        elif adatetime.tzinfo is None:
            raise Exception(
                "The time stamp in the input files does not have a time zone specified."
            )
        return adatetime

    def _timestamp_from_filename(self, s):
        """Convert a timestamp from its filename format  to its iso format

        E.g. from 2014-10-01-15-00-0100 to 2014-10-01 15:00-0100).
        """
        first_hyphen = s.find("-")
        if first_hyphen < 0:
            return s
        second_hyphen = s.find("-", first_hyphen + 1)
        if second_hyphen < 0:
            return s
        third_hyphen = s.find("-", second_hyphen + 1)
        if third_hyphen < 0:
            return s
        fourth_hyphen = s.find("-", third_hyphen + 1)
        chars = list(s)
        chars[third_hyphen] = " "
        if fourth_hyphen > 0:
            chars[fourth_hyphen] = ":"
        return "".join(chars)

    def _write_result_to_file(self, result, timestamp, adatetime):
        output_filename = self.config.evaporation_prefix + "-" + timestamp + ".tif"
        output_pathname = os.path.join(self.config.base_dir, output_filename)
        output = gdal.GetDriverByName("GTiff").Create(
            output_pathname, self.width, self.height, 1, gdal.GDT_Float32
        )
        try:
            output.SetMetadataItem("TIMESTAMP", adatetime.isoformat())
            output.SetGeoTransform(self.geo_transform)
            output.SetProjection(self.projection.ExportToWkt())
            result[result.mask] = NODATAVALUE
            output.GetRasterBand(1).SetNoDataValue(NODATAVALUE)
            output.GetRasterBand(1).WriteArray(result)
        finally:
            # Close the dataset
            output = None

    def _remove_stale_evaporation_files(self):
        """Remove evaporation files for which no input files exist.
        """
        pattern = self.config.evaporation_prefix + "-*.tif"
        saved_cwd = os.getcwd()
        try:
            os.chdir(self.config.base_dir)
            evaporation_files = glob(pattern)
            prefix_len = len(self.config.evaporation_prefix)
            for filename in evaporation_files:
                start = prefix_len + 1
                if filename[start:-4] not in self.timestamps:
                    os.unlink(filename)
        finally:
            os.chdir(saved_cwd)


@click.command()
@click.argument("configfile")
@click.version_option(version=__version__, prog_name="vaporize")
def main(configfile):
    """Calculation of evaporation and transpiration"""

    app = App(configfile)
    app.run()
