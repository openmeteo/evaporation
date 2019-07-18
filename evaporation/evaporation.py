import math
import warnings
from datetime import timedelta
from math import cos, pi, sin, tan

import numpy as np

# Note about RuntimeWarning
#
# When numpy makes calculations with masked arrays, it sometimes emits spurious
# RuntimeWarnings. This is because it occasionally does use the masked part of
# the array during the calculations (but masks the result). This is a known
# numpy bug (e.g. https://github.com/numpy/numpy/issues/4269). The numpy
# documentation, section "Operations on masked arrays", also has a related
# warning there.
#
# In order to avoid these spurious warnings, we have used, at various places in
# the code, "with warnings.catch_warnings()". We have attempted to unit test
# it, but sometimes it's hard to make the bug appear. A large array in
# production may cause the bug, but a small array in the unit test might not
# cause it, despite same python and numpy version. So the locations in which
# a fix was needed were largely located in production.


class PenmanMonteith(object):
    # Stefan-Boltzmann constant (Allen et al., 1998, p. 52)
    sigma = 4.903e-9

    def __init__(
        self,
        albedo,
        elevation,
        latitude,
        step_length,
        longitude=None,
        nighttime_solar_radiation_ratio=None,
        unit_converters={},
    ):
        self.albedo = albedo
        self.nighttime_solar_radiation_ratio = nighttime_solar_radiation_ratio
        self.elevation = elevation
        self.latitude = latitude
        self.longitude = longitude
        self.step_length = step_length
        self.unit_converters = unit_converters

    def calculate(self, **kwargs):
        if self.step_length == timedelta(minutes=60):
            return self.calculate_hourly(**kwargs)
        elif self.step_length == timedelta(days=1):
            return self.calculate_daily(**kwargs)
        else:
            raise NotImplementedError(
                "Evaporation for time steps other than hourly and daily "
                "has not been implemented."
            )

    def calculate_daily(
        self,
        temperature_max,
        temperature_min,
        humidity_max,
        humidity_min,
        wind_speed,
        adatetime,
        sunshine_duration=None,
        pressure=None,
        solar_radiation=None,
    ):
        if pressure is None:
            # Eq. 7 p. 31
            pressure = 101.3 * ((293 - 0.0065 * self.elevation) / 293) ** 5.26
        variables = self.convert_units(
            temperature_max=temperature_max,
            temperature_min=temperature_min,
            humidity_max=humidity_max,
            humidity_min=humidity_min,
            wind_speed=wind_speed,
            sunshine_duration=sunshine_duration,
            pressure=pressure,
        )

        # Radiation
        r_a, N = self.get_extraterrestrial_radiation(adatetime)
        if solar_radiation is None:
            solar_radiation = (
                0.25 + 0.50 * variables["sunshine_duration"] / N
            ) * r_a  # Eq.35 p. 50
        r_so = r_a * (0.75 + 2e-5 * self.elevation)  # Eq. 37, p. 51
        variables.update(self.convert_units(solar_radiation=solar_radiation))

        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore", RuntimeWarning)
            temperature_mean = (
                variables["temperature_max"] + variables["temperature_min"]
            ) / 2
        variables["temperature_mean"] = temperature_mean
        gamma = self.get_psychrometric_constant(temperature_mean, pressure)
        return self.penman_monteith_daily(
            incoming_solar_radiation=variables["solar_radiation"],
            clear_sky_solar_radiation=r_so,
            psychrometric_constant=gamma,
            mean_wind_speed=variables["wind_speed"],
            temperature_max=variables["temperature_max"],
            temperature_min=variables["temperature_min"],
            temperature_mean=variables["temperature_mean"],
            humidity_max=variables["humidity_max"],
            humidity_min=variables["humidity_min"],
            adate=adatetime,
        )

    def calculate_hourly(
        self,
        temperature,
        humidity,
        wind_speed,
        solar_radiation,
        adatetime,
        pressure=None,
    ):
        if pressure is None:
            # Eq. 7 p. 31
            pressure = 101.3 * ((293 - 0.0065 * self.elevation) / 293) ** 5.26
        variables = self.convert_units(
            temperature=temperature,
            humidity=humidity,
            wind_speed=wind_speed,
            pressure=pressure,
            solar_radiation=solar_radiation,
        )
        gamma = self.get_psychrometric_constant(
            variables["temperature"], variables["pressure"]
        )
        r_so = self.get_extraterrestrial_radiation(adatetime) * (
            0.75 + 2e-5 * self.elevation
        )  # Eq. 37, p. 51
        return self.penman_monteith_hourly(
            incoming_solar_radiation=variables["solar_radiation"],
            clear_sky_solar_radiation=r_so,
            psychrometric_constant=gamma,
            mean_wind_speed=variables["wind_speed"],
            mean_temperature=variables["temperature"],
            mean_relative_humidity=variables["humidity"],
            adatetime=adatetime,
        )

    def convert_units(self, **kwargs):
        result = {}
        for item in kwargs:
            varname = item
            if item.endswith("_max") or item.endswith("_min"):
                varname = item[:-4]
            converter = self.unit_converters.get(varname, lambda x: x)
            with warnings.catch_warnings():
                # See comment about RuntimeWarning on top of the file
                warnings.simplefilter("ignore", RuntimeWarning)
                result[item] = converter(kwargs[item])
        return result

    def get_extraterrestrial_radiation(self, adatetime):
        """
        Calculates the solar radiation we would receive if there were no
        atmosphere. This is a function of date, time and location.

        If adatetime is a datetime object, it merely returns the
        extraterrestrial radiation R_a; if it is a date object, it returns a
        tuple, (R_a, N), where N is the daylight hours.
        """
        j = adatetime.timetuple().tm_yday  # Day of year

        # Inverse relative distance Earth-Sun, eq. 23, p. 46.
        dr = 1 + 0.033 * cos(2 * pi * j / 365)

        # Solar declination, eq. 24, p. 46.
        decl = 0.409 * sin(2 * pi * j / 365 - 1.39)

        if self.step_length > timedelta(minutes=60):  # Daily?
            phi = self.latitude / 180.0 * pi
            omega_s = np.arccos(-np.tan(phi) * tan(decl))  # Eq. 25 p. 46

            r_a = (
                24
                * 60
                / pi
                * 0.0820
                * dr
                * (
                    omega_s * np.sin(phi) * sin(decl)
                    + np.cos(phi) * cos(decl) * np.sin(omega_s)
                )
            )  # Eq. 21 p. 46
            n = 24 / pi * omega_s  # Eq. 34 p. 48
            return r_a, n

        # Seasonal correction for solar time, eq. 32, p. 48.
        b = 2 * pi * (j - 81) / 364
        sc = 0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.025 * sin(b)

        # Longitude at the centre of the local time zone
        utc_offset = adatetime.utcoffset()
        utc_offset_hours = utc_offset.days * 24 + utc_offset.seconds / 3600.0
        lz = -utc_offset_hours * 15

        # Solar time angle at midpoint of the time period, eq. 31, p. 48.
        tm = adatetime - self.step_length / 2
        t = tm.hour + tm.minute / 60.0
        omega = pi / 12 * ((t + 0.06667 * (lz + self.longitude) + sc) - 12)

        # Solar time angles at beginning and end of the period, eqs. 29 and 30,
        # p. 48.
        t1 = self.step_length.seconds / 3600.0
        omega1 = omega - pi * t1 / 24
        omega2 = omega + pi * t1 / 24

        # Result: eq. 28, p. 47.
        phi = self.latitude / 180.0 * pi
        return (
            12
            * 60
            / pi
            * 0.0820
            * dr
            * (
                (omega2 - omega1) * np.sin(phi) * sin(decl)
                + np.cos(phi) * cos(decl) * (np.sin(omega2) - np.sin(omega1))
            )
        )

    def get_psychrometric_constant(self, temperature, pressure):
        """
        Allen et al. (1998), eq. 8, p. 32.

        This is called a "constant" because, although it is a function of
        temperature and pressure, its variations are small, and therefore it
        can be assumed constant for a location assuming standard pressure at
        that elevation and 20 degrees C. However, here we actually calculate
        it, so it isn't a constant.
        """
        lambda_ = 2.501 - (2.361e-3) * temperature  # eq. 3-1, p. 223
        return 1.013e-3 * pressure / 0.622 / lambda_

    def penman_monteith_daily(
        self,
        incoming_solar_radiation,
        clear_sky_solar_radiation,
        psychrometric_constant,
        mean_wind_speed,
        temperature_max,
        temperature_min,
        temperature_mean,
        humidity_max,
        humidity_min,
        adate,
    ):
        """
        Calculates and returns the reference evapotranspiration according
        to Allen et al. (1998), eq. 6, p. 24 & 65.
        """

        # Saturation and actual vapour pressure
        svp_max = self.get_saturation_vapour_pressure(temperature_max)
        svp_min = self.get_saturation_vapour_pressure(temperature_min)
        avp1 = svp_max * humidity_min / 100
        avp2 = svp_min * humidity_max / 100
        svp = (svp_max + svp_min) / 2  # Eq. 12 p. 36
        avp = (avp1 + avp2) / 2  # Eq. 12 p. 36

        # Saturation vapour pressure curve slope
        delta = self.get_saturation_vapour_pressure_curve_slope(temperature_mean)

        # Net incoming radiation; p. 51, eq. 38
        albedo = (
            self.albedo[adate.month - 1]
            if self.albedo.__class__.__name__ in ("tuple", "list")
            else self.albedo
        )
        rns = (1.0 - albedo) * incoming_solar_radiation

        # Net outgoing radiation
        rnl = self.get_net_outgoing_radiation(
            (temperature_min, temperature_max),
            incoming_solar_radiation,
            clear_sky_solar_radiation,
            avp,
        )

        # Net radiation at grass surface
        rn = rns - rnl

        # Soil heat flux
        g_day = 0  # Eq. 42 p. 54

        # Apply the formula
        numerator_term1 = 0.408 * delta * (rn - g_day)
        numerator_term2 = (
            psychrometric_constant
            * 900
            / (temperature_mean + 273.16)
            * mean_wind_speed
            * (svp - avp)
        )
        denominator = delta + psychrometric_constant * (1 + 0.34 * mean_wind_speed)

        return (numerator_term1 + numerator_term2) / denominator

    def penman_monteith_hourly(
        self,
        incoming_solar_radiation,
        clear_sky_solar_radiation,
        psychrometric_constant,
        mean_wind_speed,
        mean_temperature,
        mean_relative_humidity,
        adatetime,
    ):
        """
        Calculates and returns the reference evapotranspiration according
        to Allen et al. (1998), eq. 53, p. 74.

        As explained in Allen et al. (1998, p. 74), the function is
        modified in relation to the original Penman-Monteith equation, so
        that it is suitable for hourly data.
        """

        # Saturation and actual vapour pressure
        svp = self.get_saturation_vapour_pressure(mean_temperature)
        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore", RuntimeWarning)
            avp = svp * mean_relative_humidity / 100.0  # Eq. 54, p. 74

        # Net incoming radiation; p. 51, eq. 38
        albedo = (
            self.albedo[adatetime.month - 1]
            if self.albedo.__class__.__name__ in ("tuple", "list")
            else self.albedo
        )
        rns = (1.0 - albedo) * incoming_solar_radiation

        # Net outgoing radiation
        rnl = self.get_net_outgoing_radiation(
            mean_temperature, incoming_solar_radiation, clear_sky_solar_radiation, avp
        )

        # Net radiation at grass surface
        rn = rns - rnl

        # Saturation vapour pressure curve slope
        delta = self.get_saturation_vapour_pressure_curve_slope(mean_temperature)

        # Soil heat flux density
        g = self.get_soil_heat_flux_density(incoming_solar_radiation, rn)

        # Apply the formula
        numerator_term1 = 0.408 * delta * (rn - g)
        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore", RuntimeWarning)
            numerator_term2 = (
                psychrometric_constant
                * 37
                / (mean_temperature + 273.16)
                * mean_wind_speed
                * (svp - avp)
            )
        denominator = delta + psychrometric_constant * (1 + 0.34 * mean_wind_speed)

        return (numerator_term1 + numerator_term2) / denominator

    def get_net_outgoing_radiation(
        self,
        temperature,
        incoming_solar_radiation,
        clear_sky_solar_radiation,
        mean_actual_vapour_pressure,
    ):
        """
        Allen et al. (1998), p. 52, eq. 39. Temperature can be a tuple (a pair)
        of min and max, or a single value. If it is a single value, the
        equation is modified according to end of page 74.
        """
        if temperature.__class__.__name__ in ("tuple", "list"):
            with warnings.catch_warnings():
                # See comment about RuntimeWarning on top of the file
                warnings.simplefilter("ignore", RuntimeWarning)
                factor1 = (
                    self.sigma
                    * ((temperature[0] + 273.16) ** 4 + (temperature[1] + 273.16) ** 4)
                    / 2
                )
        else:
            with warnings.catch_warnings():
                # See comment about RuntimeWarning on top of the file
                warnings.simplefilter("ignore", RuntimeWarning)
                factor1 = self.sigma / 24 * (temperature + 273.16) ** 4
        factor2 = 0.34 - 0.14 * (mean_actual_vapour_pressure ** 0.5)

        # Solar radiation ratio Rs/Rs0 (Allen et al., 1998, top of p. 75).
        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore", RuntimeWarning)
            solar_radiation_ratio = np.where(
                clear_sky_solar_radiation > 0.05,
                incoming_solar_radiation / clear_sky_solar_radiation,
                self.nighttime_solar_radiation_ratio,
            )
        solar_radiation_ratio = np.maximum(solar_radiation_ratio, 0.3)
        solar_radiation_ratio = np.minimum(solar_radiation_ratio, 1.0)

        factor3 = 1.35 * solar_radiation_ratio - 0.35

        return factor1 * factor2 * factor3

    def get_saturation_vapour_pressure(self, temperature):
        "Allen et al. (1998), p. 36, eq. 11."
        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore")
            return 0.6108 * math.e ** (17.27 * temperature / (237.3 + temperature))

    def get_soil_heat_flux_density(self, incoming_solar_radiation, rn):
        "Allen et al. (1998), p. 55, eq. 45 & 46."
        coefficient = np.where(incoming_solar_radiation > 0.05, 0.1, 0.5)
        return coefficient * rn

    def get_saturation_vapour_pressure_curve_slope(self, temperature):
        "Allen et al. (1998), p. 37, eq. 13."
        numerator = 4098 * self.get_saturation_vapour_pressure(temperature)
        with warnings.catch_warnings():
            # See comment about RuntimeWarning on top of the file
            warnings.simplefilter("ignore", RuntimeWarning)
            denominator = (temperature + 237.3) ** 2
            return numerator / denominator
