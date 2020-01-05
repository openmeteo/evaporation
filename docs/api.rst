.. _api:

===
API
===

``from evaporation import PenmanMonteith``

.. class:: PenmanMonteith(albedo, elevation, latitude, time_step, longitude=None, nighttime_solar_radiation_ratio=None, unit_converters={})

   Calculates evapotranspiration according to the Penman-Monteith
   equation. The methodology used is that of Allen et al. (1998).
   Details can be found in the code, which has comments indicating
   which equations it uses.

   First the class is initialized with some parameters that are
   constant for the area of interest; then, the :meth:`calculate`
   method can be called as many times as necessary in order to
   calculate reference evapotranspiration given the date and time and
   the values of the meteorological variables.

   Evapotranspiration can be calculated either at a point or on a
   grid, so the input can be either simple scalar values or ``numpy``
   arrays. The term "scalar or array", when used below, signifies a
   parameter that can be either. You will normally either use all
   scalars or all arrays; however, when you generally use arrays, you
   may also use scalars for some of the parameters if you want the
   array to have the same value for all gridpoints; for example, you
   might want to have a single albedo value for all gridpoints.

   The class is initialized with the following parameters:

   *albedo* is either a scalar or array, or a sequence of 12 scalars
   or arrays. If it is a sequence, the first item is the albedo in
   January, the second is for February, and so on. If it is a single
   scalar or array, it is used for the entire year. The albedo is a
   number between 0 and 1.

   *elevation* is a scalar or array with the location elevation above
   sea level in meters.

   *latitude* and *longitude* are scalars or arrays, in decimal
   degrees north of the equator or east of the prime meridian
   (negative for west or south). Only *latitude* needs to be specified
   for calculating daily evaporation.

   *time_step* is a string: "D" for daily, "H" for hourly.

   In order to estimate the outgoing radiation, the ratio of incoming
   solar radiation to clear sky solar radiation is used as a
   representation of cloud cover. However, when calculating hourly
   evaporation, this does not work during the night, in which case
   *nighttime_solar_radiation_ratio* is used as a rough approximation
   of that ratio. It should be a scalar or array between 0.4 and 0.8;
   see Allen et al. (1998), top of page 75.

   The meteorological values that will be supplied after class
   initialization to the :meth:`calculate` method are supposed to be
   in the following units:
   
   ========================  =====================
   Parameter                 Unit
   ========================  =====================
   temperature               ℃
   humidity                  %
   wind speed                m/s
   pressure                  kPa
   solar radiation           MJ/m²/h
   ========================  =====================
   
   If they are in different units, *unit_converters* is a dictionary
   with functions to convert them. For example, if you have pressure 
   in hPa and solar radiation in W/m², you should specify this::

      unit_converters = {
          'pressure': lambda x: x / 10.0,
          'solar_radiation': lambda x: x * 3600 / 1e6,
      }

   Any variable whose name is not found in *unit_converters* is used
   as is, without conversion.

   .. method:: calculate(self, **kwargs)

      Calculates and returns the reference evapotranspiration in mm.

      For daily step, the keyword arguments must be *temperature_max*,
      *temperature_min*, *humidity_max*, *humidity_min*, *wind_speed*,
      *adatetime*, and one of *solar_radiation* or *sunshine_duration*.
      *adatetime* must be a :class:`~datetime.date` object, not a
      :class:`~datetime.datetime` object, but it is named *adatetime* for
      consistency with the hourly step. The result is the reference
      evapotranspiration for the given day.

      For hourly step, the keyword arguments must be *temperature*,
      *humidity*, *wind_speed*, *solar_radiation*, *adatetime*, and,
      optionally, *pressure* (if the pressure is not specified it is
      calculated from the elevation). The result is the reference
      evapotranspiration for the hour that ends at *adatetime*, which
      must be a timezone-aware :class:`~datetime.datetime` object.



References
----------

R. G. Allen, L. S. Pereira, D. Raes, and M. Smith, Crop evapotranspiration -
Guidelines for computing crop water requirements, FAO Irrigation and drainage
paper no. 56, 1998.
