=====
Usage
=====

Synopsis
========

``vaporize [--traceback] config_file``

Description and quick start
===========================

``vaporize`` calculates evapotranspiration with the Penman-Monteith
method. It works either with GeoTIFF files or with time series files.
In either case, it reads files with temperature, humidity, solar
radiation, pressure and wind speed, and produces a file with
evapotranspiration. The details of its operation are specified in the
configuration file specified on the command line.

The methodology used is that of Allen et al. (1998).  Details can be
found in :ref:`api` and in the code itself, which has comments
indicating which equations it uses.

Installation
------------

``pip install evaporation``

How to run it
-------------

First, you need to create a configuration file with a text editor such
as ``vim``, ``emacs``, ``notepad``, or whatever. Create such a file
and name it, for example, :file:`/var/tmp/vaporize.conf`, with the
following contents (the contents don't matter at this stage, just copy
and paste them from below)::

    loglevel = INFO

Then, open a command prompt and give it this command::

    vaporize /var/tmp/vaporize.conf

If you have done everything correctly, it should output an error message
complaining that something in its configuration file isn't right.

Configuration file example
--------------------------

Take a look at the following example configuration file and read the
explanatory comments that follow it:

.. code-block:: ini

   loglevel = INFO
   logfile = C:\Somewhere\vaporize.log
   base_dir = C:\Somewhere
   albedo = 0.23
   nighttime_solar_radiation_ratio = 0.8
   elevation = 8
   step_length = 60
   unit_converter_pressure = x / 10.0
   unit_converter_solar_radiation = x * 3600 / 1e6

With the above configuration file, ``vaporize`` will log information in
the file specified by :confval:`logfile`. It will calculate hourly
evaporation (:confval:`step_length`) at the specified
:confval:`elevation` with the specified :confval:`albedo` and
:confval:`nighttime_solar_radiation_ratio` (these three parameters can
be GeoTIFF files instead of numbers). For some variables, the input
files are in different units than the default ones (hPa instead of kPa
for pressure, W/m² instead of MJ/m²/h for solar radiation) and need to
be converted (:confval:`unit_converter`).

If the :confval:`base_dir` contains ``tif`` files, the calculation is
performed once for each one of the sets of files; for example, if inside
:confval:`base_dir` there are files
:file:`temperature-2014-10-12-18-00+0200.tif`,
:file:`humidity-2014-10-12-18-00+0200.tif`, and so on (including
variables named ``wind_speed``, ``pressure``, and ``solar_radiation``),
there will be a resulting file
:file:`evaporation-2014-10-12-18-00+0200.tif`; if there are files for
other dates, there will be a result for them as well.  The calculation
is performed only if the resulting file does not already exist, or if at
least one of the input files has a later modification time.  If there
are any :file:`evaporation-....tif` files without corresponding input
files, they will be deleted.

If the :confval:`base_dir` contains ``hts`` files, the calculation is
performed for these time series. For example, if inside
:confval:`base_dir` there are files :file:`temperature.hts`,
:file:`humidity.hts`, and so on, there will be a resulting file
:file:`evaporation.hts`, overwriting any previously existing such file.

Configuration file reference
============================

The configuration file has the format of INI files, but without
sections.

Parameters
----------

.. confval:: loglevel

   Optional. Can have the values ``ERROR``, ``WARNING``, ``INFO``,
   ``DEBUG``.  The default is ``WARNING``.

.. confval:: logfile

   Optional. The full pathname of a log file. If unspecified, log
   messages will go to the standard error.

.. confval:: base_dir

   The directory in which ``vaporize`` will look for input files and
   write output files.  If unspecified, it is the directory from which
   ``vaporize`` was started.

.. confval:: step_length

   An integer indicating the number of minutes in
   the time step. In this version, ``vaporize`` can only handle hourly
   (60) or daily (1440) time steps.

.. confval:: elevation

   Meters of the location above sea level; this can be either a number
   or a GeoTIFF file with a digital elevation model.

.. confval:: nighttime_solar_radiation_ratio

   (Hourly step only.)

   In order to estimate the outgoing radiation, the ratio of incoming
   solar radiation to clear sky solar radiation is used as a
   representation of cloud cover. This, however, does not work during
   the night, in which case :confval:`nighttime_solar_radiation_ratio`
   is used as a rough approximation of that ratio. It should be a
   number between 0.4 and 0.8; see Allen et al. (1998), top of page
   75. It can be a number or a GeoTIFF file.

.. confval:: albedo

   A number between 0 and 1 or a GeoTIFF file with such numbers. It
   can also be a list of twelve space-separated numbers and/or GeoTIFF
   files, where the first is for January, the second for February, and
   so on. For example::

      albedo = albedo-jan.tif albedo-feb.tif albedo-mar.tif albedo-apr.tif
               albedo-may.tif albedo-jun.tif albedo-jul.tif albedo-aug.tif
               albedo-sep.tif 0.23           albedo-nov.tif albedo-dec.tif

   Note that in the configuration file long lines can be wrapped by
   indenting the additional lines. Also note that GeoTIFF files can be
   mixed with numbers; in the above example, GeoTIFF files are
   specified for all months except for October, which has a single
   value of 0.23.

   If a single number or GeoTIFF file is specified, it is used for all
   the year.

.. confval:: unit_converter

   The meteorological values that are supplied with the input files
   of the file set sections are supposed to be in the following units:

   ========================  =====================
   Parameter                 Unit
   ========================  =====================
   temperature               ℃
   humidity                  %
   wind speed                m/s
   pressure                  kPa
   solar radiation           MJ/m²/step
   sunshine duration         h
   ========================  =====================
   
   If they are in different units,
   :confval:`unit_converter_temperature`,
   :confval:`unit_converter_humidity`, and so on, are Python
   expressions that convert the given units to the above units; in
   these expressions, the symbol ``x`` refers to the given value. For
   example, if you have temperature in ℉, specify::
   
      unit_converter_temperature = (x - 32.0) * 5.0 / 9.0
      
   Use 32.0 rather than 32, and so on, in order to ensure that the
   calculations will be performed in floating point.

   You can also use this to convert wind speed to a different height.
   Wind speed at 2 m from the ground is required. If you have wind
   speed at a different height, convert it using Eq. 47, p. 56, of
   Allen et al. (1998). For example, if you have wind speed at 10 m,
   specify this:

      unit_converter_wind_speed = x * 4.87 / math.log(67.8 * 10 - 5.42)

.. confval:: temperature_prefix
             temperature_max_prefix
             temperature_min_prefix
             humidity_prefix
             humidity_max_prefix
             humidity_min_prefix
             wind_speed_prefix
             pressure_prefix
             solar_radiation_prefix
             sunshine_duration_prefix
             evaporation_prefix

   Optional. ``vaporize`` assumes that the input files are named
   :samp:`{variable}-{date}.tif` or :samp:`{variable}.hts`, where
   *variable* one of ``temperature``, ``temperature_max``,
   ``temperature_min``, ``humidity``, ``humidity_max``,
   ``humidity_min``, ``wind_speed``, ``pressure``, ``solar_radiation``,
   and ``sunshine_duration``, and, similarly, for the output file
   *variable* is ``evaporation``. With these parameters these names can
   be changed; for example::

      humidity_prefix = hum

   In that case, the humidity files are going to have a name similar to
   :file:`hum-2014-10-12-18-00+0200.tif` (for hourly) or
   :file:`hum-2014-10-12.tif` (for daily).

   ``vaporize`` will use the pressure if it is available in the input
   files, otherwise it will calculate it from the elevation.

References
==========

R. G. Allen, L. S. Pereira, D. Raes, and M. Smith, Crop evapotranspiration -
Guidelines for computing crop water requirements, FAO Irrigation and drainage
paper no. 56, 1998.
