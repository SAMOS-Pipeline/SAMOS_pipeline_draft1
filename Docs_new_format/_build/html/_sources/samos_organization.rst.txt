=====================
Pipeline Organization
=====================

The SAMOS reduction pipeline is meant to be run within a ``Jupyter`` notebook,
as it allows for the user to track the reduction progress.  Begin each new
session with the execution of the ``SAMOSNight`` class.  This class stores
the data and information for a particular night (defined by some observation ID),
which is updated with each step of the data reduction.  The current capabilities
of the SAMOS pipeline are summarized in table, and each step is further described in the sections below.  Each section
is prefaced by the code used to execute the step.  The code represented in this
chapter is from the tutorial ``jupyter_NBtutorial/SAMOS_tutorial``.



.. _statustable:

.. list-table:: SAMOS Pipeline Status
   :header-rows: 1
   :widths: 10 20 10 20
   :stub-columns: 0

   *  -  Reduction Step
      -  Step Description
      -  Current Status
      -  Planned Updates
   *  -  `SAMOSNight`
      -  Initialize pipeline and organize raw data
      -  Complete
      -  Will eventually include the ability to make SQL queries for SAMOS data stored on a separate server.
   *  -  `ImageProcessor`
      -   Image trim, overscan, flat correction
      -   Complete
      -   Add bias level correction using bias frames
   *  -  `SlitBuckets`
      -  Trace and crop individual slits
      -  Complete
      -  Slit tracing method will need to read SAMOS slit mask patterns from the FITS.
   *  -  `WaveCalBuckets`
      -  Fits wavelength solution using comparison lamp data.
      -  Complete
      -  Linelists for first-guess solutions will refer to those obtained from calibration lamps observed with SAMOS.
   *  -  Sky Subtraction/Flux Calibration
      -  Subtract sky contribution from spectrum
      -  In progress
      -  (soon)


Major updates for this pipeline will occur when we obtain test data from the SAMOS instrument.
It should therefore be noted that some parts of the pipeline in its current state employ more
generic data reduction methods, with the plan that more sophisticated methods can be chosen
and applied once the actual parameters of the operational system are better known.  Until then,
we have adapted some working assumptions about the performance of the spectrograph which may
differ from its actual performance.  For instance, the pipeline assumes that the spectrum will
not perfectly co-align with the rows and columns of the detector, the data have an overscan
region which can be used for bias calibrations, and that other calibrations such as arc lamp
exposures can be taken in the afternoon or during the night.  The methods used for SAMOS data
anaysis are heavily influenced by the Goodman Spectroscopic Pipeline cite(GSP2019) and the
Astropy image reduction package ``ccdproc`` cite(CCDPROC).
