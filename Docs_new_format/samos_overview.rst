
=====================
Introduction to SAMOS
=====================



The SAMOS instrument is a multi-object spectrograph which uses a Digital Micromirror Device (DMD) [SMEE2018]_
to focus light towards either an imaging or spectroscopy channel.
The DMD allows for multiple slit configurations in a single observing run.
Multi-object capabilities of the *Goodman High Throughput Spectrograph* require slit masks to be made in advance,
and are installed in the afternoon leading up to the observation.
SAMOS has the unique feature of being able to create on-the-fly slit patterns, which will be saved in the FITS headers.
The instrument is not set to be commissioned until at least 2021, so the current version of the SAMOS data reduction pipeline (SRP)
uses multi-object data from the SOAR Goodman Spectrograph.
This SRP will provide the foundation for SAMOS data reduction upon its commissioning.  The test data used for the SRP was taken by \textit{SOAR Goodman} on March 19, 2014.  Eventually, the pipeline will be used to reduce multi-object spectroscopy data taken with the SOAR Adaptive-Module Optical Spectrograph (SAMOS).


There are two main parts to the SRP: basic image reduction and spectroscopic reduction.
For optical spectroscopy, we use detectors called charge-coupled devices (CCDs) to detect and convert light to a digital signal.
When sensors on the CCD are hit by incoming photons, they save the information and it gets read out to a two-dimensional pixel
array after the exposure.  The exposure is saved as a FITS (Flexible Image Transport System) file, which is able to store both the
image data and calibration information in a header.


There are a couple important systematics at play during exposure and readout, an obvious one being that CCDs are not perfect photon detectors.
This imperfection means that pixels can vary in sensitivity across the array.  We correct for this variation by normalizing, or flattening,
the data.  Additionally, CCDs are not completely cleared of their information after readout,
which means that each exposure starts with a certain amount of signal.  We call this "pre-charge" is called the readout bias.
Sometimes a region of the CCD called the "overscan" is used to correct for this readout noise.
A cartoon layout of a CCD is shown in figure \ref{fig:spec_ex}.


.. figure:: im_spec_abell.pdf
  :width: 600px
  :align: center
  :height: 500px
  :alt: alternate text
  :figclass: align-center

  Simplified depiction of multi-object spectra.  Sky light is also dispersed and appears as emission lines on either side of the spectrum.
  The lines are the same in each spectrum, but shift depending on slit position.



When dispersed photons [#]_ hit the detector, they are mapped to the 2-D array, the mapping  of light intensity,
location along the spectroscopic slit (1 dimension), and wavelength (orthogonal dimension to slit location)
is not perfectly aligned with the pixels.  We usually correct for this distortion by fitting the continuum (track of illuminated pixels)
of a standard star as a function of pixel location \citep{Marsh_1989,Horne_1986}.  This spectrum trace,
shapes measured at different slit positions, can be used to rectify the corresponding science images.
Once the correction is applied to the 2-D spectrum, we collapse it to one dimension for wavelength calibration.
Wavelength calibration involves converting pixel coordinates of the 1-D spectrum to wavelength coordinates based on
characteristics of the spectrum grating.

Because each instrument comes with its own systematics, image reduction is rarely a cookie-cutting exercise,
and spectral tracing and wavelength calibration are particularly difficult and tedious.
Thus far, the Image Reduction and Analysis Facility (IRAF) [IRAF1986]_, [IRAF1993]_ has been the only reduction package with
enough robustness and versatility to make it the go-to spectroscopic data reduction source for nearly 30 years.
However, the command language (CL) in which IRAF is written has been unable to keep up with more versatile languages such as Python.
IRAF's imminent demise is a source of worry, and there have only been a few efforts made to start the transition to a new generic
reduction package [astropy]_.
The goal of the SRP is to have a complete spectral reduction package for SAMOS data which is completely independent of IRAF.


.. [SMEE2018] https://ui.adsabs.harvard.edu/abs/2018SPIE10546E..0LS/abstract
.. [astropy] https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract
.. [IRAF1986] https://ui.adsabs.harvard.edu/abs/1986SPIE..627..733T/abstract
.. [IRAF1993] https://ui.adsabs.harvard.edu/abs/1993ASPC...52..173T/abstract

.. [#] Like a rainbow.
