# matcal

Radiocarbon (14C) age calibration using Bayesian statistics. Produces publication quality calibration plots, 1 sigma (68.27%) and 2 sigma (95.45%) calibrated age probabilities calculated using highest posterior density, as well as a probability density function of calibrated age for use in age modelling. Calibration output can be in either cal BP or BCE/CE (BC/AD), and a reservoir age can be specified if necessary. The user can choose from a number of calibration curves, including the latest version of IntCal. For Octave compatibility, see below.

For more detailed information and license information, see following manuscript:
Lougheed, B.C. and Obrochta, S.P., 2016. "MatCal: Open Source Bayesian 14C Age Calibration in MatLab." Journal of Open Research Software, 4: e42, DOI: http://dx.doi.org/10.5334/jors.130

How to install in Matlab:
-------------------------
(1) Create a directory called 'matcalfolder' somewhere on your computer and unzip the contents of the 'MatCal-master' repository to that directory. 

(2) For newer versions of Matlab (2012b and later): Go to the 'Home' tab, look under the 'Environment' section and click and click 'Set Path'. For older versions of Matlab (2012a and earlier): Go to the 'File' menu, select 'Set Path'.   
    
(3) Click 'Add Folder' to permanently add the 'matcalfolder' directory to your Matlab search path.

(4) Click 'OK' and then 'Save' to get back to the Matlab interface.

(5) Type 'help matcal' into the command window. If the install has succeeded, you will now see the matcal help documentation.

Octave compatibility:
---------------------
From version 2.2 (2017-02-20) onwards, Matcal also offers partial functionality in Octave, as tested using Octave 4.2.0 64-bit in Windows 7. Matrices containing HPD intervals and probability distribution function (PDF) will be calculated. Only very basic calibration plots are provided, due to the limitations of the Octave plotting capabilities.
