Auxiliary Material Submission for Paper 2009JB006325

Incorporating and reporting uncertainties in fault slip rates

J. Douglas Zechar
Lamont-Doherty Earth Observatory, Columbia University, Palisades, New York, USA

Kurt L. Frankel
School of Earth and Atmospheric Sciences, Georgia Institute of Technology, Atlanta, Georgia, USA


Zechar, J. D., and K. L. Frankel (2009), Incorporating and reporting uncertainties in fault 
slip rates, J. Geophys. Res., 114, B12407, doi:10.1029/2009JB006325.

Introduction
In this auxiliary material, we provide MATLAB codes used to produce the results in our 
paper, and we also provide example scripts and the corresponding outputs.  In general, 
the codes we described in our paper and that we have provided here allow you to perform 
three computations:

* compute an age p.d.f. for a given set of ages (using age*.m)
* compute a displacement p.d.f. for a given set of displacements (using displacment*.m)
* compute a slip rate p.d.f. for a given age p.d.f. and a given displacement p.d.f. (using 
slip_rate.m)

The typical workflow is to compute the age and offset p.d.f.s using your data and then 
use these results to compute a slip rate p.d.f.  As we describe in our paper, there are a few 
choices for the age and displacement uncertainty model, and the corresponding scripts 
have slightly different parameters, all of which are documented within the scripts 
themselves.  For help on any of the scripts, you can type help script_name into the 
MATLAB command window.  For example, if you want to know what parameters are 
required for the boxcar age uncertainty script, type

help age_boxcar

and a description of the function and the relevant input and output parameters will 
appear.

We have also included several examples corresponding to the cases we consider in the 
paper, and these follow the workflow we mentioned above:

example_fc: Furnace Creek example, considering two cases: all ages reported in Frankel 
et al. 2007 GRL and the "good" ages used for their analysis.  This example, and all others 
includes computation of 68.27% confidence region and 95.45% confidence region for 
age, offset, and slip rate.
example_ic: Indian Creek example
example_kaneda: S. Tenjindo example, considering two cases: treatment of reported 
displacement as boxcar, and treatment of reported displacement as Gaussian
example_rwc: Red Wall Canyon example, considering two cases as in the Furnace Creek 
example

These are the scripts that were used to generate the figures in our paper.  All of the 
functions are fairly well-commented.  If something is unclear or you discover a bug, 
please let us know by emailing zechar@ldeo.columbia.edu.

There are four auxiliary data sets, described below.

1. 2009jb006325-ds01.txt Probability distributions for offset, age, and slip rate at Furnace 
Creek site.
1.1 Column "offset", m, offset value of interest
1.2 Column "offset_pdf", probability density at the specified offset value
1.3 Column "age_all", a, age value of interest, considering all samples
1.4 Column "all_ages_pdf", probability density at the specified age_all value
1.5 Column "age_excluding_outliers", a, age value of interest, excluding outlier ages
1.6 Column "excluding_outliers_age_pdf", probability density at the specified 
age_excluding_outliers value
1.7 Column "slip_rate_all_ages", mm/a, slip rate value of interest, considering all ages
1.8 Column "all_ages_slip_rate_pdf", probability density at the specified slip rate value, 
considering all ages
1.9 Column "slip_rate_excluding_outlier_ages", mm/a, slip rate value of interest, 
excluding outlier ages
1.10 Column "excluding_outlier_ages_slip_rate_pdf", probability density at the specified 
slip rate value, excluding outlier ages

2. 2009jb006325-ds02.txt Probability distributions for offset, age, and slip rate at Indian 
Creek site.
2.1 Column "offset", m, offset value of interest
2.2 Column "offset_pdf", probability density at the specified offset value
2.3 Column "age", a, age value of interest
2.4 Column "age_pdf", probability density at the specified age value
2.5 Column "slip_rate", mm/a, slip rate value of interest
2.6 Column "slip_rate_pdf", probability density at the specified slip rate value

3. 2009jb006325-ds03.txt Probability distributions for offset, age, and slip rate for 
Neodani fault.
3.1 Column "offset_boxcar", m, offset value of interest
3.2 Column "boxcar_offset_pdf", probability density at the specified offset_boxcar value, 
assuming boxcar uncertainty
3.3 Column "offset_Gaussian", m, offset value of interest
3.4 Column "Gaussian_offset_pdf", probability density at the specified offset_Gaussian 
value, assuming Gaussian uncertainty
3.5 Column "age", a, age value of interest
3.6 Column "age_pdf", probability density at the specified age
3.7 Column "slip_rate_boxcar_offset", mm/a, slip rate of interest
3.8 Column "slip_rate_pdf_boxcar_offset", probability density at the specified slip rate 
value, assuming boxcar uncertainty
3.9 Column "slip_rate_Gaussian_offset", mm/a, slip rate of interest
3.10 Column "slip_rate_pdf_Gaussian_offset", probability density at the specified slip 
rate value, assuming Gaussian uncertainty

4. 2009jb006325-ds04.txt Probability distributions for offset, age, and slip rate at Red 
Wall Canyon site.
4.1 Column "offset", m, offset value of interest
4.2 Column "offset_pdf", probability density at the specified offset value
4.3 Column "age_all", a, age value of interest, considering all samples
4.4 Column "all_ages_pdf", probability density at the specified age_all value
4.5 Column "age_excluding_outliers", a, age value of interest, excluding outlier ages
4.6 Column "excluding_outliers_age_pdf", probability density at the specified 
age_excluding_outliers value
4.7 Column "slip_rate_all_ages", mm/a, slip rate value of interest, considering all ages
4.8 Column "all_ages_slip_rate_pdf", probability density at the specified slip rate value, 
considering all ages
4.9 Column "slip_rate_excluding_outlier_ages", mm/a, slip rate value of interest, 
excluding outlier ages
4.10 Column "excluding_outlier_ages_slip_rate_pdf", probability density at the specified 
slip rate value, excluding outlier ages
