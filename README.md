Malawi PSHA using the MSSD

Last Edited by Jack Williams 14/03/22 (jack.williams@otago.ac.nz)

The following describes the necessary files and steps to (1) generate the Malawi Seismogenic Source Database (MSSD; Williams et al 2021b), and (2) use these data to perform a PSHA for Malawi (Williams et al 2022).

Codes are written in MATLAB and require the following toolboxes:
# Image Processing
# Mapping Toolbox


Input files and folders:

#### THE FOLLOWING INPUT FILES ARE STORED ON ZENODO AND MUST BE DOWNLOADED AND STORED IN PATH AS SHOWN PRIOR TO RUNNING CODES ####

misc_functions/malawi_Vs30_active.txt: USGS slope-based Vs30 values for Malawi (Wald and Allen 2007)

syncat_MSSD/EQCAT_comb.mat: Matab file for MSSD Direct catalog for all possible rupture weightings

psha_MAP_em/GM_MSSD_em_20220302: Ground motions needed to plot PSHA maps

psha_site_em/GM_MSSD_em_20220302.mat: Ground motions needed to plot PSHA-site analysis figures

mssd_comb.mat: Matlab file for combined MSSD Direct and Adapted MSSD catalogs

syncat_adaptedMSSD/MSSD_Catalog_Adapted_em.mat: Matlab file for Adapated MSSD event catalog

syncat_bg/syncat_bg.mat: Areal source stochastic event catalog

####### The following files are included with this release on Github ######

syncat_PSHA_MSSD_input.xls: Main spreadsheet for editing MSSD and PSHA inputs. Any updates must be saved to the Matlab variable ‘syncat_PSHA_MMSD_input.’

MSSD_sourcecalcs/MSSD.xls: Spreadsheet with data for each source in the Malawi Seismogenic Source Database. Contains 7 tabs:

#1 ‘Leonard 2010’: list of sections, faults, and multi-faults where geometric data is stored to perform slip rate, earthquake magnitude, and recurrence interval calculations using the Leonard 2010 scaling relationships as described in Williams et al (2021a). Uncertainties from extreme ends of logic tree

#2-4 ‘Section Geometry,’ ‘Fault Geometry,’ ‘MultiFaultGeometry’:  list of sections with information on the geometry and source properties of each rupture type identified in the MSSD.

#5‘MSSD_AdaptedSources’: list of information necessary for PSHA using Youngs & Coppersmith 1985 approach. Each ‘section,’ ‘fault,’ or ‘multifault’ is contained within the largest system it is part of (i.e. fault or multi fault.’

#6‘BasinSpecificValues’: The plate motion vector (and 1 sigma error) for each basin in Malawi, calculated from the geodetic model described in Wedmore et al (2021). Used in slip rate calculations

#7‘FaultObliquity’: Spreadsheet to compare source moment rates in catalogs and their 'analytical' moment rate. Can also calculate their moment rate if they were optimally oriented to the regional extension azimuth and extended through the 35 km seismogenic layer in Malawi (see appendix A3 in Williams et al 2022).


The spreadsheets ‘Section Geometry,’ ‘Fault Geometry,’ ‘MultiFaultGeometry,’ and ‘MSSD_ApdatedSources’ are linked to ‘Leonard2010’ through VLOOKUP functions, and linked to Matlab through the generation of the ‘MSSD_sources’ .mat file using the ‘Generate_MATfile.m’ script

SZ_final.m: Source information for areal PSHA sources in East Africa as developed by the Global Earthquake Model (GEM; Poggi et al 2017). Additional sources for areas not covered by GEM are also included, where seismicity rates are defined by Fenton et al (2006).

Malawi Vs30_active.txt: vs30 values for Malawi as derived from topographic slope analysis (Wald and Allen 2007)

gis_files: Various GIS files to aid plotting maps such as national borders and Lake Malawi. Also includes the Malawi Active Fault Database (MAFD; Williams et al 2021c&d) and Malawi Seismogenic Source Database (MSSD) GIS files (Williams et al 2021b&e).

misc_functions: Various Matlab functions to support scripts as described below

MSSD_sourcecalcs/HangingWallFlexureInputs: Parameters for hanging-wall flexural calculations for intrarift sources in Lake Malawi basins (see Williams et al 2021b). Includes cross check with values for the Kenya Rift from Muirhead et al (2016).

MSSD_sourcecalcs/nbasin_fautdata/nb_slipratetable_id.csv: Slip rate values from offset seismic reflector for intrarrift faults in Lake Malawi (Shillington et al 2020).

GMPE/GMPE_Malawi: Function to apply GMPE to events in stochastic event catalogs during PSHA

GMPE/GMPEcoef_Malawi: Coefficients for GMPE used in PSHA

ISC_GEM_EAcatalog.mat: ISC catalog since for East Africa Seismicity (Poggi et al 2017).


#### To perform PSHA, run scripts in this order ####

1.) MSSD Source Geometry

MSSD_sourcegeom/fault_geom_1.m: Script for constructing 3D geometrical models of sources, where they are defined by 2D planes and/or grid points with an average spacing of 1 km. Points created using polygrid function (Sattari 2022) and p_poly_dist and p_poly_dist1 functions (Yoshpe 2022). The width of each source is defined by which ever value is lower out of: (1) its width as defined by an empirical scaling relation with source length (Leonard 2010) or (2) if it is intersected and cut by another source. Analysis of the latter performed using the function ‘fault intersect,’ which requires Surface Intersection function (Tuszynski 2022). The revised area of a source that is intersected is input into the MSSD.xls spreadsheet to revise source geometry. Returns variable ‘fault_geom_1,’ where grid points are sorted by fault or multifault, ‘sec_geom_1’ where grid points are sorted by section, 'and ‘fault_geometry_points_1’ where sources are defined by planes. Planes also stored as csv file MSSD_source_geometry'

MSSD_sourcegeom/fault_geom_2.m: Same as for fault_geom_1, but width of source is defined as the lower value out of: (1) extrapolation to the base of the seismogenic crust (35 km) (2) if it is intersected and cut by another source. Returns variable ‘fault_geom_2,’ where grid points are sorted by fault or multifault, and  ‘fault_geometry_points_2’ where sources are defined by planes.

MSSD_sourcegeom/Plot_MSSDpolygons.m: Option to plot basic 3D geometric model of MSSD sources using 'fault_geom_1'

MSSD_sourcegeom/source_map_figure.m: Option to plot 3D geometrical model of MSSD sources with section sources shown, SRTM-30 m Digital Elevation Model for Malawi (Sandwell 2011), Lake Malawi, and international boundaries. Requires Hillshade function (Hebeler 2022)

MSSD_sourcegeom/fault_intersect_figure.m: Option to plot 3D geometrical model of two faults intersecting

MSSD_sourcegeom/floating_rupture_figure.m: Option to plot 3D geometrical model of the distribution of a source’s events in the stochastic event catalogs once they have been randomly floated on the source. Figure plotted for source defined by 'source_id' in code. Requires source2site functions.

2.) MSSD Source Calcs

MSSD_sourcecalcs/MSSD_sourcecalcs.m: Takes fault geometry, geodetic extension rates, hanging-wall flexure calculations, and Leonard (2010) scaling relationships, and uses them to calculate fault slip rates and recurrence intervals through 10000 Monte Carlo simulations through MSSD logic tree (Williams et al 2021b). Results in terms of mean and 1 standard deviation can be saved as ‘MSSD sources_Calc’ and written into MSSD spreadsheet. For faults with slip rates estimated from the offset seismic reflector, slip rates are initially estimated from the systems-based approach, preserved in variable slip_rate_nb.sr, and then replaced with offset reflector values.

MSSD_sourcecalcs/MSSD_sourcecalcplots.m: Plots various figures to analyse MSSD_souces_calcs, such as prob density of fault slip rates, comparison between slip rates from geodesy and offset seismic reflector, histograms, and comparisons to Zomba Graben sources with the SMSSD (SMSSD.xlsx; Williams et al 2021a).

MSSD_sourcecalcs/MSSD_sourcecalcplots_nbasincomparison.m: Compares the slip rate probability distributions of intrarift faults in Lake Malawi when estimated from the offset reflector and when estimated from the systems-based approach. Requires: calc_overlap_twonormal (Kowerko 2022)

Generate_MATfile.m: With MSSD source geometry and calculations complete, store data in variable MSSD_sources

!Any changes to MSSD files should be saved to MSSD.xlsx and relevant .csv files so it can be updated for the MSSD shape files!

3.) Stochastic Event Catalogs

(a) Areal source Catalogs

syncat_bg/new_SZ_GR_Relation.xlsx: Derive a- and b-values of the G-R relationship for new areal source zones in areas not covered by Poggi et al (2017) based on their area and G-R parameters of stable carton seismicity from Fenton and Boomer (2006).

syncat_bg/new_SZ.m: Using a- and b- values from new_SZ_GR_Relation.xlsx, add new areal source zones to SZ_final.mat (Poggi et al 2017) using areal extents defined by 'GEM_SourceZones_WGS.shp'. Save in new Matlab structure 'new_SZ'

syncat_bg/syncat_bg.m: Script for simulating stochastic catalog of earthquakes as defined by the areal sources in the Matlab structure new_SZ. Returns variable ‘syncat_bg.’ Earthquakes treated as point sources, and considered to represent background/off fault seismicity in PSHA. Events >200 km from Malawi are filtered out using p_poly_distance function (Yoshpe 2022).

syncat_bg/syncat_bg_tmp.mat: 1000 year long bg event catalog for showing event distributions

syncat_bg/syncat_bg_analysis.m: Optional script to assess simulated syncat_bg: (1) find moment rate of all sources within Malawi, (2) plot mag-freq diagram for all sources and compared to theoretical rate, (3) plot distribution of events in a 1000 year long record (using syncat_bg_tmp)

syncat_bg/inMalawi-GR.xlsx: Spreadsheet for a-values of areal sources, scaled for their overlap area in Malawi

syncat_bg/GR_EastAfrica_GEM.m: Option to plot ISC catalog for recorded seismicity in Malawi (Poggi et al 2017), and saves to variable ‘ISC_GEM_EAcatalog_Malawi’. Also to plot 95% confidence intervals for Gutenberg-Richter (G-R) relationship in Malawi using approach described in Tinti & Malaria (1987) and function 'GRrelation_MLEWeichert_EQMATca'. This can be saved to variable 'GR_TintiMulargia' for plotting in MSSD_comb.


(b) MSSD Direct

syncat_MSSD/rupture_weighting.m: Script for simulating stochastic catalog of earthquakes as defined by a direct interpretation of the earthquake magnitudes and recurrence intervals in the MSSD and Poission interevent times. Runs for all source combinations of 'section', 'fault', and 'multi fault sources' to determine best-fit w combination to regional b-value. This combination should be updated in syncat_PSHA_MSSD_input.xlsx. Outputs can be saved in 'EQCAT_comb.mat' and'r_weighting_results.csv' and results plotted in figures.

syncat_MSSD/syncat_MSSD.m: Runs as above for selected source combination. Returns variable ‘eqcat_test’ where the earthquake catalog is stored under ‘EQCAT’ and to be used for future PSHA

syncat_MSSD/syncat_MSSD_analysis.m: Optional. For plotting results of EQCAT such as magnitude-frequency distribution of record and comparison of sources moment rates. Also calculates and compares analytical and event-catalog moment rate of each source and plots. Saves in variable 'syncat_comparsion.mat' which can be used on file 'MSSD_comb.'

(c) MSSD Adapted

syncat_MSSD_adapted/syncat_AdaptedMSSD_em.m:  Script for simulating stochastic catalog of earthquakes as defined by slip rate and area of sources in ‘MSSD_AdaptedSources,’ and the methods described in Youngs and Coppersmith (1985) and Goda and Sharipov, A.(2021). Uses function ‘characteristic_magnitude_YC1985’ and returns variable ‘MSSD_Catalog_Adapated_em’ where the event catalog is stored under ‘StochasticEventCatalog.’ G-R or char behaviour or length or crust limited fault widths are considered separately in distinct catalogs 

syncat_MSSD_adapted/syncat_AdaptedMSSD_analysis_em.m:  Optional. As above but for plotting analysis of StochasticEventCatalog_em. Includes plotting theoretical and catalog Y&C85 mfd curves for a single source.

(d) MSSD_comb

mssd_comb: Optional. For combining Direct MSSD, all AdaptedMSSD_em, and bg catalog together into one 10 million year long catalog for analysis. Includes plotting MFD of all catalogs and checking respective moment rates. Also for sampling MSSD-combined in 50 year intervals and and comparing it to the SSA-GEM catalog (using syncat_bg.mat and GR_TintiMulargia.mat created in syncat_bg folder) and geodetic models (as stored in misc_functions/geo_mo_rate).


4.) GMPE

This folder stores GMPE functions that are assessed during the PSHA. Note, these codes DO NOT need to be run before PSHA. GMPE selection is made in syncat_PSHA_MSSD_input.xlsx

GMPE/GMPE_Malawi: Function to run events simulated in event catalogs through Ground Motion Models (GMM) or Ground Motion Prediction Equations (GMPE) selected in 'syncat_PSHA_MSSD_input.' Uses coefficients in MATLAB function GMPEcoef_Malawi.m and variable GMPEcoef_Malawi.m. Output is median ground motion intensity and log standard deviation.

GMPE/GMPE_figure: Create figure to compare selected GMPE for a single hypothetical event 

GMPE/seismicVelocityComparisonPlot: Plot 1D seismic velocity models for Malawi from Ebinger et al 2019 and Stevens et al 2021. Also compares to Boore et al 2016 curve for active stable crust. Uses data stored in vp.csv and MalawiSeismicVelocityComparison.xlsx

5.) PSHA

***NOTE*** Require evaluating background catalog and all 5 MSSD catalogs through 4 GMPE. As a result, analysis is performed by equally dividing the catalog into 16 sub catalogs, and then evaluating these subcatalogs simultaneously on a cluster. These codes are currently set up for use on the Blue Crystal cluster at the University of Bristol (http://www.bris.ac.uk/acrc/). These codes may therefore require adaption for use elsewhere.

#### Each script is identical apart from the ‘num_par’ parameter which specifies which catalog interval it assesses. Use ‘par_opts’ variable to change how catalog is assessed ####

5a) PSHA Site

psha_site_em/PSHA_MSSD_sitexx.m: Site specific PSHA for Malawi. Returns ground motions at site for interval of stochastic event catalogs as denoted by xx. Each event in catalog is assessed using all ground motion prediction equations selected, vs30 =300 m/s or USGS value, and vs30=760 m/s and run using functions Source2Site, Source2Site3, GMPE Malawi. Runs first for background catalog, then for MSSD-based catalogs. Ground motions for each catalog are assessed over multiple spectral accelerations, and multiple sites can be assessed. Only 'num_need' ground motions are stored, with num_need determined by the length of the catalog and minimum return period assessed (currently set to 100 years, can be altered in syncat_PSHA_MSSD_input.xlsx). The num_need ground motions that are stored are iteratively updated as the catalog is assessed to stop file sizes becoming too large. They are then stored as 'AMAX_GM_xx' and AMAX_GM_bg_xx' for further analysis. In addition, the magnitude, GMPE, and source2site distances for the num_need largest ground motions are stored for use during disaggregation analysis.

psha_site_em/PSHA_MSSD_site_figures.m: Combines ground motions from all parallelisations, and then samples only num_need. Uses data for site-specific PSHA plots: (1) seismic hazard curves, (2) uniform spectra, (3) disagg plots. Since an ensemble modelling approach used, seismic hazard curves can be plotted for all event-catalogs-GMPE selection, and for a given PoE (currently 2% and 10% PoE in 50 years) the distribution of ground motions plotted using a beta distribution). *Option to plot data from previously run PSHA codes by loading GM_20220221* Otherwise new results from subcatalogs should be combined by uncommenting relevant text. Dissagg plots require functions disagg_MSSD_JW and voxel (Joel, 2022).

5b) PSHA Maps

psha_map/PSHAmap_MSSD_xx.m: PSHA code for generating Malawi PSHA Maps. Analysis only performed for one spectral acceleration (currently set as PGA). Maps created by performing site specific for multiple sites with grid spacing as defined in syncat_PSHA_MSSD_input.xls, USGS vs30 value and 760 m/s, and for PGA only. Code run over background and MSSD-based catalog interval as defined by xx and considers each GMPE. As with site-specific PSHA, only num_need ground motions are stored to prevent excessive file sizes. 

psha_map/PSHAmap_MSSD_figures.m: PSHA Maps for Malawi using results from PSHAmap_MSSD_xx.m. Maps plotted in terms of contribution from aerial (background) and MSSD sources. Also option too plot comparison for different on-fault magnitude-frequency distribution comparison (by comparing G-R to MSSD-Direct and characteristic), different fault widths, and PSHA uncertainty maps (in terms of CoV and interquartile range of 20 ground motion values calculated for each site) *Option to plot data from previously run PSHA codes by loading GM_20220302* Otherwise new results from subcatalogs should be combined by uncommenting relevant text. Also option to saves ground motions in txt files so can be compared to Poggi et al (2017) maps and Hodge et al (2015) PSHA maps. Currently set for maps 10% PoE in 50 years and 2% PoE in 50 years.

psha_map/psha_map_comparison/hodge_2015_psha_comparison.m: Optional. Script to recreate Hodge 2015 PSHA results for Malawi, and then compare ground motions with MSSD results. Requires h2015_map_values and shape file of Hodge et al (2015) sources as input. Currently set for 2% PoE in 50 years only. Saves results so can combine with GEM comparison in gem_psha_comparison. Note minimum return period is 1 in 500 year (I.e. less than 10% PoE in 50 years). To compare with MSSD files, PSHAmap_MSSD_figures.m must be run first and output seismic hazard values stored in txt files

psha_map/psha_map_comparison/gem_psha_comparison.m: Optional. Compares PSHA maps to those produced by GEM (Poggi et al 2017) and Hodge et al (2015). Data from gem stored in folder in psha_map/gem_gm.mat. Note this data is for 2% and 10% PoE only. For option to combine with figure with Hodge 2015 maps comparison, hodge_2015_psha_comparison must be run first. Also provides some quantitative comparisons between GEM and MSSD PSHA maps.


References (Literature)

Allen, T. I., & Wald, D. J. (2009). On the use of high-resolution topographic data as a proxy for seismic site conditions (VS 30). Bulletin of the Seismological Society of America, 99(2A), 935-943.

Ebinger, C. J., Oliva, S. J., Pham, T. Q., Peterson, K., Chindandali, P., Illsley‐Kemp, F., ... & Mulibo, G. (2019). Kinematics of active deformation in the Malawi rift and Rungwe Volcanic Province, Africa. Geochemistry, Geophysics, Geosystems, 20(8), 3928-3951.

Fenton, C. H., Adams, J., & Halchuk, S. (2006). Seismic hazards assessment for radioactive waste disposal sites in regions of low seismic activity. Geotechnical & Geological Engineering, 24(3), 579-592.

Goda, K., & Sharipov, A. (2021). Fault-Source-Based Probabilistic Seismic Hazard and Risk Analysis for Victoria, British Columbia, Canada: A Case of the Leech River Valley Fault and Devil’s Mountain Fault System. Sustainability, 13(3), 1440.

Hodge, M., Biggs, J., Goda, K., & Aspinall, W. (2015). Assessing infrequent large earthquakes using geomorphology and geodesy: the Malawi Rift. Natural Hazards, 76(3), 1781-1806.

Leonard, M. (2010). Earthquake fault scaling: Self-consistent relating of rupture length, width, average displacement, and moment release. Bulletin of the Seismological Society of America, 100(5A), 1971-1988.

Muirhead, J. D., Kattenhorn, S. A., Lee, H., Mana, S., Turrin, B. D., Fischer, T. P., ... & Stamps, D. S. (2016). Evolution of upper crustal faulting assisted by magmatic volatile release during early-stage continental rift development in the East African Rift. Geosphere, 12(6), 1670-1700.

Poggi, V., Durrheim, R., Tuluka, G. M., Weatherill, G., Gee, R., Pagani, M., ... & Delvaux, D. (2017). Assessing seismic hazard of the East African Rift: a pilot study from GEM and AfricaArray. Bulletin of Earthquake Engineering, 15(11), 4499-4529.

Sandwell, D., Mellors, R., Tong, X., Wei, M., & Wessel, P. (2011). Open radar interferometry software for mapping surface deformation.

Shillington, D. J., Scholz, C. A., Chindandali, P. R., Gaherty, J. B., Accardo, N. J., Onyango, E., ... & Nyblade, A. A. (2020). Controls on Rift Faulting in the North Basin of the Malawi (Nyasa) Rift, East Africa. Tectonics, 39(3), e2019TC005633.

Stevens, V. L., Sloan, R. A., Chindandali, P. R., Wedmore, L. N., Salomon, G. W., & Muir, R. A. (2021). The Entire Crust can be Seismogenic: Evidence from Southern Malawi. Tectonics, e2020TC006654.

Tinti, S., & Mulargia, F. (1987). Confidence intervals of b values for grouped magnitudes. Bulletin of the Seismological Society of America, 77(6), 2125-2134.

Wedmore, L. N., Biggs, J., Floyd, M., Fagereng, Å., Mdala, H., Chindandali, P., ... & Mphepo, F. (2021). Geodetic constraints on cratonic microplates and broad strain during rifting of thick Southern African lithosphere. Geophysical Research Letters, 48(17), e2021GL093785.

Williams, J. N., Mdala, H., Fagereng, Å., Wedmore, L. N., Biggs, J., Dulanya, Z., ... & Mphepo, F. (2021a). A systems-based approach to parameterise seismic hazard in regions with little historical or instrumental seismicity: active fault and seismogenic source databases for southern Malawi. Solid Earth, 12(1), 187-217.

Williams, J. N., Wedmore, L. N., Fagereng, Å., Werner, M. J., Mdala, H., Shillington, D. J., ... & Chindandali, P. (2021b). Geologic and geodetic constraints on the seismic hazard of Malawi’s active faults: The Malawi Seismogenic Source Database (MSSD). Natural Hazards and Earth System Sciences Discussions, 1-47.

Williams, J. N., Wedmore, L. N., Scholz, C. A., Kolawole, F., Wright, L. J., Shillington, D., ... & Werner, M. J. (2021c). The Malawi Active Fault Database: an onshore-offshore database for regional assessment of seismic hazard and tectonic evolution. Earth and Space Science Open Archive ESSOAr.

Williams, Jack, Wedmore, Luke, Scholz, Christopher A, Kolawole, Folarin, Wright, Lachlan J M, Shillington, Donna J, Fagereng, Å, Biggs, Juliet, Mdala, Hassan, Dulanya, Zuze, Mphepo, Felix, Chindandali, Patrick, & Werner, Maximilian J. (2021d). Malawi Active Fault Database (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5507190

Williams, Jack N., Wedmore, Luke N. J., Fagereng, Åke, Werner, Maximilian J., Biggs, Juliet, Mdala, Hassan, Kolawole, Folarin, Shillington, Donna J., Dulanya, Zuze, Mphepo, Felix, Chindandali, Patrick R. N., Wright, Lachlan J. M., & Scholz, Christopher A. (2021e). Malawi Seismogenic Source Database (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5599617

Williams, J. N., Werner, M. J., Goda. K., Wedmore, L. N., De Risi R., Biggs, J., Mdala, H., Dulanya, Z., Fagereng, Å., Chindandali, P.,  Mphepo, F. (2022) Fault-based probabilistic seismic hazard analysis in regions with low strain rates and a thick seismogenic layer: a case study from Malawi. Submitted to Natural Hazards

Youngs, R. R., & Coppersmith, K. J. (1985). Implications of fault slip rates and earthquake recurrence models to probabilistic seismic hazard estimates. Bulletin of the Seismological society of America, 75(4), 939-964.

References (Matlab Functions)

Crameri, F. (2018). Scientific colour-maps. Zenodo. <http://doi.org/10.5281/zenodo.1243862>

Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018.

Felix Hebeler (2022). Hillshade (https://www.mathworks.com/matlabcentral/fileexchange/14863-hillshade), MATLAB Central File Exchange. Retrieved January 9, 2022.

Suresh Joel (2022). Voxel (https://www.mathworks.com/matlabcentral/fileexchange/3280-voxel), MATLAB Central File Exchange. Retrieved January 10, 2022.

Danny Kowerko (2022). calc_overlap_twonormal(s1,s2,mu1,mu2,xstart,xend,xinterval) (https://www.mathworks.com/matlabcentral/fileexchange/49823-calc_overlap_twonormal-s1-s2-mu1-mu2-xstart-xend-xinterval), MATLAB Central File Exchange. Retrieved January 9, 2022.

Rafael Palacios (2022). deg2utm (https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm), MATLAB Central File Exchange. Retrieved January 9, 2022.

Sulimon Sattari (2022). Grid of points within a polygon (https://www.mathworks.com/matlabcentral/fileexchange/41454-grid-of-points-within-a-polygon), MATLAB Central File Exchange. Retrieved January 10, 2022.

Jaroslaw Tuszynski (2022). Surface Intersection (https://www.mathworks.com/matlabcentral/fileexchange/48613-surface-intersection), MATLAB Central File Exchange. Retrieved January 9, 2022.

Michael Yoshpe (2022). Distance from points to polyline or polygon (https://www.mathworks.com/matlabcentral/fileexchange/12744-distance-from-points-to-polyline-or-polygon), MATLAB Central File Exchange. Retrieved January 10, 2022.
