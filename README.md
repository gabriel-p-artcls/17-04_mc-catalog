

- [mc-catalog](#mc-catalog)
  - [Top level](#top-level)
    - [`ages_mass_lit/`](#ages_mass_lit)
    - [`aux_funcs/`](#aux_funcs)
    - [`extinction_MCEV/`](#extinction_mcev)
    - [`figures/`](#figures)
    - [`functions/`](#functions)
    - [`runs/`](#runs)


# mc-catalog

Set of scripts to analyze and plot the results of the **ASteCA** processing of
the 210 MC clusters in our database, observed with `C, T1` Washington
photometry.

## Top level

* `README.md`

 This file.

* `lista_unica_cumulos.ods`

 Literature data on each cluster.

* `bb_cat.dat`

 RA & DEC positions for the 3740 clusters in the
 [Bica et al. (2008)](http://cdsads.u-strasbg.fr/abs/2008MNRAS.389..678B)
 catalog ([Table 3](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/
 MNRAS/389/678/table3)).

* `mc_cat_analysis.py`

 Script that Produces the main figures by comparing parameters with each other
 (KDE maps), with literature values (1:1 plots) and distributed in space
 (RA vs DEC plots.)

 Functions are stored in the `functions/` folder.

* `asteca_output_final.dat`

 Combined final results for the entire sample.


### `ages_mass_lit/`

* `BB_10.dat`

 List of ages for 151 SMC and 539 LMC clusters in the database presented in
 [Bonato & Bica (2010)](http://cdsads.u-strasbg.fr/abs/2010MNRAS.403..996B).
 Only 285 out of the 690 clusters have an age value assigned.

* `chiosi_06.dat`

 Contains 311 SMC clusters younger than 1 Gyr (ages determined using isochrone
 fitting) from the work by [Chiosi et al. (2006)](http://adsabs.harvard.edu/
 abs/2006A%26A...452..179C). Also lists several other structures/objects, ie:
 association, supernovas, H II regions and entities in between.

* `cross_match.py`

 Matches the data in the above databases to those clusters processed by
 ASteCA.

* `glatt_10.dat`

 CMD ages taken from the [Glatt et al. (2010)](http://www.aanda.org/10.1051/
 0004-6361/201014187) catalog for 1194 LMC clusters and 322 SMC clusters.

* `H03_duplicated_clusters.ods`

 List of duplicated entries in the H03 database.

* `DBs_cross_match_analy.ods`

 Analysis of the age differences between matched clusters (not important).

* `hunter_03.dat`

  Clusters in the S/LMC from [Hunter et al. (2003)](http://adsabs.harvard.edu/
  abs/2003AJ....126.1836H), 748 belong to the LMC and 191 to the SMC for a
  total of 939 clusters.

* `matched_clusters.dat`

 Final data file with all matched clusters between **ASteCA** and the
 databases in the literature.

* `matched_clusters_G10.dat`

 Idem above but using the G10 data as input to match with the P99, P00, and
 C06 databases, and compare their ages.

* `pietrz_00_LMC.dat`

  Ages for 600 LMC clusters obtained in [Pietrzynski & Udalski (2000)](http://
  adsabs.harvard.edu/abs/2000AcA....50..337P).

* `pietrz_99_SMC.dat`

 Ages for 93 SMC clusters obtained in [Pietrzynski & Udalski (1999)](http://
 adsabs.harvard.edu/abs/1999AcA....49..157P).

* `popescu_12_LMC.dat`

 List of ages and masses for 632 clusters in the database presented in
 [Popescu et al. (2012)](http://adsabs.harvard.edu/abs/2012ApJ...751..122P)
 with ages correlated to the *Hunter et al* catalog.

* `popescu_12_LMC_glatt.dat`

 Idem above, but also adds ages taken from the *Glatt et al.* catalog, for
 288 clusters in the *Hunter et al.* catalog.

* `rafelski_05_SMC.dat`

  Ages for 195 SMC clusters from [Rafelski & Zaritsky (2005)](http://
  adsabs.harvard.edu/abs/2005AJ....129.2701R), obtained via
  integrated photometry.


### `aux_funcs/`

* `add_p_vals/`

  Contains the script and input / output files used to add KDE p_values
  to the `asteca_output.dat` file from the 1st run, since the function was
  turned off.

* `move_files_names.py`

 Script to move .png files from their `input_XX/` folders for each run, into
 a single folder containing all the analyzed clusters.
 This script uses the `asteca_output_part_XX.dat` files.

* `move_files_sizes.py`

 Script to move files into `input_XX/` folders distributed so that each
 folder has approximately the same total size.

* `print_runs_clust_values.py`

 Prints the final parameters found in each run for a given cluster.


### `extinction_MCEV/`

* `cls_exts_match.dat`

 Output of the above script.

* `extin_analysis.py`

 Gets the extinction data from the table obtained via the MCEV service
 (`ra_dec_exts_mult_matches.dat`) and produces for each of the 210 clusters a
 value of the closest, average and maximum extinction values.

 The results are stored in the `cls_exts_match.dat` file.

* `IRSA_MC_ext.tbl`

  Output of the [IRSA](http://irsa.ipac.caltech.edu/applications/DUST/) query to
  obtain the [Schlafly & Finkbeiner
  (2011)](http://adsabs.harvard.edu/abs/2011ApJ...737..103S) corrected `E(B-V)`
  values for the 210 clusters.

* `IRSA_BB_ext.tbl`

  Idem above for the 3740 clusters in the Bica et al. catalog.

* `ra_dec_exts_mult_matches.dat`

  Output of the [MCEV](http://dc.zah.uni-heidelberg.de/mcextinct/q/cone/form)
  service query. Each of the 210 cluster's position is matched with at least
  one area in the reddening maps containing a E(V-I) value.

  The cluster [OHSC28]
  (http://simbad.u-strasbg.fr/simbad/sim-id?Ident=OHSC+28)
  needed a 6.0 deg search radius to find areas with extinction in the maps,
  making it the one furthest away from the areas studied in the MCEV maps.

* `TOPCAT_instruct.dat`

 Instructions to perform a query using TOPCAT.


### `figures/`

   Output figures from main script.


### `functions/`

   `*.py`: functions called by the main script.

   `0.004.dat` & `0.008.dat`: Marigo isochrones used by the
   `CMD_obs_vs_asteca.py` script to plot the G10 fitted isochrones.


### `runs/`

 Folder that contains the results for all the runs.

 `large_diff_age_clusts.ods`: analysis of the obtained values in different runs
 for those clusters that display a large age difference with the literature
 values.

#### `1st_run/`

  First batch of data obtained using the following parameters:

  - Semi mode.
  - Center found with 100px search area.
  - Auto radius.
  - Field regions: 10.
  - Decontamination algorithm: auto.
  - Reduced membership: Bayesian blocks binning
  - Restricted range in extinction: MCEV max + 0.1.
  - Restricted distance modulus:  SMC, [18.8, 19.201, 0.05] ;
    LMC, [18.3, 18.701, 0.05].
  - Best fit: Dolphin + Knuth binning.

#### `2nd_run/`

  Second batch of data for 83 clusters (60 LMC, 23 SMC) from the 1st run,
  marked to be re-processed. The data was obtained using the same extinction
  and distance modulus range, and changing the following parameters:

  * Center, radius, number of field regions, max E(B-V) and binning of red_memb
   and best_fit, as described in `README.dat`.

  The following clusters were left un-processed (no reason, just didn't bother):

  * NGC1917, SL579, SL588, LW54, NGC419, SL244

#### `3rd_run/`

  Third batch of data obtained using the following parameters, for all
  clusters:

  * Semi mode with center and radius fixed for several clusters.
  * Number of field regions set individually for some clusters.
  * Restricted range in extinction: **MCEV max for all clusters**.
  * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
    LMC, [18.4, 18.601, 0.02].

#### `4th_run/`

  Fourth run, 44 clusters from 3rd run re-processed. Same parameters as in
  3rd run but increased number of generations and mutation rate, and skipped
  the removal of stars post-DA.
  Some clusters had their parameters adjusted in `semi_input.dat`.

  * skip red_memb
  * gens=3000
  * mut_rate=0.2

#### `5th_run/`

  25 clusters from 4th run re-process. Increased mutation rate, three blocks
  of clusters processed with the following options:

  * mut_rate=0.25 (all blocks)
  * skip DA + red_memb: local + scott
  * skip DA + red_memb: local + blocks
  * skip DA + skip red_memb <-- Will fit mostly field stars if n_memb is low.

#### `6th_run/`

  Sixth run, 16 problematic clusters processed in three blocks with:

  * DA + red_memb=scott
  * skip DA + red_memb=scott
  * DA + red_memb=blocks (same as 3rd run)

#### `7th_run/`

  Seventh run, 20 high mass clusters re-process with higher mass max limit.
  Parameters used:

  * max mass = 30000
  * n_gen = 3000
  * mut = 0.25

#### `8th_run/`

  Re-process 8 clusters from 1st and 2nd run with extinction and distance
  modulus limits as in 3rd run, more generations and an increased mutation
  rate.

  * Restricted range in extinction: MCEV max for all clusters.
  * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
    LMC, [18.4, 18.601, 0.02]
  * n_gen = 3000
  * mut = 0.25

#### `9th_run/`

  Three clusters re-processed from the 8th run that showed large age
  differences, with the following parameters changed:

  * n_gen = 5000
  * n_el = 10
  * n_ei = 100
  * n_es = 25

#### `10th_run/`

  Re-process 13 cluster. Two are the ones were parameter values are still
  taken from the 2nd run, and the remaining 11 are those with
  `\delta log(age)>= 0.5`.

  Use the following parameters:

  * DA + red_memb = scott
  * n_gen = 5000
  * n_el = 1
  * n_ei = 25
  * n_es = 150

#### `11th_run/`

 Re-process the 2 clusters still dragged from the 2nd run (SL218 & BSDL654)
 plus 4 more clusters with large age differences (H88-131, BSDL631, L35, SL579).

 Use the following parameters:

 * DA + red_memb = different values per cluster
 * n_gen = 3000
 * p_mut = 0.1
 * n_el = 1
 * n_ei = 50
 * n_es = 30

#### `12th_run/`

 Re-process SL218 using the following parameters:

 * field regions = 25
 * DA + red_memb = auto + local (blocks)
 * n_pop = 500
 * n_gen = 2000
 * p_mut = 0.5
 * n_el = 1
 * n_ei = 25
 * n_es = 70

#### `13th_run`

 Re-process all clusters using a binary fraction of 0.

 Eight clusters from the SMC are missing (got lost somewhere in the FTP moving).
 These clusters will not be present in the runs: 14, 15, 16.

#### `14th_run`

  Re-process all clusters using fixed distance moduli (18.49/18.96, S/LMC) and
  increasing the maximum extinction by 0.05. Maximum mass increased to 30000.

#### `15th_run`

  Re-process all clusters using:

  * Range in extinction: **MCEV max + 0.1** for all clusters.
  * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
    LMC, [18.4, 18.601, 0.02]
  * Maximum mass: 30000.

#### `16th_run`

  Equal to 15th run but using a fixed mass of 1500 Mo and `Tolstoy` likelihood
  instead of `Dolphin`.
  **Correction**: I used  *whitout noticing* a mass interval of:
  `TM 10 1500` in the `params_input.dat` file instead of fixing it to a value
  of `1500`.
  This caused the mass values in the range to be defined by
  `np.arange(10, 1500)` as `[  10,   11,   12, ..., 1497, 1498, 1499]`, ie:
  a step of 1 solar mass.
