

- [mc-catalog](#mc-catalog)
  - [Top level](#top-level)
    - [`databases/`](#databases)
    - [`AMRs/`](#AMRs)
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


### `databases/`

* `BB_10.dat`

 List of ages for 151 SMC and 539 LMC clusters in the database presented in
 [Bonato & Bica (2010)](http://cdsads.u-strasbg.fr/abs/2010MNRAS.403..996B).
 Only 285 out of the 690 clusters have an age value assigned.

* `McL_vdM05.dat`

 Masses from [McLaughlin & van der Marel (2005)](http://adsabs.harvard.edu/
 abs/2005ApJS..161..304M). Only three OCs could be cross matched: LMC-NGC1860,
 LMC-SL663, and SMC-NGC339.

* `pietrz_99_SMC.dat`

 Ages for 93 SMC clusters obtained in [Pietrzynski & Udalski (1999)](http://
 adsabs.harvard.edu/abs/1999AcA....49..157P).

* `pietrz_00_LMC.dat`

  Ages for 600 LMC clusters obtained in [Pietrzynski & Udalski (2000)](http://
  adsabs.harvard.edu/abs/2000AcA....50..337P).

* `hunter_03.dat`

  Clusters in the S/LMC from [Hunter et al. (2003)](http://adsabs.harvard.edu/
  abs/2003AJ....126.1836H), 748 belong to the LMC and 191 to the SMC for a
  total of 939 clusters.

* `rafelski_05_SMC.dat`

  Ages for 195 SMC clusters from [Rafelski & Zaritsky (2005)](http://
  adsabs.harvard.edu/abs/2005AJ....129.2701R), obtained via
  integrated photometry.

* `chiosi_06.dat`

 Contains 311 SMC clusters younger than 1 Gyr (ages determined using isochrone
 fitting) from the work by [Chiosi et al. (2006)](http://adsabs.harvard.edu/
 abs/2006A%26A...452..179C). Also lists several other structures/objects, ie:
 association, supernovas, H II regions and entities in between.

* `glatt_10.dat`

 CMD ages taken from the [Glatt et al. (2010)](http://www.aanda.org/10.1051/
 0004-6361/201014187) catalog for 1194 LMC clusters and 322 SMC clusters.

* `popescu_12_LMC.dat`

 List of ages and masses for 632 clusters in the database presented in
 [Popescu et al. (2012)](http://adsabs.harvard.edu/abs/2012ApJ...751..122P)
 with ages correlated to the *Hunter et al* catalog.

* `popescu_12_LMC_glatt.dat`

 Idem above, but also adds ages taken from the *Glatt et al.* catalog, for
 288 clusters in the *Hunter et al.* catalog.

* `cross_match.py`

 Matches the data in the above databases to those clusters processed by
 ASteCA.

* `matched_clusters.dat`

 Final data file with all matched clusters between **ASteCA** and the
 databases in the literature.

* `matched_clusters_G10.dat`

 Idem above but using the G10 data as input to match with the P99, P00, and
 C06 databases, and compare their ages.

* `cross_h03_p12.py`

 Matches the data in the H03 and P12 databases.

* `matched_H03_P12.dat`

 Crossed-matched OCs from H03 and P12.


### `AMRs/`

  Data files with literature AMR values.

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

* `IRSA_MC_ext.tbl`

  Output of the [IRSA](http://irsa.ipac.caltech.edu/applications/DUST/) query to
  obtain the [Schlafly & Finkbeiner
  (2011)](http://adsabs.harvard.edu/abs/2011ApJ...737..103S) corrected `E(B-V)`
  values for all the clusters.

* `IRSA_BB_ext.tbl`

  Idem above for the 3740 clusters in the Bica et al. catalog.

* `TOPCAT_instruct.dat`

 Instructions to perform a
 [MCEV](http://dc.zah.uni-heidelberg.de/mcextinct/q/cone/form) query using
 TOPCAT.

* `topcat-lite.jnlp`

 File that downloads and runs the lite version of TOPCAT (see instructions to
 run).

* `ra_dec.dat`

 RA & DEC data for all clusters in CSV format.

* `ra_dec_exts_mult_matches.dat`

 Output of the MCEV service query with TOPCAT. Each of the cluster's position is
 matched with at least one area in the reddening maps containing a E(V-I) value.

 The cluster [OHSC28](http://simbad.u-strasbg.fr/simbad/sim-id?Ident=OHSC+28)
 needed a 6.0 deg search radius to find areas with extinction in the maps,
 making it the one furthest away from the areas studied in the MCEV maps.

* `extin_analysis.py`

  Gets the extinction data from the table obtained via the MCEV service
  (`ra_dec_exts_mult_matches.dat`) and produces for each of the clusters a
  value of the closest, average and maximum E(B-V) extinction values.
  The transformation equation used is:

  `E(V-I) = 1.38 * E(B-V)`

  according to Tammann et al. 2003, A&A, 404, 423.

  The results are stored in the `cls_exts_match.dat` file.

* `cls_exts_match.dat`

   Output of the above script.


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

  For NGC339 (input_05) the max mass limit was set to 30000.

#### `6th_run/`

  Sixth run, 16 problematic clusters processed in three blocks with:

  * DA + red_memb=scott
  * skip DA + red_memb=scott
  * DA + red_memb=blocks (same as 3rd run)

#### `7th_run/`

  Seventh run, 20 high mass clusters re-processed with higher mass max limit.
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
  * No bootstrap process is applied.

#### `15th_run`

  Re-process all clusters using:

  * No bootstrap.
  * Range in extinction: **MCEV max + 0.1** for all clusters.
  * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
    LMC, [18.4, 18.601, 0.02]
  * Maximum mass: 30000.

#### `16th_run`

  Equal to 15th run but using a fixed mass of 1500 Mo and `Tolstoy` likelihood
  instead of `Dolphin`.
  **Correction**: I used  *without noticing* a mass interval of:
  `TM 10 1500` in the `params_input.dat` file instead of fixing it to a value
  of `1500`.
  This caused the mass values in the range to be defined by
  `np.arange(10, 1500)` as `[  10,   11,   12, ..., 1497, 1498, 1499]`, ie:
  a step of 1 solar mass.

#### `17th_run`

 Process 8 clusters with large age differences: BSDL631, H88-131, H88-316,
 KMHK975, KMHK979, SL218, SL579, L35.

 * Range in extinction: **MCEV max** for all clusters.
 * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
   LMC, [18.4, 18.601, 0.02]
 * Fixed mass: 1500.
 * `Tolstoy`.

#### `18th_run`

 *Important*: for all these clusters' data files the x,y columns are inverted.

 Process the 29 SMC clusters from Maia et al. (2014). Use the same parameters
 as in the 3rd run:

 * Semi mode with center and radius fixed for several clusters.
 * Number of field regions set individually for some clusters.
 * Restricted range in extinction: **MCEV max for all clusters**.
 * Restricted distance modulus:  SMC, [18.86, 19.061, 0.02] ;
   LMC, [18.4, 18.601, 0.02].
 * Total mass: 10000

Per cluster DA/GA parameters:

##### 1st run (29)
 1. 3000 generations
  * H86-97 (crashed, got not results)
 1. DA: auto + local + knuth
  * B103, H86-190, K47, NGC241, NGC242
 1. DA: auto + local + scott
  * B99
 1. DA: skip + local + knuth
  * B48, H86-76, H86-85, K61, SOGLE196
 1. DA: skip + local + scott
  * B55, H86-87, H86-90
 1. DA: skip + local + blocks
  * B134

###### Clusters kept from 1st run (12 + H86-188)
 * B103, B111, BS75, BS80, H86-76, H86-174, HW52, K43, K55, K57, K63, SOGLE196,
 H86-188 (it is processed several more times, but this one is kept)

##### 2nd run (17)
1. Defaults (DA: auto + local + blocks)
 * B124 (03), B55 (06), B99 (07), H86-188 (11), HW32 (18), NGC241 (27)
1. DA: skip + local + knuth
 * B134 (04)
1. DA: auto + local + knuth
 * B48 (05), H86-85 (14), H86-87 (15), H86-90 (16), H86-97 (17), K61 (24)
1. DA: auto + local + sqrt
 * H86-190 (12)
1. DA: auto + skip
 * K47 (21), L39 (26), NGC242 (28)

###### Clusters kept from 2nd run (6)
 * B55, B99, HW32, L39, K61, K47

##### 3rd run (11)
**All**: p_mut=0.35
1. Defaults params for DA
 * B48 (05), B124 (03), H86-188 (11), NGC241 (27), H86-97 (17), H86-190 (12),
   NGC242 (28)
1. DA: skip + local + knuth
 * H86-85 (14), H86-90 (16)
1. DA: skip + local + scott
 * B134 (04)
1. DA: skip + local + bb
 * H86-87 (15)

###### Clusters kept from 3rd run (2)
 * H86-190, NGC241

##### 4th run (9 + NGC241 is run again but 3rd run is kept)
**All**: p_mut=0.5, n_ei=25, n_es=60
1. Defaults params for DA
 * B124 (03), B134 (04), B48 (05), H86-188 (11), H86-87 (15)
1. DA: auto + local + knuth
 * H86-85 (14)
1. DA: auto + skip
 * H86-90 (16), H86-97 (17), NGC241 (27), NGC242 (28)

###### Clusters kept from 4th run (3)
 * B48, H86-87, B124

##### 5th run (6)
**All**: n_gen=3000, p_mut=0.5, n_ei=50, n_es=40
1. Defaults params for DA
 * H86-188 (11), H86-90 (16)
1. DA: auto + local + knuth
 * B134 (04)
1. DA: auto + skip
* H86-85 (14), H86-97 (17)
1. DA: skip + local + knuth
 * NGC242 (28)

###### Clusters kept from 5th run (4)
 * H86-85, NGC242, H86-97, B134

##### 6th run (3, B134 is processed again)
1. Defaults params for DA; GA: n_pop=200, n_gen=2000, p_mut=0.85, n_el=5
 * B134 (04), H86-188 (11), H86-90 (16)

###### Clusters kept from 6th run (1)
 * H86-90

##### 7th run (1)
1. Defaults params for DA; GA: n_pop=50, n_gen=2000, p_mut=0.5, n_el=1
 * H86-188 (11)

##### 8th run (1)
1. Defaults params for DA; GA: n_pop=200, n_gen=2000, p_mut=0.95, n_el=10
 * H86-188 (11)

##### 9th run (1)
1. DA: skip + local + knuth; GA: n_pop=100, n_gen=2000, p_mut=0.35, n_el=1
 * H86-188 (11)

##### 10th run (8)
Re-run 8 clusters that needed an offset applied to their photometry (which I
found out after processing them in the runs above): B48, H86-76, H86-90,
H86-97, L39, NGC241, NGC242, SOGLE196.
1. p_mut=0.35: NGC241
2. p_mut=0.5, DA skip, local+knuth: H86-76; p_mut=0.5, DA auto, skip: H86-97
3. p_muy=0.5, n_gen=3000, 

#### `19th_run`

Re-process these clusters to try to improve their fits. H88-131 could be better
fitted for the age and the other two for the metallicity.

1. DA: skip + local + blocks; GA: p_mut=0.35
* H88-131
1. DA: defaults; GA: p_mut=0.35
* NGC294, HW85

#### `20th_run`

Re-process L45 and L50 with better radii estimates. Originally, the values used
where:
rad_L45 = 162.5 (semi); rad_L50 = 170. (auto)

new values are:
rad_L45 = 100. (auto); rad_L50 = 100 (semi)

Since no substantial change is made, these new values are incorporated into the
3rd and 4th run where L50 and L45 where processed, respectively.
KDE probs are borrowed from the previous run, since they were not obtained
here.

Re-process NGC294 and HW85 fixing all parameters except metallicity, with
binary fractions of [0., 0.25, 0.5, 0.75, 1.]
Metallicity of NGC294 changed 0.001 --> 0.0005, due to issue #248.
These runs are only used to check the metal-binarity dependence, and are not
added to the final output file.

#### `21st_run`

Re-process 5 OCs whose masses differ greatly from the H03 and P12 estimates.
We use:
* no DA + no local removal of field stars
* no bootstrap
* max mass limit of 400000 with a step of 5000 (Mo)

All other parameter ranges are kept from 3rd run.

