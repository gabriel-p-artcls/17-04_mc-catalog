# mc-catalog

Set of scripts to analyze and plot the results of the **ASteCA** processing of
the 210 MC clusters in our database, observed with Washington photometry.


## Scripts and data files description

* `bb_cat.dat`

  RA & DEC positions for the 3740 clusters in the
  [Bica et al. (2008)](http://cdsads.u-strasbg.fr/abs/2008MNRAS.389..678B)
  catalog ([Table 3](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/MNRAS/389/678/table3)).

* `mc_cat_analysis.py`

 Script that Produces the main figures by comparing parameters with each other
 (KDE maps), with literature values (1:1 plots) and distributed in space
 (RA vs DEC plots.)

 Functions are stored in the `functions/` folder.

* `extinction_MCEV/IRSA_MC_ext.tbl`

  Output of the [IRSA](http://irsa.ipac.caltech.edu/applications/DUST/) query to
  obtain the [Schlafly & Finkbeiner
  (2011)](http://adsabs.harvard.edu/abs/2011ApJ...737..103S) corrected `E(B-V)`
  values for the 210 clusters.

* `extinction_MCEV/IRSA_BB_ext.tbl`

  Idem above for the 3740 clusters in the Bica et al. catalog.

* `extinction_MCEV/ra_dec_exts_mult_matches.dat`

  Output of the [MCEV](http://dc.zah.uni-heidelberg.de/mcextinct/q/cone/form)
  service query. Each of the 210 cluster's position is matched with at least
  one area in the reddening maps containing a E(V-I) value.

  The cluster [OHSC28]
  (http://simbad.u-strasbg.fr/simbad/sim-id?Ident=OHSC+28)
  needed a 6.0 deg search radius to find areas with extinction in the maps,
  making it the one furthest away from the areas studied in the MCEV maps.

* `extinction_MCEV/extin_analysis.py`

 Gets the extinction data from the table obtained via the MCEV service
 (`ra_dec_exts_mult_matches.dat`) and produces for each of the 210 clusters a
 value of the closest, average and maximum extinction values.

 The results are stored in the `cls_exts_match.dat` file.

* `extinction_MCEV/cls_exts_match.dat`

  Output of the above script.
