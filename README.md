# mc-catalog

Set of scripts to analyze and plot the results of the ASteCA processing of the
210 MC clusters observed with Washington photometry.


## mc_cat_analysis

Produces the main figures comparing parameters with each other (KDE maps), 
with literature values (1:1 plots) and distributed in space (RA vs DEC plots.)


## extinction_MCEV/extin_analysis

Gets the extinction data from the table obtained via the MCEV service and
produces for each of the 210 clusters a value of the closest, average and
maximum extinction values.

The results are stored in the `cls_exts_match.dat` file.