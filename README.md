# spt3g_cutter
Software for the SPT-3G thumbnail cutout service.

Setup path
----------
```
source /opt/miniconda3/bin/activate base
source ~/spt3g-devel/spt3g_cutter/setpath.sh ~/spt3g-devel/spt3g_cutter
```

Example 1:
----------
```
spt3g_cutter ~/spt3g-dummy.csv --dbname /data/spt3g/dblib/spt3g.db --outdir .  --date_start 2020-01-01 --date_end 2021-06-30  --bands 150GHz --np 16 --yearly yearly_winter_2020
```
