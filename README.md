# spt3g_cutter
Software for the SPT-3G thumbnail cutout service.

To install:
-----------
```
INSTALL_PATH=/opt/spt
SPT3G_CUTTER_VERSION=0.4.6
cd /tmp
git clone https://github.com/CMB-NCSA/spt3g_cutter -b $SPT3G_CUTTER_VERSION
sudo mv -v /tmp/spt3g_cutter $INSTALL_PATH
```

To install in `/usr/local`
```
sudo pip install .  --prefix=/usr/local
```

Setup path
----------
```
INSTALL_PATH=/opt/spt
source $INSTALL_PATH/spt3g_cutter/setpath.sh $INSTALL_PATH/spt3g_cutter
```

Example 1:
----------
```
spt3g_cutter ~/spt3g-dummy.csv --dbname /data/spt3g/dblib/spt3g.db --outdir .  --date_start 2020-01-01 --date_end 2021-06-30  --bands 150GHz --np 16 --yearly yearly_winter_2020
```
