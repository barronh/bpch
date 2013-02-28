bpch
====

GEOS-Chem Binary Punch File Reader/Plotter/ncdumper

The python-based library for Binary Punch files is designed to provide NetCDF-like access to data in Binary Punch files. In addition, it can be used to produce plots.

Prerequisites
-------------

 * Linux<sup>`*`<sup>
 * Python >= 2.7
 * numpy >= 1.5
 * matplotlib >= 1.0<sup>`**`<sup>
 * mpl_toolkits.basemap >= 1.0<sup>`***`<sup>

<sup>`*`</sup>Linux is not actually required; Windows and Mac versions have been used. They are not supported.
<sup>`**`</sup>matplotlib is only required for plotting
<sup>`***`</sup>basemap is used only if available to add geopolitical boundaries to lat-lon plots.

Prerequisites via EPD
---------------------
This is the easiest way to get going. All of these prerequisites are available via the <a href=http://www.enthought.com/products/epd.php>Enthought Python Distribution</a> (EPD). EPD is free for academics, but the edX version is not sufficient (unless combined with virtualenv instructions below). 


Prerequisites via virtualenv
----------------------------
If you have a working version of Python 2.7 (works with 2.5, but not supported) but you don't have root access to install, then I recommend you use <a href=https://pypi.python.org/pypi/virtualenv>virtualenv</a> by downloading <a href=https://raw.github.com/pypa/virtualenv/master/virtualenv.py>virtualenv.py</a>. The commands below setting up virtualenv and installing prerequisites:

```
cd ~
curl -LO https://raw.github.com/pypa/virtualenv/master/virtualenv.py
python virtualenv.py aqc
source aqc/bin/activate
pip install numpy
pip install matplotlib
pip install basemap
```

Any time you want to use bpch, you'll need to activate virtualenv by running `source bpch/bin/activate`. You'll need to `source bpch/bin/activate`.

Example Plotting
----------------

1. This example produce a Lat-Lon, time average, mean layer with a log color-scale.

    ```
    $ python bpch.py -g IJ-AVG-$ -v O3  -t mean -l mean --log ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png
    ```

2. This example produces a zonal mean, time average plot stopping at 20km (if BOXHEIGHT available) or the 21st layer.

    ```
    $ python bpch.py -g IJ-AVG-$ -v O3 -t mean -c mean --ymax 20 ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png
    ```

3. This example produces a 1st layer Latitude-Time Hovmoller Diagram. 

    ```
    $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -c mean ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowall_colmean.png
    ```

4. This example produces a 1st layer Longitude-Time Hovmoller Diagram. 

    ```
    $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -r mean ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png
    ```

5. This example would produce two Ox figures. Both from time 1 (default), but the first file from layer 1 and the second file from layer 2. The first figure has a minimum (max) value of 20 (60) and the second has a minimum (maximum) of 25 (65).

    ```
    $ python bpch.py -g IJ-AVG-$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2
    Successfully created ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png
    Successfully created ctm.bpch2_IJ-AVG_O3_time0_layer1_rowall_colall.png
    ```
    
6. This example would produce one Ox difference figure with a minimum of -2 and a maximum of 2.
    
    ```
    $ python bpch.py -d -g IJ-AVG-$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2
    Successfully created ctm.bpch-ctm.bpch2-diff_IJ-AVG_O3_time0_layer0_rowall_colall.png
    ```

Example Python
--------------
```
from bpch import bpch

bcfile = bpch('ctm.bpch')
print gcfile.groups.keys()
group = bcfile.groups['IJ-AVG-$']
print group.variables.keys()
ozone = group.variables['O3']
print ozone.dimensions
print ozone.units
print ozone.mean(0).reshape(ozone.shape[1], -1).mean(1)
...

```