bpch
====

GEOS-Chem Binary Punch File Reader/Plotter/ncdumper

Example Plotting
----------------

Examples:  
    
1. This example produce a Lat-Lon, time average, mean layer with a log color-scale.

   $ python bpch.py -g IJ-AVG-$ -v O3  -t mean -l mean --log ctm.bpch
   Successfully created ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png


2. This example produces a zonal mean, time average plot stopping at 20km (if BOXHEIGHT available) or the 21st layer.

   $ python bpch.py -g IJ-AVG-$ -v O3 -t mean -c mean --ymax 20 ctm.bpch
   Successfully created ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png


3. This example produces a 1st layer Latitude-Time Hovmoller Diagram. 

    $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -c mean ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowall_colmean.png


4. This example produces a 1st layer Longitude-Time Hovmoller Diagram. 

    $ python bpch.py -g IJ-AVG-$ -v O3 -l 0 -r mean ctm.bpch
    Successfully created ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png


5. This example would produce two Ox figures. Both from time 1 (default), but the first file from layer 1 and the second file from layer 2. The first figure has a minimum (max) value of 20 (60) and the second has a minimum (maximum) of 25 (65).

    $ python bpch.py -g IJ-AVG-$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2
    Successfully created ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png
    Successfully created ctm.bpch2_IJ-AVG_O3_time0_layer1_rowall_colall.png

6. This example would produce one Ox difference figure with a minimum of -2 and a maximum of 2.
    
    $ python bpch.py -d -g IJ-AVG-$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2
    Successfully created ctm.bpch-ctm.bpch2-diff_IJ-AVG_O3_time0_layer0_rowall_colall.png
