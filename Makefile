all: \
ctm.bpch_IJ-AVG_O3_timemean_layer0_rowall_colall.png \
ctm.bpch_IJ-AVG_O3_timeanimate_layer0_rowall_colall.mp4 \
ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png \
ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png \
ctm.bpch_IJ-AVG_O3_time_layer0_rowall_colmean.png \
ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png \
ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png \
ctm.bpch_IJ-AVG_O3_time0_layer1_rowall_colall.png \
ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png \
ctm.bpch2_IJ-AVG_O3_time0_layer0_rowall_colall.png \
ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colsum.png \
ctm.bpch_IJ-AVG_O3_timesum_layersum_rowall_colall.png


ctm.bpch_IJ-AVG_O3_timemean_layer0_rowall_colall.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -t mean -l 0 ctm.bpch

ctm.bpch_IJ-AVG_O3_timeanimate_layer0_rowall_colall.mp4: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -t animate -l 0 ctm.bpch

ctm.bpch_IJ-AVG_O3_timemean_layermean_rowall_colall.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3  -t mean -l mean --log ctm.bpch

ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colmean.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -t mean -c mean --ymax 20 ctm.bpch

ctm.bpch_IJ-AVG_O3_time_layer0_rowall_colmean.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -l 0 -c mean ctm.bpch

ctm.bpch_IJ-AVG_O3_timeall_layer0_rowmean_colall.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -l 0 -r mean ctm.bpch

ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png: bpch.py ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2

ctm.bpch_IJ-AVG_O3_time0_layer1_rowall_colall.png: ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png

ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png: bpch.py ctm.bpch
	python bpch.py -d -g IJ-AVG-$$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2

ctm.bpch2_IJ-AVG_O3_time0_layer0_rowall_colall.png: ctm.bpch_IJ-AVG_O3_time0_layer0_rowall_colall.png

ctm.bpch_IJ-AVG_O3_timemean_layerall_rowall_colsum.png: bpch.py ctm.bpch
	python bpch.py -g NOX-LI-$$ -v NOx -n 1e6 -x 1e12 -t mean -c sum ctm.bpch --log --ymax 20

ctm.bpch_IJ-AVG_O3_timesum_layersum_rowall_colall.png: bpch.py ctm.bpch
	python bpch.py -g BIOGSRCE -v ISOP -n 1e6 -x 4e12 -t sum -l sum ctm.bpch --log
