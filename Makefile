all:
	python bpch.py -g IJ-AVG-$$ -v O3 -t mean -l 0 ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3  -t mean -l mean --log ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -t mean -c mean --ymax 20 ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -l 0 -c mean ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -l 0 -r mean ctm.bpch
	python bpch.py -g IJ-AVG-$$ -v O3 -n 20 -x 60 -t 0 -l 0 -n 25 -x 65 -t 0 -l 1 ctm.bpch ctm.bpch2
	python bpch.py -d -g IJ-AVG-$$ -v O3 -n -2 -x 2 -t 0 -l 0 ctm.bpch ctm.bpch2
	python bpch.py -g NOX-LI-$ -v NOx -n 1e6 -x 1e12 -t mean -c sum ctm.bpch --log --ymax 20
	python bpch.py -g BIOGSRCE -v ISOP -n 1e6 -x 4e12 -l sum ctm.bpch --log