all:
	robodoc
	# remove docs directory to remove recursions
	rm -rf docs

realclean:
	find . -name "*.html" -exec rm {} \; -print
	find . -mindepth 1 -maxdepth 1 -type d -regex "./[0-9,a-z,A-Z]*" -exec rm -rf {} \; -print
	rm -rf dl-find newscf_f77
