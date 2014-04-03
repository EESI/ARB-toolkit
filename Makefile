DESTDIR = /usr/local/

all: 

install: 
	@cp -v src/arb_import_pipeline.sh $(DESTDIR)/bin/
	@cp -v src/addUIDtoFasta.py $(DESTDIR)/bin/
	@cp -v src/build_ift_from_metalabels.py $(DESTDIR)/bin/
	@cp -v src/extract_leaf_names.py $(DESTDIR)/bin/
	@cp -v src/getAccession.py $(DESTDIR)/bin/
	@cp -v src/rename_tree_leaves.py $(DESTDIR)/bin/
	@echo you need to deal with your src/GDEmenu
