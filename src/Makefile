SHELL := bash

headers := $(wildcard *.h)
utilities := $(wildcard */*.py)

.PHONY : clean clean_doc devHTML devPDF doc pdf test userHTML userPDF

################################################################################
# Build tests.

test :
	$(MAKE) -C test/

################################################################################
# Documentation building.

devHTML := doc/developer/html/index.html
userHTML := doc/user/html/index.html
devPDF := doc/developer/latex/refman.pdf
userPDF := doc/user/latex/refman.pdf

# Make Doxygen HTML manuals for users and developers.
doc : $(devHTML) $(userHTML)
devHTML: $(devHTML)
userHTML: $(userHTML)
$(devHTML) : doc/Doxyfile_developer doc/mainpage.dox $(headers) $(utilities)
	cd doc && doxygen Doxyfile_developer
$(userHTML) : doc/Doxyfile_user doc/mainpage.dox $(headers) $(utilities)
	cd doc && doxygen Doxyfile_user
# Convert README.org into Doxygen format.  Insert a begin-comment and \mainpage
# command; change Org-mode sections to Doxygen sections; remove Org-mode source
# block delimiters; insert an end-comment.
doc/mainpage.dox : README.org doc/org2dox.py
	doc/org2dox.py < $< > $@

# Also make PDF versions of the manuals.
pdf : $(devPDF) $(userPDF)
devPDF : $(devPDF)
userPDF : $(userPDF)
$(devPDF) : $(devHTML)
	$(MAKE) -C doc/developer/latex/
$(userPDF) : $(userHTML)
	$(MAKE) -C doc/user/latex/

# Clean Doxygen outputs.
clean_doc :
	$(RM) -r doc/mainpage.dox doc/developer doc/user

################################################################################
# Code parsing.

BROWSE : $(headers)
	ebrowse $(headers)
TAGS : $(headers)
	etags $(headers)

################################################################################
# Remove everything.

clean : clean_doc
	$(RM) BROWSE TAGS
	$(MAKE) -C test/ clean
