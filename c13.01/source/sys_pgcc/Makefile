# It looks like the Portland Group compilers v11.2 and 11.9 cannot handle parallel builds
# so disable them here by explicitly including -j 1 on the $(MAKE) command
FORCE:
%::
	$(MAKE) -j 1 -f ../Makefile SRCDIR=.. $(MAKECMDGOALS)
