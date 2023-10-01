#!/usr/bin/perl -w
# compile Quick Start Guide and all three parts of Hazy with pdflatex
# do multiple times to resolve bibtex entries
#
chdir( "hazy1" ) or die "invalid directory hazy1\n";
system( "pdflatex hazy1" );
system( "bibtex hazy1" );
system( "pdflatex hazy1" );
system( "pdflatex hazy1" );
#
chdir( "../hazy2" ) or die "invalid directory hazy2\n";
system( "pdflatex hazy2" );
system( "bibtex hazy2" );
system( "pdflatex hazy2" );
system( "pdflatex hazy2" );
#
chdir( "../hazy3" ) or die "invalid directory hazy3\n";
system( "pdflatex hazy3" );
system( "bibtex hazy3" );
system( "pdflatex hazy3" );
system( "pdflatex hazy3" );
#
chdir( "../QuickStart" ) or die "invalid directory hazy3\n";
system( "pdflatex QuickStart" );
system( "bibtex QuickStart" );
system( "pdflatex QuickStart" );
system( "pdflatex QuickStart" );
