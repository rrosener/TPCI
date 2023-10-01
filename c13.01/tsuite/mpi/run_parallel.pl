#!/usr/bin/perl

# This is a little perl script that runs the mpi test suite on suitable
# machines or clusters. See the file auto/run_parallel.pl for a full
# description of how this script should be used.
#
# Peter van Hoof

do "../run_service.pm";

# this sets the path to various executables, sets up the list of CPU slots,
# and determines the number of CPUs to be used
$nproc = RunService::initialize( ".." );

# clean up the test suite
system "./clean_tsuite.pl";

# create list of input scripts to be executed
chomp( $jobList = `./run_sequence.pl` );

print "\nNow I will run the test suite using $nproc processors\n\n";

# this does the heavy lifting...
RunService::run_jobs_mpi( $jobList, $nproc );

# now do the checking
system "./checkall.pl";
