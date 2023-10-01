#!/usr/bin/perl

# This script tests for identical output among individual sims in "grid repeat" runs.

$out = $ARGV[0];
$delimiter = "GRID_DELIMIT -- grid";

# grep "Cloudy was called" and grab the number of times
$num_sims = qx/grep "$delimiter" $out | wc -l/;
chomp($num_sims);

# delete any old split files
while( defined( $input = glob("xx*") ) )
{
	unlink "$input";
}

# get rid of "duplicate reaction", "reverse reaction" warnings and "Using STOUT model" messages
# strip the "xxxxxxxxx" from "GRID_DELIMIT -- gridxxxxxxxxx"
system("egrep -v 'Warning: duplicate reaction|Warning! No reverse reaction|Using STOUT model' $out | sed 's/grid[0-9]\\{9\\}/grid/' > tmpfile");
# split output file and count number of resultant files
# should be number of sims plus 2 for header and footer
$command = "csplit -k tmpfile \'\/GRID_DELIMIT\/+1\' \'{$num_sims}\' | wc -l";
my $num_split = qx/$command/;
`unlink tmpfile`;

# now substract the 2
$num_split -= 2;

# make sure $num_sims and $num_split agree
if( $num_sims != $num_split )
{
	die "PROBLEM These should agree: ".$num_sims." and ".$num_split."\n";
}

# find how many xx* files there are
`unlink xx00`;
# generate filename for footer
$footer = "xx".sprintf( "%02i",$num_sims+1 );
`unlink $footer`;

# generate the filename for the last sim
$lastsim = "xx".sprintf( "%02i",$num_sims );

$res = 0;
while( defined( $input = glob("xx*") ) )
{
	# Compare everything to the last sim.
	# This minimizes output if first is different from all others which are themselves identical.
	# That is the most common pathology.
	$res += system( "diff -q $lastsim $input" );
}
if( $res == 0 )
{
    print "all sims are equal\n";
}
else
{
    print "PROBLEM: some sims did not repeat exactly\n";
}
