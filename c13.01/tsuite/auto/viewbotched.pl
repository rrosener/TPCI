#!/usr/bin/perl -w

# view all models with botched asserts  

# editor
$ED = "vi";

$nMod = 0;
while(defined($output= glob("*.out")) )
{
	#system "grep 'Cloudy ends'  $output >checkend.txt";
        system "grep Botched  $output >checkend.txt";
        system "grep 'W-'  $output >>checkend.txt";

#	if non-zero length, model had botched asserts
#       -s option returns number of bytes if file exists and has finite size
	if(-s "checkend.txt")
	{
           ++$nMod;

#         actually view the output
          system "$ED $output";
	}
}
printf( "\n %i models had botched asserts, and were viewed.\n ",$nMod);

