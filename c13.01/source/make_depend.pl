#!/usr/bin/perl
$command = $ARGV[0];
$sed_arg = $ARGV[1];
$sed_arg2 = $ARGV[2];
$outfile = $ARGV[3];
$temp_file = "$outfile.tmp";

# make sure that double quotes are escaped when passed to the shell
$command =~ s/"/\\"/g;

system( "$command > $temp_file" ) && &error_exit();
# this produces the dependencies for the object file
system( "cat $temp_file | sed \"$sed_arg\" > $outfile" ) && &error_exit();
# this produces the dependencies for the dependency file itself
# the full dependencies are needed to make sure that the dependency file
# is updated when included header files change, but not the source file itself...
system( "cat $temp_file | sed \"$sed_arg2\" >> $outfile" ) && &error_exit();
unlink $temp_file;

sub error_exit {
    unlink $temp_file, $outfile;
    exit 1;
}
