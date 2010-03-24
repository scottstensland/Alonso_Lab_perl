#!/usr/bin/perl

# helper code to implement cdd bash alias found in ~/.bashrc 
# to allow this syntax at command line :
#
#           cdd   run03 run04
#
#     to effect a pwd change from current pwd which contains string run03
#     into a cd change into a regex of run04
#
#     for example:  our pwd is 
#
#  /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run03
#
#  and we desire to cd into 
#
# /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run04
#
#      so we just issue    cdd run03 run04
#
#      or more cheeky just issue    cd 03 04  #  since we have mulitple 3 in source dir 03 uniquely specifies desired 3
#
#      scott stensland - march 2010
#
#   here are upstream code N logic
#
#   cat /home/scott/cdd_inner
#   cd `/home/scott/workspace_eclipse_helios/Alonso_Lab_perl/cdd.pl $1 $2 `
#
#   and in ~/.bashrc we have
#
#   alias cdd=' source /home/scott/cdd_inner $1 $2 '



my $new_dir_to_cd_into = $ENV{PWD};  #  copy current dir into preliminary target dir to cd into

$new_dir_to_cd_into =~ s/$ARGV[0]/$ARGV[1]/;

if ( -d $new_dir_to_cd_into ) {
	
	print "$new_dir_to_cd_into";
		
} else {
	
	print $ENV{PWD};
}