use strict;
use warnings;
use Data::Dumper;
use IO::File;
use Benchmark;
use Time::HiRes qw(usleep ualarm gettimeofday tv_interval);

use Getopt::Long;
use File::Basename;
use String::Diff;


#  perl statistics_bootstrapping.pl -input_dir_prefix /home/scott/Documents/data/alonso_lab/repeat_v3_fast_35_50_144_run| grep '^SBS'  > /home/scott/Documents/data/alonso_lab/repeat_v3_fast_35_50_144_run01/synth_op.txt


#			SBS FBgn0029588 FBgn0040038 10
#			SBS FBgn0029588 FBgn0040368 10
#			SBS FBgn0029588 FBgn0040892 10
#			SBS FBgn0029589 FBgn0029588 10
#			SBS FBgn0029589 FBgn0040038 10
#			SBS FBgn0029589 FBgn0040368 10
#			SBS FBgn0029589 FBgn0040892 10
#			SBS FBgn0040038 FBgn0029588 10
#			SBS FBgn0040038 FBgn0029589 10
#			SBS FBgn0040038 FBgn0040368 10
#			SBS FBgn0040038 FBgn0040892 10
#			SBS FBgn0040345 FBgn0029580 10
#			SBS FBgn0040347 FBgn0000826 10
#			SBS FBgn0040348 FBgn0023534 10
#			SBS FBgn0040364 FBgn0025393 10


#
#  here is how to run this script from command line
#
#  perl statistics_bootstrapping.pl -input_dir_prefix /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run
#
#  where we have previously executed predict_microRNA_CoRegulated_Genes.pl
#  using this mode :  #my @preset_testing_fast   = ( 35, 50, 144,  0, 0, 0,   $NO,  $NO,  $YES, $YES );  #    $YES    $NO

&main();

# <><><>         end of mainline logic           <><><> #


#		ls /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run*/stand*
#		
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run01/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run02/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run03/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run04/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run05/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run06/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run07/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run08/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run09/standard_out.txt
#		/home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run10/standard_out.txt

# how is this possible

                                                                                           
#          'FBgn0040344' => {                                                                                             
#                             'FBgn0040351' => 1,                                                                         
#                             'FBgn0040373' => 5,                                                                         
#                             'FBgn0039955' => 1,                                                                         
#                             'FBgn0029580' => 1,                                                                         
#                             'FBgn0027794' => 6,                                                                         
#                             'FBgn0040022' => 10,                                                                        
#                             'FBgn0000022' => 1,                                                                         
#                             'FBgn0020381' => 3,                                                                         
#                             'FBgn0040371' => 1,                                                                         
#                             'FBgn0040348' => 2,                                                                         
#                             'FBgn0025383' => 1,                                                                         
#                             'FBgn0004648' => 3,                                                                         
#                             'FBgn0040365' => 1,                                                                         
#                             'FBgn0040350' => 1,
#                             'FBgn0029525' => 17,   <-- where is count of 17 coming from 
#                             'FBgn0026876' => 5,
#                             'FBgn0016038' => 1,
#                             'FBgn0021764' => 1,
#                             'FBgn0000210' => 2,
#                             'FBgn0025616' => 1,
#                             'FBgn0025394' => 1,
#                             'FBgn0040345' => 1,
#                             'FBgn0023534' => 2
#                           }




#		cat  /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run*/stand* | grep '^CRA' |grep FBgn0040344
#		CRA curr_cluster 0 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 20 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 2 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 31 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 28 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 31 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 9 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 12 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 25 KP00131-405 FBgn0040344 curr_gene 79
#		CRA curr_cluster 29 KP00131-405 FBgn0040344 curr_gene 79



#		 cat  /home/scott/Documents/data/alonso_lab/repeat_fast_35_50_144_run*/stand* | grep '^CRA' |grep FBgn0029525
#		CRA curr_cluster 0 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 0 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 20 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 20 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 2 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 2 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 31 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 31 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 26 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 28 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 31 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 31 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 9 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 9 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 12 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 12 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 25 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 25 KP00076-750 FBgn0029525 curr_gene 36
#		CRA curr_cluster 27 KP00075-750 FBgn0029525 curr_gene 35
#		CRA curr_cluster 27 KP00076-750 FBgn0029525 curr_gene 36


#  these are good pairings :
#
#         'FBgn0028697' => {                                                                                                                       
#                             'FBgn0025637' => 9,                                                                                                   
#                             'FBgn0002579' => 10,                                                                                                  
#                             'FBgn0015288' => 10,                                                                                                  
#                             'FBgn0003517' => 8                                                                                                    
#                           },   

sub main {
	
	my $time0 = new Benchmark;
    my $verbose = 0;                    # frequently referred
    my $debug = 0;                      # frequently referred
    my %h = ('verbose' => \$verbose, 'debug' => \$debug);
    GetOptions (\%h, 'verbose', 'debug', 'filter', 'size=i', 'input_dir_prefix=s');

    if ( ! exists $h{input_dir_prefix} ) {
		
		print "ERROR - must supply input parm 'input_dir_prefix' \n";
	}
	my $input_dir_prefix = $h{input_dir_prefix};
	
	my %hash_file_contents;
	
	my $input_dirname = dirname( $input_dir_prefix);
#	print "input_dirname $input_dirname\n";
	
	my $input_basename = basename( $input_dir_prefix);
#	print "input_basename $input_basename\n";
	
	opendir(RUNDATA, $input_dirname) || die("ERROR - cannot open directory $input_dirname $!\n"); 	
	foreach my $curr_dir (sort readdir(RUNDATA)) {
		
#		print "curr_dir $curr_dir\n";
		
		if ( $curr_dir =~ m/$input_basename/ ) {
			
#			print "seeing similar dir $curr_dir to given prefix of $input_basename\n";
			
			my $full_pathname = "$input_dirname/$curr_dir/standard_out.txt";
			
			print "$full_pathname\n";
			
			my $fh_full_pathname = IO::File->new( "< $full_pathname" )
				or die "ERROR - cannot open file full_pathname " .
					"$full_pathname : $!\n";	
					
			read_file_contents( $fh_full_pathname, $curr_dir, $input_basename, \%hash_file_contents);
		}
	}
	closedir(RUNDATA); 
	
	print_Dumper ( "hash_file_contents", \%hash_file_contents);
	
	print_hash_side_by_side("hash_file_contents", \%hash_file_contents);
				
	my $time1 = new Benchmark;
	my $interval_in_secs = tv_interval($time0, $time1);
	
	print "\n\nEnd of Computation ... completed in $interval_in_secs seconds\n";
	
}  #  main

sub print_hash_side_by_side {
	
	my ( $label, $ref_hash_file_contents ) = @_;
	
	foreach my $curr_key ( sort keys %{$ref_hash_file_contents}) {
		
		my $ref_hash_curr_value = $ref_hash_file_contents->{$curr_key};
				
		foreach my $curr_inner_key ( sort keys %{$ref_hash_curr_value}) {
			
			my $inner_value = $ref_hash_curr_value->{$curr_inner_key};
			
			print "SBS $curr_key $curr_inner_key $inner_value\n";	
		}
	}
}


sub read_file_contents {

	my ($fh_full_pathname, $curr_file, $input_basename, $ref_hash_file_contents) = @_;

	my $diff = String::Diff::diff_fully($curr_file,$input_basename);
	
	my $num_run = '';
	
	for my $line (@{ $diff->[0] }) {
		$num_run =  "$line->[1]\n";
	}
	chomp($num_run);
	print "num_run $num_run\n";
	
	my $prev_num_cluster = -1;
	my @array_sibling_genes_curr_cluster = ();  #  initialize array 
	
	while (<$fh_full_pathname>) {
		
		chomp;
		my $curr_line = $_;
		
		if ( $curr_line =~ m/^CRA /) {
		
#			CRA curr_cluster 27 KP00194-750 FBgn0027794 curr_gene 127
#			CRA curr_cluster 28 KP00048-228 FBgn0039993 curr_gene 13 
		
		
#			print "$curr_line\n";
		
		my ( $stubCRA,  $stub_curr_cluster, $num_cluster, $KP_id, $FB_id ) = split(/ /, $curr_line);
		
		if ( $num_cluster != $prev_num_cluster ) {
			
			print "seeing new cluster $num_cluster .........\n";
			
			process_this_cluster( $ref_hash_file_contents, \@array_sibling_genes_curr_cluster, $num_run);
			
			@array_sibling_genes_curr_cluster = ();  #  initialize array to emptiness
			
			$prev_num_cluster = $num_cluster;
		}


		print " $stubCRA  $stub_curr_cluster  $num_cluster  $KP_id  $FB_id\n";
#		print " $stubCRA  $stub_curr_cluster\n";
			
			push @array_sibling_genes_curr_cluster, $FB_id;
			
		}	
	}
	process_this_cluster( $ref_hash_file_contents, \@array_sibling_genes_curr_cluster);
	
	
}  #  read_file_contents

sub process_this_cluster {
	
	my ( $ref_hash_file_contents, $ref_array_sibling_genes_curr_cluster, $num_run) = @_;
	
	print "inside process with num_run $num_run and no more \n";
	
#	print_Dumper("ref_array_sibling_genes_curr_cluster",  $ref_array_sibling_genes_curr_cluster);
	
	foreach my $curr_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
		
			print "curr_FB_id $curr_FB_id\n";
				
			foreach my $curr_sibling_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
				
				next if ( $curr_sibling_FB_id eq $curr_FB_id);  #  skip if seeing self

				$ref_hash_file_contents->{$curr_FB_id}->{$curr_sibling_FB_id}++;
			}
	}
		
}  #  process_this_cluster


#
#sub process_this_cluster {
#	
#	my ( $ref_hash_file_contents, $ref_array_sibling_genes_curr_cluster, $num_run) = @_;
#	
#	print "inside process with num_run $num_run and no more \n";
#	
##	print_Dumper("ref_array_sibling_genes_curr_cluster",  $ref_array_sibling_genes_curr_cluster);
#	
#	foreach my $curr_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
#		
#			print "curr_FB_id $curr_FB_id\n";
#			
#			my %hash_sibling_FB_id_counts;
#			if ( exists $ref_hash_file_contents->{$curr_FB_id}) {
#				
#				%hash_sibling_FB_id_counts = %{$ref_hash_file_contents->{$curr_FB_id}};
#			}
#			
#			foreach my $curr_sibling_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
#				
#				next if ( $curr_sibling_FB_id eq $curr_FB_id);  #  skip if seeing self
#				
#				my $count_saw_this_sibling = 0;
#				if (exists $hash_sibling_FB_id_counts{$curr_sibling_FB_id}) {
#					
#					$count_saw_this_sibling = $hash_sibling_FB_id_counts{$curr_sibling_FB_id};
#				}
#				$hash_sibling_FB_id_counts{$curr_sibling_FB_id} = ++$count_saw_this_sibling;
#			}
#			$ref_hash_file_contents->{$curr_FB_id} = \%hash_sibling_FB_id_counts;
#	}
#		
#}  #  process_this_cluster



sub print_Dumper {
	
	my ($label, $to_dump ) = @_;
	
	print "$label\n";
	print Dumper ( $to_dump);
}


