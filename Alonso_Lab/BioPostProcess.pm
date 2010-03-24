
package BioPostProcess;

use strict;
use warnings;
use Data::Dumper;
use IO::File;
use Getopt::Long;
use File::Basename;
use String::Diff;

use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(post_process );
our %EXPORT_TAGS = ( DEFAULT => [qw(&post_process)],
					);

my $YES = 'YES';
my $NO  = 'NO';

my $fh_standard_out;     #  standard out is send into this supplied file handle
my $send_print_to_file = $NO;  #  flag to determine whether to send standard out to output file

#   $| = 1;  #  turn buffering off

    sub new {

    	my $proto = shift;
        my $class = ref($proto) || $proto;
    	
    	my $self = {
    		
#    		_PART_NUM => undef,
#    		_QTY_ON_HAND => undef,
#    		_firstName => undef,
#    		_lastName => undef,
    		
    	};
		
		bless ($self, $class);
		return $self;
    }
    
    sub set_standard_out_fh {
        my $self = shift;
        my $value = shift;
        $fh_standard_out = $value;
	}
	
    sub set_send_print_to_file {
        my $self = shift;
        my $value = shift;
        $send_print_to_file = $value;
	}
	
	
	
	

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

sub post_process {
	
	my ($self, $input_dir_prefix, $fh_gene_pairs_count_file ) = @_;
	
	my %hash_file_contents;
	my %hash_file_contents_per_run;
	
	my %hash_genes_2_count_sibling_genes_per_cluster;  #  key gene/value count of each sibling gene sharing cluster
	
	
	opendir(RUNDATA, $input_dir_prefix) || die("ERROR - cannot open input_dir_prefix $input_dir_prefix $!\n");
	print "just opened dir $input_dir_prefix\n";
	
	foreach my $which_run_num_dir (sort readdir(RUNDATA)) {
		
		next if ($which_run_num_dir eq '.' || $which_run_num_dir eq '..' );  #  skip over . and ..
		next if ( ! -d "$input_dir_prefix/$which_run_num_dir");  #  skip if not a dir
		
		my $curr_outer_dir = "$input_dir_prefix/$which_run_num_dir";
		
		print "raw  curr_outer_dir $curr_outer_dir\n";
		
		my $full_pathname = "$curr_outer_dir/standard_out.txt";
		
		if ( ! -f $full_pathname ) { die "ERROR do not see $full_pathname $full_pathname\n"; }
		
		print "$full_pathname\n";
		
		my $fh_full_pathname = IO::File->new( "< $full_pathname" )
			or die "ERROR - cannot open file full_pathname " .
				"$full_pathname : $!\n";	
				
		read_file_contents( $fh_full_pathname, $curr_outer_dir, \%hash_file_contents,
							\%hash_file_contents_per_run,
							$full_pathname, $which_run_num_dir);
	}
	closedir(RUNDATA); 
	
	print_Dumper ( "hash_file_contents_per_run", \%hash_file_contents_per_run);
	#	print_hash_side_by_side("hash_file_contents", \%hash_file_contents, $fh_gene_pairs_count_file);
	
	do_bootstrapping_judgement( \%hash_file_contents_per_run, \%hash_genes_2_count_sibling_genes_per_cluster );
	
	print_Dumper("hash_genes_2_count_sibling_genes_per_cluster", \%hash_genes_2_count_sibling_genes_per_cluster);
	
}  #  post_process

sub do_bootstrapping_judgement {
	
	my ( $ref_hash_file_contents_per_run, $ref_hash_genes_2_count_sibling_genes_per_cluster ) = @_;
		
	foreach my $which_run_num_dir (keys %{$ref_hash_file_contents_per_run}) {
		
		my $ref_hash_all_clusters_this_run = $ref_hash_file_contents_per_run->{$which_run_num_dir};
		
		print_Dumper("ref_hash_all_clusters_this_run", $ref_hash_all_clusters_this_run);
		
		foreach my $num_cluster (keys %{$ref_hash_all_clusters_this_run}) {
			
			my @array_sibling_genes_curr_cluster = @{$ref_hash_all_clusters_this_run->{$num_cluster}};
			
			print "$which_run_num_dir $num_cluster array_sibling_genes @array_sibling_genes_curr_cluster\n";
			
			foreach my $parent_gene (@array_sibling_genes_curr_cluster) {
				
				my $ref_hash_sibling_genes_2_count;  #  keys sibling genes/value count time in same cluster as parent key
				
				if (exists $ref_hash_genes_2_count_sibling_genes_per_cluster->{$parent_gene}) {
					
					$ref_hash_sibling_genes_2_count = $ref_hash_genes_2_count_sibling_genes_per_cluster->{$parent_gene};
				}
				
				foreach my $sibling_gene (@array_sibling_genes_curr_cluster) {
						
					$ref_hash_sibling_genes_2_count->{$sibling_gene}++;
				}
				$ref_hash_genes_2_count_sibling_genes_per_cluster->{$parent_gene} = $ref_hash_sibling_genes_2_count;
			}
		}
	}
}


sub print_hash_side_by_side {
	
	my ( $label, $ref_hash_file_contents, $fh_gene_pairs_count_file ) = @_;
	
	foreach my $curr_key ( sort keys %{$ref_hash_file_contents}) {
		
		my $ref_hash_curr_value = $ref_hash_file_contents->{$curr_key};
				
		foreach my $curr_inner_key ( sort keys %{$ref_hash_curr_value}) {
			
			my $inner_value = $ref_hash_curr_value->{$curr_inner_key};
			
			print "SBS $curr_key $curr_inner_key $inner_value\n";	
			print $fh_gene_pairs_count_file "SBS $curr_key $curr_inner_key $inner_value\n";	
		}
	}
}

sub read_file_contents {

	my ($fh_full_pathname, $curr_file, $ref_hash_file_contents, $ref_hash_file_contents_per_run, 
		$full_pathname, $which_run_num_dir) = @_;

	print "\n\ninside read_file_contents\n";
	
	print "    full_pathname $full_pathname\n";
	print "        curr_file $curr_file\n";
	print "which_run_num_dir $which_run_num_dir\n";

#	print_Dumper("ref_hash_file_contents", $ref_hash_file_contents);
		
	my $num_run = 'bogus';
	
#	for my $line (@{ $diff->[0] }) {
#		$num_run =  "$line->[1]\n";
#	}
#	chomp($num_run);
#	print "num_run $num_run\n";
	
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
				
				process_this_cluster( \@array_sibling_genes_curr_cluster, $which_run_num_dir, $ref_hash_file_contents);
				
#				$ref_hash_file_contents_per_run->{"${which_run_num_dir}:${num_cluster}"} = \@array_sibling_genes_curr_cluster;
				
				my $ref_hash_inner_num_cluster_2_ref_array_genes;
				if (exists $ref_hash_file_contents_per_run->{$which_run_num_dir}) {
					
					$ref_hash_inner_num_cluster_2_ref_array_genes = $ref_hash_file_contents_per_run->{$which_run_num_dir};
				}
				push @{$ref_hash_inner_num_cluster_2_ref_array_genes->{$num_cluster}}, @array_sibling_genes_curr_cluster;
				$ref_hash_file_contents_per_run->{$which_run_num_dir} = $ref_hash_inner_num_cluster_2_ref_array_genes;
				
				print "which_run_num_dir $which_run_num_dir  num_cluster $num_cluster  @array_sibling_genes_curr_cluster\n";
				
				#  ---- post current cluster ----  #
				
				@array_sibling_genes_curr_cluster = ();  #  initialize array to emptiness
				
				$prev_num_cluster = $num_cluster;
			}
	
			print " $stubCRA  $stub_curr_cluster  $num_cluster  $KP_id  $FB_id\n";
				
			push @array_sibling_genes_curr_cluster, $FB_id;
		}	
	}
	process_this_cluster( \@array_sibling_genes_curr_cluster, $num_run, $ref_hash_file_contents);
	
}  #  read_file_contents

sub process_this_cluster {
	
	my ( $ref_array_sibling_genes_curr_cluster, $num_run, $ref_hash_file_contents) = @_;
	
	print "inside process with num_run $num_run and no more \n";
	
	print_Dumper("ref_array_sibling_genes_curr_cluster",  $ref_array_sibling_genes_curr_cluster);
	
	foreach my $curr_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
		
			print "curr_FB_id $curr_FB_id\n";
				
			foreach my $curr_sibling_FB_id ( @{$ref_array_sibling_genes_curr_cluster}) {
				
				next if ( $curr_sibling_FB_id eq $curr_FB_id);  #  skip if seeing self

				$ref_hash_file_contents->{$curr_FB_id}->{$curr_sibling_FB_id}++;
			}
	}
		
}  #  process_this_cluster



sub print_Dumper {
	
	my ( $label, $stuff_to_Dump ) = @_;
	
	print "label $label\n";
	print Dumper $stuff_to_Dump;
}

sub print_it {
	
	my ( $stuff_to_print ) = @_;
	
	if ( $send_print_to_file eq $YES ) {
		
		print $stuff_to_print;  #  also send to standard out as well as to file
		print $fh_standard_out $stuff_to_print;
		
	} else {

		print $stuff_to_print;
	}
}


1;  #  loaded OK