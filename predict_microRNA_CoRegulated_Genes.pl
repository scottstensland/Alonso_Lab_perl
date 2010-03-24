#!/usr/bin/perl -w

#  hello Kristiania

use strict;
use warnings;
use Data::Dumper;
use IO::File;
use Algorithm::Cluster qw/kcluster  clusterdistance/;
use Bio::SeqIO;
# use PDL;
#  use PDL::Graphics::PGPLOT; # http://pdl.sourceforge.net/PDLdocs/Graphics/PGPLOT/Window.html
# use Text::Levenshtein qw(distance fastdistance); # stens todo replace with compiled c version LevenshteinXS
use Algorithm::HowSimilar qw(compare);
use Benchmark;
use Time::HiRes qw(usleep ualarm gettimeofday tv_interval);

use Alonso_Lab::BioMotifDiscovery;
use Alonso_Lab::BioGraphics;
use Alonso_Lab::BioIO;
use Alonso_Lab::BioStrategy;
use Alonso_Lab::BioUtility;
use Alonso_Lab::BioPostProcess;

my $YES = 'YES';
my $NO  = 'NO';

$| = 1;  #  turn buffering off

# ------------    filter  <-- points of important data limiting logic - search for word: filter to see these

my $windows_dir = "C:/Users/Scott/Documents/data";
my $linux_dir = "/home/scott/Documents/data";
my $mac_OSX_dir = "/Users/scottstensland/Documents/data";

my $data_dir = "IGNORE_DIR";
my $to_plot_or_not = $NO;
my $vary_cluster_windowing = $NO;  #  flag whether we hardcode or dynamically calculate number of clusters  RRRRRRRRR
my $do_microRNA_sequence_match = $NO;
my $send_print_to_file = $NO;
#my $vary_timepoints_in_window = "vary_timepoints_in_window";  #  cluster optimization technique
my $percentage_cutoff_best_clusters = 0.25;	#  will keep genes from best X % of clusters as determined 
											#  by given fitness algo for further consideration - dumped to file
my $mask_value = 1;    #  until necesary just use value 1 for all points
my $value_weight_per_this_timepoint = 1;    # until necessary use value 1
	                                        # as weight for each time point
my $pad_len = 4;  #  determines total num chars in output string - used to in routine to pad 0 to left of given num

print "OS is saying $^O\n";

if ($^O eq "MSWin32" ) { $data_dir = $windows_dir;

} elsif ($^O eq "linux"  ) { $data_dir = $linux_dir;
	
} elsif ($^O eq "darwin" ) { $data_dir = $mac_OSX_dir;	} else {

	die "ERROR - did not recognize your OS $^O so cannot define data_dir\n";
}

my $num_clusers         = 4; # choose number of clusters
my $number_of_clustering_passes = 1; # typical value is 100 - must determine min good value
									 # for speed in testing use small value like 2
my $length_KP_gene_id = 11;  #  typical value would be   KP08680-316  if not seeing this 11 char string error out
my $length_FlyBase_id = 11;  #  typical value would be   FBgn0051188  if not seeing this 11 char string error out
my $flag_filter_dupe_mappings_KP_into_FB = $YES;  #  if YES then ignore time series KP curves which are duplicate
												  #  mappings from KP into same FlyBase gene identifer

#            <><><>  input files <><><>            #

#  my $gene_id_conversion_file = "$data_dir/GPL4455_id_conversion.txt";
my $gene_id_conversion_file = "$data_dir/GPL4455_id_conversion_no_M_control.txt";

# this time series dataset is taken from
#
# Identification of tightly regulated groups of genes during Drosophila melanogaster embryogenesis
# Sean D Hooper, Stephanie Bou, Roland Krause, Lars J Jensen, Christopher E Mason, Murad Ghanim,
# Kevin P White, Eileen EM Furlong, Peer Bork
# Molecular Systems Biology 3, (16 January 2007) doi:10.1038/msb4100112
# http://www.nature.com/msb/journal/v3/n1/full/msb4100112.html
#  time series microarray Hooper & White 2007 http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4455
my $ncbi_dataset = "$data_dir/GSE6186_series_matrix.txt";# time series GEO dataset

my $num_data_lines      = 100; # value 0 is unlimited, else value limits num rows of input data (genes)
my $init_num_timepoints = 10; # use value of 0 for unlimited time points, else value limits num time points

my $three_prime_UTR_sequence_file = "$data_dir/dmel-all-three_prime_UTR-r5.24.fasta";
my $num_3_prime_UTR_sequence_lines = 0; # use value of 0 to pull in all 3' UTR rows else value limits num rows to read

my $microRNA_file = "$data_dir/dmel-all-miRNA-r5.24.fasta";
my $num_microRNA_lines = 20; # use value of 0 to pull in all microRNA rows else value limits num rows to read

#  taken from http://www.russell.embl.de/nob/miRNA/StarkBrennecke2005_All.csv
#             which hangs off from  http://www.mirbase.org/help/targets.shtml
my $microRNA_target_site_prediction_matches_file_input = "$data_dir/20100301_mirbase_org_www_russel_embl_de_StarkBrennecke2005_All.csv";

#my $gene_pairs_counts_across_multi_runs = "/home/scott/Documents/data/alonso_lab/synth_op_dist2_copy.txt";

#            <><><>  output files <><><>            #

my $input_dir_prefix = "$data_dir/alonso_lab/play_OSX"; #  used to name set of output data files for plotting in pyxplot 8-)))

my $curr_num_run = pad_number(1, $pad_len);  #  initialize multiple run number indicator - used in file names
my $max_num_run  = 5;  #  if $multi_run == $YES then do this number of runs
my $bootstrapping_percentage_criteria = 0.80;


#my $run_name_dir = "$data_dir/alonso_lab/repeat_v2_fast_35_50_144_run10";  #  used to name set of output data files for plotting & post proc
#my $run_name_dir = "$data_dir/alonso_lab/play_20100321";  #  used to name set of output data files for plotting & post proc
my $run_name_dir = "${input_dir_prefix}/${curr_num_run}";  #  used to name set of output data files for plotting & post proc

my $standard_out_filename = "$run_name_dir/standard_out.txt";

#my $match_fitness_btw_all_genes_N_all_microRNA_self_output    = "$run_name_dir/match_fitness_btw_all_genes_N_all_microRNA_self.txt";
my $match_fitness_btw_all_genes_N_all_microRNA_miRBase_output = "$run_name_dir/match_fitness_btw_all_genes_N_all_microRNA_miRBase.txt";

my $genes_in_best_clusters_fitness_best_curve_output = "genes_in_best_clusters_fitness_best_curve.txt";

my $histogram_filename_output = "histogram_output.dat";
my $histogram_filename_pyxplot_output = "$run_name_dir/histogram_filename_pyxplot_output.ppl";

my $gene_pairs_counts_across_multi_runs = "${input_dir_prefix}/post_process.txt";  #  post processing aggregation intermediate file

my $multi_run = $NO;	#  flag indicating to do multiple runs or not - useful for multiple $RUN_NORMAL
						#  this auto incrments 
#  cluster output files prefix
my $cluster_parm_label = "cparml";  #  num_clusers / num_timepoints  /  start_timepoint  / stop_timepoint

my $RUN_NORMAL          = $NO;  #  read time series input file and prepare datastructures for file outputs
my $RUN_POST_PROCESSING = $NO;  #  aggregate output generated from RUN_NORMAL to prepare for RUN_DE_NOVO
my $RUN_DE_NOVO         = $NO;  #  read precalculated gene pair counts then de novo motif discovery

#     here is repeat where dir looks like    repeat_fast_35_50_144_run01 <--> run10
my @preset_testing_fast   = ( 85, 30, 1440,  0, 0, 0,  $NO,  $NO,  $NO, $YES,   $YES, $YES, $YES, $NO );
#my @preset_testing_fast   = ( 35, 50, 144,  0, 0, 0,  $NO,  $NO,  $YES, $YES,  $NO, $NO, $NO, $NO );

my @preset_testing_medium = ( 81, 50, 0,    0, 0, 0,   $NO,  $NO,  $NO, $YES,   $NO, $NO, $NO, $NO );  # 144,   14^2 = 196
my @preset_semi_prod	  = ( 4, 1, 0,      2, 0, 0,   $YES, $NO,  $NO, $NO,    $NO, $NO, $NO, $NO );

(	$num_clusers,			$number_of_clustering_passes,	  $num_data_lines,
	$init_num_timepoints,	$num_3_prime_UTR_sequence_lines,  $num_microRNA_lines, 
	$vary_cluster_windowing,   $do_microRNA_sequence_match,   $to_plot_or_not,  $send_print_to_file,
	$multi_run,    $RUN_NORMAL, $RUN_POST_PROCESSING, $RUN_DE_NOVO )

		= @preset_testing_fast;

print "input_dir_prefix $input_dir_prefix\n";

if (! -d ${input_dir_prefix} ) {

	mkdir(${input_dir_prefix}) or die "ERROR - cannot create input_dir_prefix ${input_dir_prefix} : $!\n";
}

if (! -d ${run_name_dir} ) {

	mkdir(${run_name_dir})  or die "ERROR - cannot create run_name_dir ${run_name_dir} : $!\n";
}
		
my $fh_standard_out = IO::File->new( "> $standard_out_filename" )
			or die "ERROR - cannot open file standard_out_filename " .
					"$standard_out_filename : $!\n";
		
sub print_it {
	
	my ( $stuff_to_print ) = @_;
	
	if ( $send_print_to_file eq $YES ) {
		
		print $stuff_to_print;  #  also send to standard out as well as to file
		print $fh_standard_out $stuff_to_print;
		
	} else {

		print $stuff_to_print;
	}
}

print_it "num_clusers $num_clusers\n";
print_it "number_of_clustering_passes $number_of_clustering_passes\n";
print_it "num_data_lines $num_data_lines\n";
print_it "init_num_timepoints $init_num_timepoints\n";
print_it "num_3_prime_UTR_sequence_lines $num_3_prime_UTR_sequence_lines\n";
print_it "num_microRNA_lines $num_microRNA_lines\n";


sub print_Dumper {
	
	my ( $label, $stuff_to_Dump ) = @_;
	
	print "label $label\n";
	print Dumper $stuff_to_Dump;
}

# die "did direct setting work\n";

&main(); # entry point of entire script - from start to finish

# <><><>  <><><>  <><><>  end of logic flow <><><>  <><><>  <><><>

sub main {
	
	my $debug = $YES;
	
	print_it "start of computation\n";
	
	my $motif_discovery = new BioMotifDiscovery();
	my $bioUtility      = new BioUtility();
	my $bioIO    		= new BioIO();
	my $bioStrategy     = new BioStrategy();
	my $bioPostProcess  = new BioPostProcess();
	
	$bioIO->set_standard_out_fh($fh_standard_out);
	
	$motif_discovery->set_send_print_to_file($send_print_to_file);
	     $bioUtility->set_send_print_to_file($send_print_to_file);
	          $bioIO->set_send_print_to_file($send_print_to_file);
	    $bioStrategy->set_send_print_to_file($send_print_to_file);
	 $bioPostProcess->set_send_print_to_file($send_print_to_file);	
	
	my $time0 = new Benchmark;
	
	my %hash_UTR_motif_count_2_ref_array_motifs_success;
	my %hash_UTR_motif_count_2_ref_array_motifs_mismatch;
	
	my %hash_UTR_segment_freq_success;  #  count of each sequence motif which matches btw parent N child
	my %hash_UTR_segment_freq_mismatch; #  count of such mismatches, counts each so double qty of entries
			
	my $ref_hash_gene_KB_to_3_prime_UTR_s; # hash to store 3' UTR sequence(s) for each gene identifer
	my $ref_hash_KP_to_FB; # hash to map between KP ids found in time series and FB found in sequence file
	my $ref_hash_FlyBase_to_array_of_KP;  #  inverse of ref_hash_KP_to_FB which maintains dupe maps KP to FB using array
	
	my $ref_hash_all_microRNA;  #  hash of all microRNA by gene id
	
	# FB to TGAATTAACCATACTATACAACTATATGT ...
	$ref_hash_gene_KB_to_3_prime_UTR_s = 
		$bioIO->parse_sequence_file_fasta_format
			($three_prime_UTR_sequence_file, $num_3_prime_UTR_sequence_lines, "parent");
	
#	print_Dumper ("ref_hash_gene_KB_to_3_prime_UTR_s", $ref_hash_gene_KB_to_3_prime_UTR_s );
	
	
	$ref_hash_all_microRNA = 
		$bioIO->parse_sequence_file_fasta_format($microRNA_file, $num_microRNA_lines, "name");
	
#	print_it "to show ref_hash_all_microRNA\n";
#	print_Dumper ( $ref_hash_all_microRNA );
	
#	print_it "end of show post read sequence files\n";
#	return;
	
	#  //////////////  \\\\\\\\\\\\\\           ////////////// \\\\\\\\\\\\\	
		
	( $ref_hash_KP_to_FB, $ref_hash_FlyBase_to_array_of_KP) = 
			$bioIO->parse_conversion_file
				($gene_id_conversion_file, $length_FlyBase_id,
				$flag_filter_dupe_mappings_KP_into_FB); # KP=KP01054-747	FBgn=FBgn0002576	CG=CG1689
	
	my $how_many_KP_ids = scalar(keys %{$ref_hash_KP_to_FB});
	my $how_many_FlyBase_ids = scalar(keys %{$ref_hash_FlyBase_to_array_of_KP});
	
	print "how_many_KP_ids $how_many_KP_ids    \nhow_many_FlyBase_ids $how_many_FlyBase_ids\n";
	
#	print_Dumper ( "ref_hash_KP_to_FB", $ref_hash_KP_to_FB );
#	print_Dumper ( "ref_hash_FlyBase_to_array_of_KP", $ref_hash_FlyBase_to_array_of_KP );
	
#	print_it "\nend of parse_conversion_file\n";
#	return;

	if ($RUN_NORMAL eq $YES) {
		
		    my $continue_running = $YES;
		while ($continue_running eq $YES) {
			
			if (! -e ${run_name_dir} ) {
		
				mkdir(${run_name_dir}) or die "ERROR - cannot create run_name_dir ${run_name_dir} : $!\n";
			}
		
			close $fh_standard_out;
			
			$fh_standard_out = IO::File->new( "> $standard_out_filename" )
				or die "ERROR - cannot open file standard_out_filename " .
					"$standard_out_filename : $!\n";
		
	        $motif_discovery->set_standard_out_fh($fh_standard_out);
	             $bioUtility->set_standard_out_fh($fh_standard_out);
	                  $bioIO->set_standard_out_fh($fh_standard_out);
	            $bioStrategy->set_standard_out_fh($fh_standard_out);
	         $bioPostProcess->set_standard_out_fh($fh_standard_out);
		
			print "run_name_dir $run_name_dir\n";
					
			my (	$ref_array_time_series, $mask1, $weight1, 
					$num_timepoints, $ref_array_KP_gene_identifers ) = 
				$bioIO->parse_NCBI_timeseries_dataset(
					$ncbi_dataset, $init_num_timepoints, $ref_hash_KP_to_FB,
					$ref_hash_gene_KB_to_3_prime_UTR_s, $num_data_lines, $length_KP_gene_id,
					$mask_value, $value_weight_per_this_timepoint );
		            
			
			$bioStrategy->prosecute_strategy(	
					$num_clusers, $ref_array_time_series, $ref_array_KP_gene_identifers,
					$number_of_clustering_passes, $num_timepoints, $ref_hash_KP_to_FB,
					$ref_hash_gene_KB_to_3_prime_UTR_s, $ref_hash_all_microRNA,
					$vary_cluster_windowing, $do_microRNA_sequence_match, $run_name_dir,
					$debug, $cluster_parm_label, $genes_in_best_clusters_fitness_best_curve_output,
					$value_weight_per_this_timepoint, $mask_value,
					$percentage_cutoff_best_clusters, $to_plot_or_not,
					$microRNA_target_site_prediction_matches_file_input,
					$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output,
					$curr_num_run, $max_num_run,
			);
			
			if ( $vary_cluster_windowing eq $YES ) {
				
				#  RRRR need to revisit this to handle set of files each run 
				#       maybe put inside prosecute_strategy
			
	#			plot_histogram($genes_in_best_clusters_fitness_best_curve_output, $run_name_dir,
	#							$histogram_filename_output, $histogram_filename_pyxplot_output,
	#							$cluster_parm_label);
			}

			if ($multi_run eq $YES && $curr_num_run < $max_num_run) {
				
				$curr_num_run++;				
				$curr_num_run = pad_number($curr_num_run, $pad_len);

				print "curr_num_run $curr_num_run\n";
								
				$run_name_dir = "${input_dir_prefix}/${curr_num_run}";
				$standard_out_filename = "$run_name_dir/standard_out.txt";
				$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output = "$run_name_dir/match_fitness_btw_all_genes_N_all_microRNA_miRBase.txt";
				$histogram_filename_pyxplot_output = "$run_name_dir/histogram_filename_pyxplot_output.ppl";
	
			} else {
				
				$continue_running = $NO;
			}
		}
	}
	
	if ($RUN_POST_PROCESSING eq $YES) {
		
		print_it "about to call post processing\n";
		
		my $fh_gene_pairs_count_file = IO::File->new( "> $gene_pairs_counts_across_multi_runs" )
					or die "ERROR - cannot open file gene_pairs_count_file " .
							"$gene_pairs_counts_across_multi_runs : $!\n";
		
		$bioPostProcess->post_process($input_dir_prefix, $fh_gene_pairs_count_file);
		
		print "\ngene_pairs_counts_across_multi_runs \n$gene_pairs_counts_across_multi_runs\n";
	} 
	
	if ($RUN_DE_NOVO eq $YES) {
		
		print "gene_pairs_counts_across_multi_runs $gene_pairs_counts_across_multi_runs\n";
		
#		print_Dumper("ref_hash_gene_KB_to_3_prime_UTR_s", $ref_hash_gene_KB_to_3_prime_UTR_s);
		
		$motif_discovery->de_novo_motif_discovery($gene_pairs_counts_across_multi_runs,
								$ref_hash_gene_KB_to_3_prime_UTR_s,
								\%hash_UTR_segment_freq_success,
								\%hash_UTR_segment_freq_mismatch);

		print_Dumper("hash_UTR_segment_freq_success", \%hash_UTR_segment_freq_success);
		print_Dumper("hash_UTR_segment_freq_mismatch", \%hash_UTR_segment_freq_mismatch);

		print "number of unique success  motif matches: ". scalar(keys %hash_UTR_segment_freq_success) ."\n";
		print "number of unique mismatch motif matches: ". scalar(keys %hash_UTR_segment_freq_mismatch) ."\n";
	
		$bioUtility->flip_hash_into_hash_values_of_ref_arrays_of_keys(
				\%hash_UTR_segment_freq_success, 
				\%hash_UTR_motif_count_2_ref_array_motifs_success);
		
#		print_Dumper("hash_UTR_segment_freq_success", \%hash_UTR_segment_freq_success);
#		print_Dumper("hash_UTR_motif_count_2_ref_array_motifs_success", \%hash_UTR_motif_count_2_ref_array_motifs_success);

		my $total_motifs;
		
		$total_motifs = BioMotifDiscovery::sum_total_motifs_found(\%hash_UTR_motif_count_2_ref_array_motifs_success);
		print "\n\ntotal_motifs success $total_motifs\n";

		$bioUtility->flip_hash_into_hash_values_of_ref_arrays_of_keys(
				\%hash_UTR_segment_freq_mismatch, 
				\%hash_UTR_motif_count_2_ref_array_motifs_mismatch);

#		print_Dumper("hash_UTR_motif_count_2_ref_array_motifs_mismatch", \%hash_UTR_motif_count_2_ref_array_motifs_mismatch);

		$total_motifs = BioMotifDiscovery::sum_total_motifs_found(\%hash_UTR_motif_count_2_ref_array_motifs_mismatch);
		print "\n\ntotal_motifs mismatch $total_motifs\n";
	}
	
	print_it "\n\nrun_name_dir $run_name_dir\n";
	
	my $time1 = new Benchmark;
	my $time_synth_output = timestr(timediff($time1,$time0));
	
	my $interval_in_secs = tv_interval($time0, $time1);
	
	print "interval_in_secs $interval_in_secs\n";
		
#	print_it "\nScript Benchmark: ".timestr(timediff($time1,$time0));
	print_it "\nScript Benchmark: " . $time_synth_output;
	
	print_it "\n\nEnd of Computation ... completed in $interval_in_secs seconds\n";
	
}  #  end of main

sub pad_number {
	
	my ($given_num_to_pad_with_zeros_to_left, $pad_len) = @_;
	
	my  $padded = sprintf("%0${pad_len}d", $given_num_to_pad_with_zeros_to_left);
#	print "padded $padded\n";
	
	return $padded;	
}

#  -------------
