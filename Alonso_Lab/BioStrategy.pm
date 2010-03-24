
package BioStrategy;

use warnings;
use Data::Dumper;
use strict;
use Exporter qw(import);

use POSIX qw/ceil/;

# use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(prosecute_strategy flip_hash_into_hash_values_of_ref_arrays_of_keys);
our %EXPORT_TAGS = ( DEFAULT => [qw(&prosecute_strategy)],
                 Both    => [qw(&prosecute_strategy &flip_hash_into_hash_values_of_ref_arrays_of_keys)]);

my $YES = 'YES';
my $NO  = 'NO';

my $fh_standard_out;  #  standard out is send into this supplied file handle
my $send_print_to_file = $NO;  #  flag to determine whether to send standard out to output file

#   $| = 1;  #  turn buffering off

    sub new {

    	my $proto = shift;
        my $class = ref($proto) || $proto;
    	
    	my $self = {
    		
    		_fh_standard_out => undef,
    		
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
	

sub print_it {
	
	my ( $stuff_to_print ) = @_;
	
	if ( $send_print_to_file eq $YES ) {
		
		print $stuff_to_print;  #  also send to standard out as well as to file
		print $fh_standard_out $stuff_to_print;
		
	} else {

		print $stuff_to_print;
	}
}

sub print_Dumper {
	
	my ( $label, $stuff_to_Dump ) = @_;
	
	print "label $label\n";
	print Dumper $stuff_to_Dump;
}
sub prosecute_strategy {
	
	my (	$self, $num_clusers, $ref_array_time_series_canonical, $ref_array_KP_gene_identifers,
			$number_of_clustering_passes, $num_timepoints, $ref_hash_KP_to_FB,
			$ref_hash_gene_KB_to_3_prime_UTR_s, $ref_hash_all_microRNA,
			$vary_cluster_windowing, $do_microRNA_sequence_match, $run_name_dir, 
			$debug, $cluster_parm_label, $genes_in_best_clusters_fitness_best_curve_output,
			$value_weight_per_this_timepoint, $mask_value, $percentage_cutoff_best_clusters,
			$to_plot_or_not, $microRNA_target_site_prediction_matches_file_input,
			$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output,
			$curr_num_run, $max_num_run,
	) = @_;
	
	$debug = $YES;
	
#	print_Dumper ( "ref_array_time_series_canonical",  $ref_array_time_series_canonical );
	
	my $speed_up_factor = 1;  #  default is 1 - see below for speed up if desired
	my $curr_num_timepoints;
	my $minimum_num_timepoints = $num_timepoints - 1;  #  only used for normal mode; is overwritten o/w
	my $maximum_num_timepoints = $num_timepoints;
		
	# if $NO then use given var $num_clusers,  else if $YES then calculate $num_clusers
	my $flag_alter_num_clusers_based_on_num_timepoints_in_time_series_window = $NO; 
		
#	my $biostrategy     = new BioStrategy();
		
	if ( $vary_cluster_windowing eq $YES ) {
		
		#  if $NO then we do one run of all available time points in time series
		#  and do cluster on all time points all at once 
		
		#  else if $YES then do clustering many times --- repeatedly cluster inside
		#  a double nested pair of loops where outer loop varies number of time points
		#  to look at and inner loop varies which time point we start from
		#
		#  For example:  outer loop is variable below called: $minimum_num_timepoints
		#                inner loop is variable below called: X - put troubleshooting flag into below
		#                                                     so we can see these outer and inner loop values printed
		
	   $speed_up_factor = 10;  #  no speed up == 1,  speed up factor determines size of increment of
							  #                     window size and left-to-right movement increment, useful
							  #                     when iterating across hundreds of runs in below loops
							  
		$minimum_num_timepoints = 5;
		$flag_alter_num_clusers_based_on_num_timepoints_in_time_series_window = $YES; 
	}
	
	while ($minimum_num_timepoints < $maximum_num_timepoints) {
		
		if ($flag_alter_num_clusers_based_on_num_timepoints_in_time_series_window eq $YES) {
			
			# in this mode where we dynamically calc num of clusters
			#  will vary $num_clusers by formula 
			
			my $minimum_number_of_clusters = 25;  # set baseline minimum number of clusters.
												  # Especially important for low values 
												  # of $minimum_num_timepoints
			
			#  number of clusters is twice size of TS window
			$num_clusers = $minimum_number_of_clusters + $minimum_num_timepoints * 2.0;
		}
		
		my $curr_start_timepoint = 0;
		my $curr_stop_timepoint = $curr_start_timepoint + $minimum_num_timepoints;
		while ($curr_stop_timepoint < $maximum_num_timepoints) {
				
			print "MMM num_clusers $num_clusers num_timepoints $minimum_num_timepoints start $curr_start_timepoint  stop $curr_stop_timepoint\n";
			
			my $curr_run_cluster_parms = "${cluster_parm_label}_${num_clusers}_${minimum_num_timepoints}_" .
						"${curr_start_timepoint}_${curr_stop_timepoint}";
						
						
			my $cluster_op_file = "${run_name_dir}/${curr_run_cluster_parms}_${genes_in_best_clusters_fitness_best_curve_output}";
						
						
			
			my @array_weight_per_timepoint;				
			my @array_mask;
			my @array_time_series;
		
			foreach my $ref_array_one_row_data ( @{$ref_array_time_series_canonical}) {

			# print_it ( $curr_gene );
			
				# my $curr_gene_time_series = $ref_array_time_series_canonical->[$curr_gene];
			
				@array_weight_per_timepoint = ();  #  empty for current window size
				my @array_one_row_mask;
				my @array_one_row_data;
				my $curr_timepoint;
									
#				print "minimum_num_timepoints $minimum_num_timepoints " if ($debug eq $YES);
				
				for ($curr_timepoint = $curr_start_timepoint; $curr_timepoint <= $curr_stop_timepoint; $curr_timepoint++) {

					my $one_timepoint_Y_value = $ref_array_one_row_data->[$curr_timepoint];
					
					push @array_weight_per_timepoint, $value_weight_per_this_timepoint;
					push @array_one_row_mask, $mask_value;
					push @array_one_row_data, $one_timepoint_Y_value;
					
#					print "time $curr_timepoint " if ($debug eq $YES);
				}
#				print "\n" if ($debug eq $YES);
				
				
				push @array_mask,  \@array_one_row_mask;
				push @array_time_series, \@array_one_row_data;	
			}

		# <><><> below is called for each window size of number of cluster timepoints <><><> #
			
		#	print_it "ref_array_time_series\n";
		#	print_Dumper ( $ref_array_time_series );
		
			my $ref_hash_cluster_fitness_2_cluster;  #  measure of euclidean distance fitness of all clusters
			my $have_we_already_loaded_gene_FB_2_microRNA = $NO;
			my $ref_hash_gene_FB_id_2_microRNA_match_scores;
			
			print_it "\nrun $curr_num_run of $max_num_run    calculating clusters ...\n";
			my ( $ref_array_mapping_gene_to_cluster, $gen_error, 
				$gen_found, $ref_hash_mapping_cluster_to_genes ) =
					&do_Algorithm_Cluster( $num_clusers, \@array_time_series, \@array_mask, 
										\@array_weight_per_timepoint, $ref_array_KP_gene_identifers, 
										$number_of_clustering_passes, $debug );
			
			# calculate tightness per cluster
			
			print "num_timepoints $minimum_num_timepoints start $curr_start_timepoint  stop $curr_stop_timepoint\n";
						
			$ref_hash_cluster_fitness_2_cluster = 
				&identify_cluster_fitness_based_on_curve_shape(
					$num_clusers, $ref_array_mapping_gene_to_cluster, \@array_time_series, 
					$num_timepoints, $ref_array_KP_gene_identifers, $run_name_dir,
					$genes_in_best_clusters_fitness_best_curve_output,
					$ref_hash_mapping_cluster_to_genes, $percentage_cutoff_best_clusters,
					$minimum_num_timepoints, $curr_start_timepoint, $curr_stop_timepoint,
					$cluster_op_file, $ref_hash_KP_to_FB );						
		
			BioGraphics::do_plot($to_plot_or_not, $num_clusers, $ref_array_mapping_gene_to_cluster,
					\@array_time_series, $num_timepoints, $ref_array_KP_gene_identifers,
					$run_name_dir );

			if ( $do_microRNA_sequence_match eq $YES ) {
					  
				if ( $have_we_already_loaded_gene_FB_2_microRNA eq $NO ) {
					
					$ref_hash_gene_FB_id_2_microRNA_match_scores =
						&do_3_prime_UTR_motif_pattern_match_miRBase(
								$microRNA_target_site_prediction_matches_file_input,
								$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output );
		
					$have_we_already_loaded_gene_FB_2_microRNA = $YES;
				}
				
				&judge_cluster_fitness_based_on_microRNA_match_scoring(
							$ref_array_mapping_gene_to_cluster, 
							$ref_hash_KP_to_FB,
							$ref_array_KP_gene_identifers, 
							$ref_hash_gene_KB_to_3_prime_UTR_s,
							$ref_hash_mapping_cluster_to_genes, 
							$num_clusers, 
							$ref_hash_all_microRNA,
							$ref_hash_gene_FB_id_2_microRNA_match_scores );
			}
						
#			execute_de_novo_motif_discovery($ref_hash_cluster_fitness_2_cluster);
			
			$curr_start_timepoint += $speed_up_factor;
			$curr_stop_timepoint = $curr_start_timepoint + $minimum_num_timepoints;
		}
		$minimum_num_timepoints += $speed_up_factor;  
	}
	#  vary size of timeseries window 

}  #  prosecute_strategy



sub judge_cluster_fitness_based_on_microRNA_match_scoring {
	
	my (
		$ref_array_mapping_gene_to_cluster, 
		$ref_hash_KP_to_FB,
		$ref_array_KP_gene_identifers, 
		$ref_hash_gene_KB_to_3_prime_UTR_s,
		$ref_hash_mapping_cluster_to_genes, 
		$num_clusers, 
		$ref_hash_all_microRNA,
		$ref_hash_gene_FB_id_2_microRNA_match_scores
	) = @_;
	
	print_it "inside judge_cluster_fitness_based_on_microRNA_match_scoring\n";

	# print_it "about to dump ref_hash_all_microRNA\n";
	# print_Dumper ( $ref_hash_all_microRNA );
	
	print_it "about to dump ref_hash_mapping_cluster_to_genes\n";
	print_Dumper ( "ref_hash_mapping_cluster_to_genes", $ref_hash_mapping_cluster_to_genes );
	
	print_it "about to see ref_hash_gene_FB_id_2_microRNA_match_scores\n";
	print_Dumper ( "ref_hash_gene_FB_id_2_microRNA_match_scores", $ref_hash_gene_FB_id_2_microRNA_match_scores );
	
	my $count_num_FlyBase_ids_seen = 0;
	my %hash_unique_FlyBase_id_seen;
	my %hash_cluster_fitness;  #  storage for info supporting fitness per cluster - hash key is cluster number
	
	foreach my $curr_cluster (sort { $a <=> $b } keys %{$ref_hash_mapping_cluster_to_genes}) {
		
		#  sort { $a <=> $b }  #  this is numeric sort (0 1 2 .. 10 11) not alphabetic (0 1 10 11 2 3 4 ...)
		
		print_it "\n\n\n\n<><><> new cluster " .(1 + $curr_cluster). " of $num_clusers <><><>\n\n\n\n";
		
		my @array_scalar_fitness_related_to_this_cluster;
		my %hash_microRNA_to_FB_ID_this_cluster;
		
#		my %hash_all_3_prime_UTRs_for_all_genes_in_this_cluster;
		
		foreach my $curr_gene (@{$ref_hash_mapping_cluster_to_genes->{$curr_cluster}}) {

			my $KP_id = $ref_array_KP_gene_identifers->[$curr_gene];
			my $FlyBase_id = $ref_hash_KP_to_FB->{$KP_id};
			
			my %hash_microRNA_to_match_score;
			if (exists $ref_hash_gene_FB_id_2_microRNA_match_scores->{$FlyBase_id}) {
				
				%hash_microRNA_to_match_score = %{$ref_hash_gene_FB_id_2_microRNA_match_scores->{$FlyBase_id}};

#				print_it "OOOOOOOOOOKKKKKKKKKKKKKKK FlyBase_id $FlyBase_id  hash_microRNA_to_match_score\n";
#				print_Dumper ( \%hash_microRNA_to_match_score );
				
				while (my ($key_microRNA, $value_match_score) = each(%hash_microRNA_to_match_score)){
					
#					print_it "key_microRNA $key_microRNA    value_match_score $value_match_score\n";
					
					my $ref_array_FlyBase_id_N_match_score_this_microRNA_this_cluster;
					
					if (exists $hash_microRNA_to_FB_ID_this_cluster{$key_microRNA}) {
						
						$ref_array_FlyBase_id_N_match_score_this_microRNA_this_cluster = $hash_microRNA_to_FB_ID_this_cluster{$key_microRNA};
					}
					
					push @{$ref_array_FlyBase_id_N_match_score_this_microRNA_this_cluster}, "${FlyBase_id}:${value_match_score}";
					
					$hash_microRNA_to_FB_ID_this_cluster{$key_microRNA} = $ref_array_FlyBase_id_N_match_score_this_microRNA_this_cluster;
				}
				
				#  %hash_microRNA_to_FB_ID_this_cluster
				
				$hash_unique_FlyBase_id_seen{$FlyBase_id} = 1;
				$count_num_FlyBase_ids_seen++;

			} else {
				
				# print_it "NOTICE - failed to find FlyBase_id $FlyBase_id inside ref_hash_gene_FB_id_2_microRNA_match_scores\n";
			}			
			
#			print_it "judge curr_cluster $curr_cluster looking at KP_id $KP_id   FlyBase_id $FlyBase_id\n";	
		}
		
		print_it "\nhash_microRNA_to_FB_ID_this_cluster " .(1 + $curr_cluster). " of $num_clusers\n";
		print_Dumper ( "hash_microRNA_to_FB_ID_this_cluster", \%hash_microRNA_to_FB_ID_this_cluster );
		
	}  #  foreach $curr_cluster
	
	print_it "count unique FlyBase ids seen using hash " . scalar(keys %hash_unique_FlyBase_id_seen). "\n";
	print_it "count_num_FlyBase_ids_seen $count_num_FlyBase_ids_seen\n";
	
}  #  judge_cluster_fitness_based_on_microRNA_match_scoring

sub do_3_prime_UTR_motif_pattern_match_miRBase {
	
	my (	
	
#	$ref_array_KP_gene_identifers, 
#	$ref_hash_gene_KB_to_3_prime_UTR_s,
#			$ref_hash_mapping_cluster_to_genes, 
#			$num_clusers, 
#			$ref_hash_all_microRNA,
			$microRNA_target_site_prediction_matches_file_input,
			$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output
	) = @_;	
	
	my %hash_gene_FB_id_2_microRNA_match_scores;
	  			
	my $fh_microRNA_target_site_prediction_matches_file_input = 
		IO::File->new( "< $microRNA_target_site_prediction_matches_file_input")
	  	or die "ERROR - cannot open file microRNA_target_site_prediction_matches_file_input " .
	  			"$microRNA_target_site_prediction_matches_file_input : $!\n";
	
	my $fh_match_fitness_btw_all_genes_N_all_microRNA_miRBase_output = 
		IO::File->new(  ">> $match_fitness_btw_all_genes_N_all_microRNA_miRBase_output")
	  	or die "ERROR - cannot open file match_fitness_btw_all_genes_N_all_microRNA_miRBase_output " .
	  			"$match_fitness_btw_all_genes_N_all_microRNA_miRBase_output\n";

#	print_it "ref_array_mapping_gene_to_cluster\n";
#	print_Dumper ( $ref_array_mapping_gene_to_cluster );
	
#	print_it "ref_hash_mapping_cluster_to_genes\n";
#	print_Dumper ( $ref_hash_mapping_cluster_to_genes );
	
#	print_it "ref_array_KP_gene_identifers\n";
#	print_Dumper ( $ref_array_KP_gene_identifers );
	
	
#	print_it "ref_hash_KP_to_FB\n";
#	print_Dumper ( $ref_hash_KP_to_FB );
	
	
	#  print $fh_match_fitness_btw_all_genes_N_all_microRNA_miRBase_output "curr_cluster $curr_cluster $curr_FlyBase_id\t$curr_microRNA\t$curr_accuracy_percentage\n";

	#  my $line_count = 0;
	while (<$fh_microRNA_target_site_prediction_matches_file_input>) {

		# $curr_artificially_limit_num_timepoints = 0;
		# $_ =~ s/"//g;    #  remove the double quote as in "
		
		next if ( $_ =~ m/^\#/ );  #  skip row if starts with comment symbol #
		
		chomp;  #  strip off trailing newline char
		
		# if (/series_matrix_table_begin/)

		# miRNA, CG-ID, gene name, Flybase-ID, validated (Real) or predicted (PRED) UTR, score of best site, total score of all sites, number of sites

		my ( $microRNA, $CG_ID, $gene_name, $flybaseID, $validated_or_predicted,
			$score_of_best_site, $total_score_of_all_sites,
			$num_sites_with_a_match ) = split('\t');
				
		print_it "flybaseID $flybaseID  microRNA $microRNA CG_ID $CG_ID " .
				"total_score_of_all_sites $total_score_of_all_sites gene_name $gene_name " .
				"validated_or_predicted $validated_or_predicted score_of_best_site $score_of_best_site " .
				"num_sites_with_a_match $num_sites_with_a_match\n";

		$hash_gene_FB_id_2_microRNA_match_scores{$flybaseID}{$microRNA} = $total_score_of_all_sites;
	}
	
#	seeing 4219 unique flybaseID in hash_gene_FB_id_2_microRNA_match_scores
#	my $num_FB_ids = scalar(keys %hash_gene_FB_id_2_microRNA_match_scores);
#	print_it "num_FB_ids $num_FB_ids\n";
	
	print_it "hash_gene_FB_id_2_microRNA_match_scores\n";
	print_Dumper ( "hash_gene_FB_id_2_microRNA_match_scores", \%hash_gene_FB_id_2_microRNA_match_scores );
	
	return ( \%hash_gene_FB_id_2_microRNA_match_scores );
		
}  #  do_3_prime_UTR_motif_pattern_match_miRBase
				
sub identify_cluster_fitness_based_on_curve_shape {

	my ( $num_clusers, $ref_array_mapping_gene_to_cluster, $ref_array_time_series, 
	     $num_timepoints, $ref_array_KP_gene_identifers, $run_name_dir,
	     $genes_in_best_clusters_fitness_best_curve_output,
	     $ref_hash_mapping_cluster_to_genes, $percentage_cutoff_best_clusters,
	     $minimum_num_timepoints, $curr_start_timepoint, $curr_stop_timepoint,
	     $cluster_op_file, $ref_hash_KP_to_FB ) = @_;
	     
#	print_Dumper ( "ref_array_mapping_gene_to_cluster", $ref_array_mapping_gene_to_cluster);
#	print_Dumper ( "ref_array_KP_gene_identifers", $ref_array_KP_gene_identifers);

	my $ref_hash_cluster_fitness_2_cluster;

	my $num_per_side = ceil( sqrt($num_clusers) );

	my $num_x_panel;
	my $num_y_panel;
#	my $curr_cluster;
	my $curr_output_filename = "";
	my $curr_gene = 0;


#ref_hash_mapping_cluster_to_genes


	foreach my $curr_cluster (sort {$a<=>$b} keys %{$ref_hash_mapping_cluster_to_genes}) {

#		print_it "\n\n\n\n<><><> new cluster " .(1 + $curr_cluster). " of $num_clusers <><><>\n\n\n\n";
		
		# Here is new data structure which will store all Y values for a given X time point
		# across all X time points
		
		my %hash_timepoints_2_gene_values_this_cluster;

		my @array_all_genes_this_cluster = @{$ref_hash_mapping_cluster_to_genes->{$curr_cluster}};
		
		if ( scalar(@array_all_genes_this_cluster) == 1) {  # filter = skip over single gene clusters
		
			my $single_gene_in_cluster = $array_all_genes_this_cluster[0];
			my $KP_id = $ref_array_KP_gene_identifers->[$single_gene_in_cluster];
			my $KB_id = $ref_hash_KP_to_FB->{$KP_id};
			
			print_it "NOTICE - filter on single gene @array_all_genes_this_cluster in cluster $curr_cluster " .
				"KP_id $KP_id    KB_id $KB_id\n";
				
			next; 	#  skip this cluster as it contains only a single gene hence is cannot be
					#  judged in this subroutine as it values are the mean so error is 0 and meaningless
		}

		my $count_num_genes_this_cluster = 0;
		foreach my $curr_gene ( @array_all_genes_this_cluster ) {

			my $curr_gene_time_series = $ref_array_time_series->[$curr_gene];
					
			# print_it "\n curr_gene_time_series \n";
			# print_Dumper ( $curr_gene_time_series );
						
			print_it "CRA curr_cluster $curr_cluster " .
					"$ref_array_KP_gene_identifers->[$curr_gene] " .
					"$ref_hash_KP_to_FB->{$ref_array_KP_gene_identifers->[$curr_gene]} ".
					"curr_gene $curr_gene\n";
	
#  here is cmd to view output of above print into file				
#  cat /home/scott/Documents/data/alonso_lab/repeat_v3_fast_35_50_144_run06/standard_out.txt | grep '^CRA'|head
#  which is prosecuted by 
#  perl statistics_bootstrapping.pl -input_dir_prefix /home/scott/Documents/data/alonso_lab/repeat_v3_fast_35_50_144_run| grep '^SBS'  > /home/scott/Documents/data/alonso_lab/repeat_v3_fast_35_50_144_run01/synth_op.txt

					
															
			my $index_timepoint = 0;
			foreach my $curr_Y_value ( @{$curr_gene_time_series}) {
				
				#  iterate across all time points this gene
				
#				print_it "index_timepoint $index_timepoint    curr_Y_value   $curr_Y_value and no other\n";
				
				my $ref_array_gene_values_this_timepoint;
				
				if (exists $hash_timepoints_2_gene_values_this_cluster{$index_timepoint}) {
					
					$ref_array_gene_values_this_timepoint = $hash_timepoints_2_gene_values_this_cluster{$index_timepoint};
				}
				push @{$ref_array_gene_values_this_timepoint}, $curr_Y_value;
				
				$hash_timepoints_2_gene_values_this_cluster{$index_timepoint} = $ref_array_gene_values_this_timepoint;
				
				$index_timepoint++;
			}
			$count_num_genes_this_cluster++;
						
		}  #  $curr_gene
		
#		print_it "curr_cluster $curr_cluster hash_timepoints_2_gene_values_this_cluster\n";
#		print_Dumper ("hash_timepoints_2_gene_values_this_cluster", \%hash_timepoints_2_gene_values_this_cluster );
		
		my $ref_array_medians;
		
		my $total_euclidean_distance = 0;
		my $total_time_points_per_cluster = 0;
				
		foreach my $curr_timepoint (sort { $a <=> $b } keys %hash_timepoints_2_gene_values_this_cluster) {
			
			my $curr_median = Algorithm::Cluster::median($hash_timepoints_2_gene_values_this_cluster{$curr_timepoint});
			
			$ref_array_medians->[$curr_timepoint] = $curr_median;
			
#			print_it "curr_timepoint $curr_timepoint  median_value $ref_array_medians->[$curr_timepoint]\n";
			
#			print_it "hash_timepoints_2_gene_values_this_cluster\n";
#			print_Dumper ( $hash_timepoints_2_gene_values_this_cluster{$curr_timepoint} );
			
			my $abs_diff_median_2_this_gene_curr_timepoint = 0;
			
			foreach my $curr_Y_value ( @{$hash_timepoints_2_gene_values_this_cluster{$curr_timepoint}} ) {
								
				$abs_diff_median_2_this_gene_curr_timepoint += abs($curr_median - $curr_Y_value);

#				print_it "median $curr_median Y $curr_Y_value abs_dif $abs_diff_median_2_this_gene_curr_timepoint\n";
			}
			my $avg_diff_this_timepoint = $abs_diff_median_2_this_gene_curr_timepoint / $count_num_genes_this_cluster;
			
#			print_it "avg_diff_this_timepoint $avg_diff_this_timepoint\n";
			
			$total_euclidean_distance += $avg_diff_this_timepoint;	
			$total_time_points_per_cluster++;
		}	
#		print_it "ref_array_medians\n";
#		print_Dumper ( $ref_array_medians );


#		print "PRE curr_cluster $curr_cluster  total_euclidean_distance $total_euclidean_distance  " .
#			"total_time_points_per_cluster $total_time_points_per_cluster\n";
			
		
		#  divide by number of timepoints so fitness is independent of num of timepoints
		$total_euclidean_distance /= $total_time_points_per_cluster;
		
#		print "POST curr_cluster $curr_cluster  total_euclidean_distance $total_euclidean_distance  " .
#			"total_time_points_per_cluster $total_time_points_per_cluster\n";

		$ref_hash_cluster_fitness_2_cluster = 
			add_hash_of_array_elements($ref_hash_cluster_fitness_2_cluster,
										$total_euclidean_distance, $curr_cluster);  #  stens 20100319
										
#		print_Dumper("ref_hash_cluster_fitness_2_cluster", $ref_hash_cluster_fitness_2_cluster);
		
	}  #  $curr_cluster
	
#	print_Dumper ("ref_hash_cluster_fitness_2_cluster", $ref_hash_cluster_fitness_2_cluster );
	
#	show_numeric_ordered_ref_hash_key_values_array( "ref_hash_cluster_fitness_2_cluster", 
#													$ref_hash_cluster_fitness_2_cluster);
													
	my $num_clusters_to_keep = int($num_clusers * $percentage_cutoff_best_clusters + 1);
	
	print "num_clusters_to_keep $num_clusters_to_keep\n";
													
	#  $ref_hash_mapping_cluster_to_genes
	
	
#	RRRRRRR   need to insert current flavor of which outer + inner loop label into filename of curve fitness
#             so we can keep separate each of these iterations for downstream processing
		
	my $fh_genes_in_best_clusters_fitness_best_curve_output = 
		IO::File->new( ">> $cluster_op_file")
			or die "ERROR - cannot open file cluster_op_file $cluster_op_file $!\n";

	my $curr_num_cluster = 0;
	foreach my $curr_key_fitness (sort { $a <=> $b } keys %{$ref_hash_cluster_fitness_2_cluster} ) {
		
		if ( $curr_num_cluster < $num_clusters_to_keep ) {
			
			print "top cluster fitness " . $curr_key_fitness . " from ";
			
			foreach my $curr_cluster (@{$ref_hash_cluster_fitness_2_cluster->{$curr_key_fitness}} ) {
				
				print "curr_cluster $curr_cluster \n";
				
				foreach my $curr_gene (@{$ref_hash_mapping_cluster_to_genes->{$curr_cluster}}) {
					
#					print "fit ${curr_key_fitness} " .
#						"cluster ${curr_cluster}  gene ${curr_gene} " .
#						$ref_array_KP_gene_identifers->[$curr_gene] ."\n";

					# gene  num_clusers    num_time    start_time    stop_time    fit   cluster

					print $fh_genes_in_best_clusters_fitness_best_curve_output 
						"$ref_array_KP_gene_identifers->[$curr_gene]\t$num_clusers\t" .
						"$minimum_num_timepoints\t$curr_start_timepoint\t$curr_stop_timepoint\t" .
						"${curr_key_fitness}\t${curr_cluster}\n";
				}
			}
#			print "\n";
		}
		
		$curr_num_cluster++;
	}
	return $ref_hash_cluster_fitness_2_cluster;
	
}  #  end of identify_cluster_fitness_based_on_curve_shape



sub do_Algorithm_Cluster {
	
	my ( $num_clusers, $ref_array_time_series, $given_mask, 
		 $given_weight, $ref_array_KP_gene_identifers, $number_of_clustering_passes, $debug ) = @_;

	$debug = $NO;  #  if $NO will turn off debug print
	
	# my %hash_KP_
	my %params = (
				   nclusters => $num_clusers,
				   transpose => 0,
				   npass     => $number_of_clustering_passes,
				   method    => 'a',
				   dist      => 'e',
	);
	my %hash_mapping_cluster_to_genes;
	my ( $ref_array_mapping_gene_to_cluster, $gen_error, $gen_found ) =
	  Algorithm::Cluster::kcluster(
				%params,
				data   => $ref_array_time_series,
				mask   => $given_mask,
				weight => $given_weight,
	);
	
#	print_it "ref_array_mapping_gene_to_cluster Dumper\n";

#	print_Dumper("ref_array_mapping_gene_to_cluster", $ref_array_mapping_gene_to_cluster);
	
#	print_it "\n" ;
#	print_it "Clustering - which Gene in which cluster :\n\n";

	my $curr_gene = 0;
	foreach my $curr_cluster ( @{$ref_array_mapping_gene_to_cluster} ) {
		
		my @array_genes_in_this_cluster;
		if (exists $hash_mapping_cluster_to_genes{$curr_cluster}) {
			
			@array_genes_in_this_cluster = @{$hash_mapping_cluster_to_genes{$curr_cluster}};
		}
		push @array_genes_in_this_cluster, $curr_gene;
		$hash_mapping_cluster_to_genes{$curr_cluster} = \@array_genes_in_this_cluster;
		print_it "gene $curr_gene " . $ref_array_KP_gene_identifers->[$curr_gene] .
				" belongs to cluster $curr_cluster\n" if ( $debug eq $YES);
		$curr_gene++;
	}
	
	print_it "Within-cluster sum of distances is $gen_error\n";
	
	
	return (	$ref_array_mapping_gene_to_cluster, $gen_error, 
				$gen_found, \%hash_mapping_cluster_to_genes );
	
}  #  end of do_Algorithm_Cluster


sub add_hash_of_array_elements {
	
	my ($given_ref_hash_of_array_elements, $curr_key, $curr_value) = @_;
	
	my $ref_array_curr_elements;
	if (exists $given_ref_hash_of_array_elements->{$curr_key}) {
		
		$ref_array_curr_elements = $given_ref_hash_of_array_elements->{$curr_key};
	}
	push @{$ref_array_curr_elements}, $curr_value;
	$given_ref_hash_of_array_elements->{$curr_key} = $ref_array_curr_elements;
	
	return $given_ref_hash_of_array_elements;
}

sub show_numeric_ordered_ref_hash_key_values_array {
	
	my ( $given_label, $given_ref_hash_key_values_array ) = @_;
	
	foreach my $curr_key (sort { $a <=> $b } keys %{$given_ref_hash_key_values_array} ) {
		
		print_it "$given_label $curr_key @{$given_ref_hash_key_values_array->{$curr_key}}\n";
	}
	
}


1;  #  loaded OK