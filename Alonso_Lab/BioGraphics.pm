
package BioGraphics;

use strict;
use warnings;
use Data::Dumper;
use Exporter  qw(import);

use POSIX qw/ceil/;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(  do_plot   plot_histogram);
our %EXPORT_TAGS = ( DEFAULT => [qw(&do_plot)],
                  );

my $YES = 'YES';
my $NO  = 'NO';

my $fh_standard_out;  #  standard out is send into this supplied file handle
my $send_print_to_file = $NO;  #  flag to determine whether to send standard out to output file

#  $| = 1;  #  turn buffering off

    use Carp;
    my $Debugging = 0;  #  http://www.perl.com/doc/manual/html/pod/perltoot.html

    sub new {
    	
    	my $proto = shift;
        my $class = ref($proto) || $proto;
    	
    	my $self = {
    		
    		send_print_to_file => undef,
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
	
   
        sub debug {
        my $class = shift;
        if (ref $class)  { confess "Class method called as object method" }
        unless (@_ == 1) { confess "usage: CLASSNAME->debug(level)" }
        $Debugging = shift;
    }
    
        sub DESTROY {
        my $self = shift;
        if ($Debugging) { carp "Destroying $self " . $self->name }
        -- ${ $self->{"_CENSUS"} };
    }
    
sub do_plot {
	
	my ($to_plot_or_not, $num_clusers, $ref_array_mapping_gene_to_cluster,
		$ref_array_time_series, $num_timepoints, $ref_array_KP_gene_identifers,
		$run_name_dir ) = @_;


	if ( $to_plot_or_not eq $YES && $^O eq "MSWin32" ) {
		
		&do_plot_PC( $num_clusers, $ref_array_mapping_gene_to_cluster, $ref_array_time_series, 
					$num_timepoints, $ref_array_KP_gene_identifers );
		
	} elsif ( $to_plot_or_not eq $YES && $^O eq "linux" ) {
		
		# print_it "-------- start of dumping --------- \n";
		
		# print_it "$num_clusers $num_clusers";
		
		# print_it "ref_array_mapping_gene_to_cluster\n";
		# print_Dumper ( $ref_array_mapping_gene_to_cluster );
		
		# print_it "ref_array_time_series\n";
		# print_Dumper ( $ref_array_time_series );
	
		# print_it "--------   end of dumping --------- \n";
		
		
		my $pyxplot_ppl_file_prefix = <<HERE_DOC_PREFIX
set multiplot                                                          
set nodisplay                                                          
#  width=17                                                  
width=8                                     
gold_ratio = 1/((1+sqrt(5))/2)         

         
#  set terminal png
#  set dpi 300                       


set width width
set nokey
HERE_DOC_PREFIX
;

		# print_it "pyxplot_ppl_file_prefix $pyxplot_ppl_file_prefix\n";
		
		my $pyxplot_ppl_file_contents = $pyxplot_ppl_file_prefix;
		
		$pyxplot_ppl_file_contents .= "\n" .
				&do_plot_linux( $num_clusers, $ref_array_mapping_gene_to_cluster, 
					$ref_array_time_series, $num_timepoints, 
					$ref_array_KP_gene_identifers, $run_name_dir ) . "\n";
					
					
		my $pyxplot_ppl_file_suffix = <<HERE_DOC_SUFFIX

# Now that we are finished preparing multiplot,
# turn display on
set display
refresh
HERE_DOC_SUFFIX
;
		$pyxplot_ppl_file_contents .= "\n" . $pyxplot_ppl_file_suffix;
		
		my $output_filename_pyxplot_contents = "${run_name_dir}/body_pyxplot.ppl";

		my $fh_output_filename_pyxplot_contents = IO::File->new( "> $output_filename_pyxplot_contents" )
			or die "ERROR - cannot open file output_filename_pyxplot_contents " .
					"$output_filename_pyxplot_contents\n";
	  
		print $fh_output_filename_pyxplot_contents "$pyxplot_ppl_file_contents";
		print  "\npyxplot formatted file created as :\n$output_filename_pyxplot_contents\n\n";
					
	} elsif ( $to_plot_or_not eq $YES ) {
		
		die "ERROR - your OS is not recognized so cannot do plot\n";
	}
	
}  #  do_plot

sub plot_histogram {
	
	my ($self,
	    $genes_in_best_clusters_fitness_best_curve_output, $run_name_dir,
		$histogram_filename_output, $histogram_filename_pyxplot_output,
		$cluster_parm_label) = @_;	
	
	#    need to muster up mechanism to discover all cluster output files
	#        all of which will have filename prefix $cluster_parm_label
	
	my %hash_KP_id_2_count;

	my $fh_genes_in_best_clusters_fitness_best_curve_output = 
		IO::File->new( "< $genes_in_best_clusters_fitness_best_curve_output" )
		or die "ERROR - cannot open file genes_in_best_clusters_fitness_best_curve_output " .
				"$genes_in_best_clusters_fitness_best_curve_output\n";
				
	my $fh_histogram_filename_output = IO::File->new( "> $run_name_dir/$histogram_filename_output" )
		or die "ERROR - cannot open file histogram_filename_output " .
				"$run_name_dir/$histogram_filename_output\n";				
				
	my $fh_histogram_filename_pyxplot_output = IO::File->new( "> $histogram_filename_pyxplot_output" )
		or die "ERROR - cannot open file histogram_filename_pyxplot_output " .
				"$histogram_filename_pyxplot_output\n";						
				
				
	while(<$fh_genes_in_best_clusters_fitness_best_curve_output>) {
		
		chomp;  #  strip off trailing newline char
		
		# if (/series_matrix_table_begin/)

		my ( $KP_id, $minimum_num_timepoints, $curr_start_timepoint, $curr_stop_timepoint,
		 ) = split('\t');
		
#		print_it "KP_id $KP_id   minimum_num_timepoints $minimum_num_timepoints   " .
#				"curr_start_timepoint $curr_start_timepoint   curr_stop_timepoint $curr_stop_timepoint\n";
				
		$hash_KP_id_2_count{$KP_id}++;
	}
	
#	print_Dumper ( "hash_KP_id_2_count", \%hash_KP_id_2_count);
	
	my %hash_count_to_KP_ids_with_this_count;  #  reverse of %hash_KP_id_2_count in array to handle dupes
					
	while (my ($KP_id, $count_per_KP_id) = each(%hash_KP_id_2_count)) {
		
		push @{$hash_count_to_KP_ids_with_this_count{$count_per_KP_id}}, $KP_id;
	}
			
	foreach my $count_per_KP_id (sort { $a <=> $b } keys %hash_count_to_KP_ids_with_this_count ) {
		
		my @array_KP_ids_at_this_count = @{$hash_count_to_KP_ids_with_this_count{$count_per_KP_id}};
		
		print "count_per_KP_id $count_per_KP_id  " . scalar(@array_KP_ids_at_this_count) ."\n";

		my $count_KP_ids_at_this_count = scalar(@array_KP_ids_at_this_count);
		
		my $curr_result=sprintf("%.1f %.1f\n",$count_per_KP_id, $count_KP_ids_at_this_count);
		
		print $fh_histogram_filename_output $curr_result;
	}
				
	my $histogram_pyxplot_ppl_contents = <<HERE_DOC_HISTOGRAM
set multiplot
set nodisplay
# width=5.4
# gold_ratio = 1/((1+sqrt(5))/2)

# set width width
set xrange [0.1:30.4]
# set yrange [0:55.1]
set yrange [220.0:940.0]
set nokey

# Plot 0 (bottom left)
# set origin 0*width, 0*width*gold_ratio
set xlabel 'x'
set ylabel 'y'
set label 1 '(c)' 8.2,0.9
# plot 'barchart2.dat' with wboxes
# plot 'barchart2_var5.dat' with boxes
plot '$histogram_filename_output' with boxes

# Now that we are finished preparing multiplot,
# turn display on
set display
refresh
HERE_DOC_HISTOGRAM
;

	
	
#		print "pyxplot_ppl_file_prefix $pyxplot_ppl_file_prefix\n";

	print $fh_histogram_filename_pyxplot_output $histogram_pyxplot_ppl_contents;
	
}  #  plot_histogram

sub do_plot_linux {

	my ( $num_clusers, $ref_array_mapping_gene_to_cluster, $ref_array_time_series, 
	     $num_timepoints, $ref_array_KP_gene_identifers, $run_name_dir ) = @_;

	my %hash_cluster_to_array_gene;
	my %hash_gene_to_cluster;
	my $body_pyxplot_contents = "";
	
	my $curr_gene = 0;
	foreach my $curr_num_cluster ( @{$ref_array_mapping_gene_to_cluster} ) {

		# print_it "curr_num_cluster $curr_num_cluster\n";
		$hash_gene_to_cluster{$curr_gene} = $curr_num_cluster;
		my @array_genes;
		if ( exists $hash_cluster_to_array_gene{$curr_num_cluster} )
		{
			@array_genes = @{ $hash_cluster_to_array_gene{$curr_num_cluster} };
		}
		push @array_genes, $curr_gene++;
		$hash_cluster_to_array_gene{$curr_num_cluster} = \@array_genes;
	}
	
	# print_it "hash_cluster_to_array_gene\n";
	# print_Dumper ( \%hash_cluster_to_array_gene );

	my $num_per_side = ceil( sqrt($num_clusers) );

	my $num_x_panel;
	my $num_y_panel;
	my $curr_cluster;
	my $curr_output_filename = "";
	$curr_gene = 0;

	foreach my $curr_cluster (sort {$a<=>$b} keys %hash_cluster_to_array_gene) {

		# print_it "\n ----------- curr_cluster $curr_cluster  -----------\n";

		my @array_all_genes_this_cluster = @{$hash_cluster_to_array_gene{$curr_cluster}};

		# print_it "array_all_genes_this_cluster ", @array_all_genes_this_cluster, "\n";
		
		( $num_x_panel, $num_y_panel ) =
				calc_panel_XY_positions( $num_clusers, $curr_cluster, $num_per_side );
		  
		$num_x_panel--;  #  decrement by 1 to become base 0
		$num_y_panel--;  #  decrement by 1 to become base 0
			
		$body_pyxplot_contents .= "\n\n\# plot $curr_cluster\n" .		
			"set title $curr_cluster\n" .
			"set origin $num_x_panel*width, $num_y_panel*width*gold_ratio\n";
		
		$body_pyxplot_contents .= "plot ";
		

		my $count_num_genes_this_cluster = 0;
		foreach my $curr_gene ( @array_all_genes_this_cluster ) {

			# print_it "OOOOOO curr_cluster $curr_cluster curr_gene $curr_gene\n";

			my $curr_gene_time_series = $ref_array_time_series->[$curr_gene];

			# print_Dumper ( $curr_gene_time_series );
			
			# ~~~~~ #
					
			# print_it "\n curr_gene_time_series \n";
			# print_Dumper ( $curr_gene_time_series );
			
			#  $curr_cluster = $hash_gene_to_cluster{$curr_gene};

	
			my $synth_seq = gen_seq($num_timepoints);
			
			# print_it "synth_seq\n";
			# print_Dumper ( $synth_seq );
			
			#   $body_pyxplot_contents .= 
			
			# RRRRRRRRRR
			
			$curr_output_filename = gen_output_file_this_plot( $run_name_dir, $curr_cluster, $synth_seq, 
										$curr_gene_time_series, $num_x_panel, $num_y_panel, $curr_gene);
										
			$body_pyxplot_contents .= ", " if ($count_num_genes_this_cluster > 0);  #  insert a comma if more than 1 gene this cluster
			$body_pyxplot_contents .= "\'${curr_output_filename}\' with line";
							 
			# print_it "this subplot num_x_panel ", $num_x_panel, "  num_y_panel ", $num_y_panel, " \n";
			
			$count_num_genes_this_cluster++;
			# ~~~~~ #
			
		}  #  $curr_gene
		
	}  #  $curr_cluster
	
	return $body_pyxplot_contents;

}  #  end of do_plot_linux


sub gen_output_file_this_plot {
	
	my ( $run_name_dir, $curr_cluster, $synth_seq, $curr_gene_time_series,
			$num_x_panel, $num_y_panel, $curr_gene ) = @_;
					
	my $output_file_this_plot = "cluster_${curr_cluster}_gene_${curr_gene}_subplot_${num_x_panel}_${num_y_panel}.dat";
	# print_it "output_file_this_plot $output_file_this_plot\n";
	
	my $full_pathname_file_this_plot = "${run_name_dir}/${output_file_this_plot}";
	my $fh_pyxplot_this_subplot_file = IO::File->new( "> $full_pathname_file_this_plot")
	  or die "ERROR - cannot open file full_pathname_file_this_plot $full_pathname_file_this_plot\n";
	
	my $index = 0;
	foreach my $curr_Y_value ( @{$curr_gene_time_series}) {
		
		# print_it "X $index  Y $curr_Y_value\n";
		
		print $fh_pyxplot_this_subplot_file "$index $curr_Y_value\n";
		$index++;
	}
	return $output_file_this_plot;

}  #  gen_output_file_this_plot


sub do_plot_PC {
	
	my ( $num_clusers, $ref_array_mapping_gene_to_cluster, $ref_array_time_series, 
	     $num_timepoints, $ref_array_KP_gene_identifers ) = @_;

	my %hash_cluster_to_array_gene;
	my %hash_gene_to_cluster;
	# print_it "ref_array_mapping_gene_to_cluster\n";
	# print_Dumper ($ref_array_mapping_gene_to_cluster);
	
	# print_it "ref_array_time_series\n";
	# print_Dumper ($ref_array_time_series);
	
	# print_it "ref_array_KP_gene_identifers\n";
	# print_Dumper ($ref_array_KP_gene_identifers);
	
	my $curr_gene = 0;
	foreach my $curr_num_cluster ( @{$ref_array_mapping_gene_to_cluster} ) {

		# print_it "curr_num_cluster $curr_num_cluster\n";
		$hash_gene_to_cluster{$curr_gene} = $curr_num_cluster;
		my @array_genes;
		if ( exists $hash_cluster_to_array_gene{$curr_num_cluster} )
		{
			@array_genes = @{ $hash_cluster_to_array_gene{$curr_num_cluster} };
		}
		push @array_genes, $curr_gene++;
		$hash_cluster_to_array_gene{$curr_num_cluster} = \@array_genes;
	}

	# print_it "seeing post %hash_cluster_to_array_gene\n";
	# print_Dumper ( \%hash_cluster_to_array_gene );
	# print_it "ref_array_time_series\n";
	# print_Dumper ( $ref_array_time_series );
	my $num_per_side = ceil( sqrt($num_clusers) );

	# print_it "num_per_side $num_per_side\n";
	
	my $rate_win = "dummy";
#	my $rate_win =
#	  PDL::Graphics::PGPLOT::Window->new(
#										  Device      => '/gw',
#										  Aspect      => 1,
#										  WindowWidth => 5,
#										  NXPanel     => $num_per_side,
#										  NYPanel     => $num_per_side
#	  );
	my $num_x_panel;
	my $num_y_panel;
	my $curr_cluster;
	$curr_gene = 0;
	foreach my $curr_gene_time_series ( @{$ref_array_time_series} ) {
		
		$curr_cluster = $hash_gene_to_cluster{$curr_gene};
		( $num_x_panel, $num_y_panel ) =
		  calc_panel_XY_positions( $num_clusers, $curr_cluster, $num_per_side );

# print_it "curr_gene $curr_gene curr_cluster $curr_cluster @{$curr_gene_time_series} " .
# 		"num_x_panel $num_x_panel num_y_panel $num_y_panel\n";
# print_it "num_per_side $num_per_side\n";

		my $synth_seq = gen_seq($num_timepoints);

		#print_it "num_timepoints $num_timepoints\n";
		#print_it "synth_seq @{$synth_seq} as array\n";
		#print_it "curr_gene $curr_gene\n";
		#print_it "ref_array_KP_gene_identifers ", $ref_array_KP_gene_identifers->[$curr_gene],"\n";
		
		#print_it "here is X of size ", scalar(@{$synth_seq}), "  ", pdl(@{$synth_seq}), " eodX\n";
		#print_it "here is Y ", scalar(@{$curr_gene_time_series}), "  ", pdl( @{$curr_gene_time_series} ), " eodY\n";
		#print_it "num_x_panel  $num_x_panel   num_y_panel $num_y_panel\n";
		
		$rate_win->line( pdl( @{$synth_seq} ),
						 pdl( @{$curr_gene_time_series} ),
						 { Panel => [ $num_x_panel, $num_y_panel ] } );
		$curr_gene++;
	}

}  #  end of do_plot_PC

sub gen_seq {
	
	my ($given_count) = @_;
	my $curr_count = 0;
	my @array_building;
	while ( $curr_count < $given_count )
	{
		push @array_building, $curr_count;
		$curr_count++;
	}

	# print_it "array is @array_building\n";
	return \@array_building;
}

sub calc_panel_XY_positions {

	# use POSIX qw(ceil floor);
	# use POSIX qw(ceil);
	# use POSIX;
	# given total number of clusters and current cluster number,
	# calculate both X and Y panel numbers for plotting current cluster
	# IE. if we have 4 possible clusters (plotted 2 by 2 )
	#     where XY would be {1,1}, {1,2}
	#                       {2,1}  {2,2}
	#     and we are looking at cluster 3 of 4, then XY would {2,1}
	my ( $total_num_clusers, $curr_cluster, $num_per_side ) = @_;
	$curr_cluster++;    # to start from base of 1 instead of from 0

   # print_it "total_num_clusers $total_num_clusers  curr_cluster $curr_cluster " .
   # 		" num_per_side $num_per_side\n";
   # my $which_row_tmp = 0.0 + ($curr_cluster / $num_per_side);
   # my $which_row = ceil($which_row_tmp);
	my $which_row = ceil( $curr_cluster / $num_per_side );
	my $which_col = $curr_cluster - ( ( $which_row - 1 ) * $num_per_side );

	# print_it "which_row $which_row  which_col $which_col\n\n";
	return ( $which_row, $which_col );
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

1;  #   loaded OK