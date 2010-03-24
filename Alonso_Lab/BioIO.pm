
package BioIO;

use warnings;
use Data::Dumper;
use strict;
use Exporter qw(import);


# use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(reorganize_data_structure flip_hash_into_hash_values_of_ref_arrays_of_keys);
our %EXPORT_TAGS = ( DEFAULT => [qw(&reorganize_data_structure)],
                 Both    => [qw(&reorganize_data_structure &flip_hash_into_hash_values_of_ref_arrays_of_keys)]);

my $YES = 'YES';
my $NO  = 'NO';

my $fh_standard_out;  #  standard out is send into this supplied file handle
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
	
sub print_it {
	
	my ( $stuff_to_print ) = @_;
	
	if ( $send_print_to_file eq $YES ) {
		
		print $stuff_to_print;  #  also send to standard out as well as to file
		print $fh_standard_out $stuff_to_print;
		
	} else {

		print $stuff_to_print;
	}
}

sub parse_NCBI_timeseries_dataset {
	
	my (	$self, $ncbi_dataset, $init_num_timepoints, $ref_hash_KP_to_FB, 
			$ref_hash_gene_KB_to_3_prime_UTR_s, $num_data_lines, $length_KP_gene_id,
			$mask_value, $value_weight_per_this_timepoint) = @_;
	
	my ( $ref_array_time_series, $mask1, $weight1, $num_timepoints, 
			$ref_array_KP_gene_identifers );
			
	my $flag_have_we_already_done_KP_validation = $NO;
	my %hash_validated_KP_genes_with_sequence;
	
	
	#  ref_array_time_series        - Time Series of Y values per gene in same sorted order as 
	#                                 ref_array_KP_gene_identifers
	#  ref_array_KP_gene_identifers - Ordered list of all genes using KP gene identifiers
	
	( $ref_array_time_series, $mask1, $weight1, $num_timepoints, 
		$ref_array_KP_gene_identifers ) =
			&parse_ncbi_geo( $ncbi_dataset, 0, $init_num_timepoints,
						$flag_have_we_already_done_KP_validation, 
						\%hash_validated_KP_genes_with_sequence, $length_KP_gene_id, $mask_value,
						$value_weight_per_this_timepoint );

	if ( $init_num_timepoints != 0 ) {
		
		$num_timepoints--;
	}
	# print_it "post run of parse_ncbi_geo calculated num_timepoints $num_timepoints\n";
	# print_it "to show ref_array_time_series ", scalar(@{$ref_array_time_series}), "\n";
	# print_it "to show ref_array_KP_gene_identifers ", scalar(@{$ref_array_KP_gene_identifers}), "\n";
	
#	print_Dumper ( "ref_array_KP_gene_identifers", $ref_array_KP_gene_identifers );  
	
	# --- confirm KP gene identifiers which have time series points 
	#     are found on conversion as well as have sequence data
	#     All genes which fail this criteria are ignored in all further processing.
	
	# $ref_array_KP_gene_identifers --> $ref_hash_KP_to_FB --> $ref_hash_gene_KB_to_3_prime_UTR_s
		  
	my $index_num_valid_genes = 0;
	foreach my $KP_id ( @{$ref_array_KP_gene_identifers}) {
		
		if (exists $ref_hash_KP_to_FB->{$KP_id}) {
			
			# OK does contain KP_id
			
			my $FlyBase_id = $ref_hash_KP_to_FB->{$KP_id};
			
			# print_it "FlyBase_id $FlyBase_id\n";
			
			if (exists $ref_hash_gene_KB_to_3_prime_UTR_s->{$FlyBase_id}) {
				
				# OK does contain sequence data for FlyBase_id
				
				$hash_validated_KP_genes_with_sequence{$KP_id} = 1; # these KP ids have been validated
				
				my @array_sequence = @{$ref_hash_gene_KB_to_3_prime_UTR_s->{$FlyBase_id}};
				
				# print_it "array_sequence @array_sequence\n";
				
				$index_num_valid_genes++;
				if ( $num_data_lines > 0 && $index_num_valid_genes > $num_data_lines) {
					
					# print_it "found desired number of genes $index_num_valid_genes\n";
					last;
				}
				
			} else {
				
				# print_it "WARNING - did NOT find sequence data for FlyBase_id $FlyBase_id in ref_hash_gene_KB_to_3_prime_UTR_s\n";
			}
			
		} else {
			
			# print_it "WARNING KP_id $KP_id is NOT found in ref_hash_KP_to_FB\n";
		}
	}
	# --- end of confirm
	
	#  number of KP in time series     : 14064
	#  number of KP with sequence data : 10581 # seeing 4,000 genes unable to resolve FIX stens todo

	print_it "\n";
	print_it "index_num_valid_genes $index_num_valid_genes\n";
	print_it "number of KP in time series     : ". scalar(@{$ref_array_KP_gene_identifers}). "\n";
	print_it "number of KP with sequence data : ". scalar(keys %hash_validated_KP_genes_with_sequence). "\n\n";
	
	# now RE-run the time series file read to ignore KP not correctly validated
		
	$flag_have_we_already_done_KP_validation = $YES; # validation done now repopulate data structures
	
	( $ref_array_time_series, $mask1, $weight1, $num_timepoints, 
		$ref_array_KP_gene_identifers ) =
			&parse_ncbi_geo( $ncbi_dataset, $num_data_lines, $init_num_timepoints,
						$flag_have_we_already_done_KP_validation, 
						\%hash_validated_KP_genes_with_sequence, $length_KP_gene_id, 
						$mask_value, $value_weight_per_this_timepoint );
	
	if ( $init_num_timepoints != 0 ) {
		
		$num_timepoints--;
	}
	
#	print_Dumper ( "ref_array_time_series", $ref_array_time_series);
	
	return ( $ref_array_time_series, $mask1, $weight1, $num_timepoints, 
		$ref_array_KP_gene_identifers );
	
}  #  parse_NCBI_timeseries_dataset


sub parse_sequence_file_fasta_format {

	my ($self, $three_prime_sequence_file, $num_sequence_lines, $tag_search_pattern ) = @_;
	
	print_it "\ntag_search_pattern $tag_search_pattern inside parse_sequence_file_fasta_format\n\n";
	
	my %hash_gene_id_2_array_of_sequences; # hash to store 3' UTR sequence for each gene identifer

	# http://www.bioperl.org/wiki/HOWTO:SeqIO
	my $seqio_obj = Bio::SeqIO->new(-file => $three_prime_sequence_file, -format => "fasta" );
	my $seq_obj;
	my $curr_table_line_count = 0;
	my $did_we_find_parent; # this assertion will alert when we do NOT find parent gene identifer
	while ($seq_obj = $seqio_obj->next_seq) {
		
		$curr_table_line_count++;
		$did_we_find_parent = $NO;

	    # print_it the sequence   
	    # print_it "new sequence\n";

		#  print_Dumper ( $seq_obj ); # priceless to see definition of this !!!
			    
		# print_it "display_id ", $seq_obj->display_id,"\n";
		# print_it "primary_id ", $seq_obj->primary_id,"\n";
		
		# type=three_prime_untranslated_region; loc=X:complement(15246390..15247138); name=acj6-RE; MD5=d5df61ce90957c5d23b784c66d3d57a0; length=749; parent=FBgn0000028; release=r5.24; species=Dmel;
		
#		print_it "desc ", $seq_obj->desc,"\n";
#		print_it "seq ", $seq_obj->seq,"\n";
		# print_it "accession_number ", $seq_obj->accession_number,"\n";
		# print_it "alphabet ", $seq_obj->alphabet,"\n";
		
		# 'desc' => 'type=three_prime_untranslated_region; loc=2R:complement(3207059..3208269,3209088..3209154,3210065..3210223,3210281..3210461,3210527..3210837,3211324..3211473,3212087..3212242,3213657..3213869,3216794..3216960,3217018..3217194,3217248..3217356,3217422..3218105,3218164..3218319,3218481..3219707,3219766..3219885,3230088..3230138); name=Dscam-RU; MD5=095978902ea0dc6bee5e611606efc5db; length=5139; parent=FBgn0033159; release=r5.24; species=Dmel; ',
		
		
		# now parse value of desc as it contains parent gene id which links to outside world
		
		my ( @array_desc_elements ) = split(/;/, $seq_obj->desc); # pop array with each element of desc
		
		my $parent_gene_identifer;
		my $stub;
		my @array_sequence;
		foreach my $curr_element ( @array_desc_elements ) {
			
			if ($curr_element =~ /$tag_search_pattern/) {
				
				( $stub, $parent_gene_identifer) = split('=', $curr_element);
				
#				print_it "\nseeing tag_search_pattern $tag_search_pattern $parent_gene_identifer\n";
				
				if ( defined $parent_gene_identifer && length($parent_gene_identifer) > 0) {
					
					$did_we_find_parent = $YES;
					my $does_sequence_already_exist_for_this_gene_id = $NO;
					if (exists $hash_gene_id_2_array_of_sequences{$parent_gene_identifer}) {
						
						@array_sequence = @{$hash_gene_id_2_array_of_sequences{$parent_gene_identifer}};
						
						# --- check whether current seq is already stowed for this gene id
						
						foreach my $curr_sequence (@array_sequence) {
							
							if ($curr_sequence eq  $seq_obj->seq) {
								
								$does_sequence_already_exist_for_this_gene_id = $YES;
								
#								print_it "NOTICE - seeing same sequence previously stored for this gene id " .
#									"$parent_gene_identifer\n" .
#									"curr_sequence $curr_sequence\n";
							}
						}
					}
					if ( $does_sequence_already_exist_for_this_gene_id eq $NO) {
						
						my $curr_sequence = $seq_obj->seq;
						push @array_sequence, $curr_sequence;
						
#						print_it "curr_sequence $curr_sequence\n";
					}
					$hash_gene_id_2_array_of_sequences{$parent_gene_identifer} = \@array_sequence;
				}
			}
		}
		
		# --- following is error checking and assertions --- #
		
		if ( $did_we_find_parent eq $NO ) {

			print_Dumper ("seq_obj", $seq_obj );
			
			die "ERROR - failed to find parent gene identifier on row :\n";
		}
		if ( $num_sequence_lines > 0 ) {    #  do we care about limiting number of rows to visit
		
			if ( $curr_table_line_count > $num_sequence_lines ) {
				
				print_it "stopping while loop since saw limit of num_sequence_lines " .
					"$num_sequence_lines with curr_table_line_count $curr_table_line_count\n";
				last;
			}
		}
	}
	
	######### stens todo cleanup logic in parse_conversion_file to handle this :
	# FlyBase_id CG17023:FBgn0024804,CG40120:FBgn0058120,CG40006:FBgn0058006
	# WARNING - did NOT find sequence data for FlyBase_id CG17023:FBgn0024804,CG40120:FBgn0058120,CG40006:FBgn0058006 in ref_hash_gene_KB_to_3_prime_UTR_s
	
	# print_it "about to show hash_gene_id_2_array_of_sequences\n";
	# print_Dumper ( \%hash_gene_id_2_array_of_sequences );
	
	return ( \%hash_gene_id_2_array_of_sequences );
	
}  #  end of parse_sequence_file_fasta_format

sub parse_conversion_file {
	
	my ($self, $gene_id_conversion_file, $length_FlyBase_id,
			$flag_filter_dupe_mappings_KP_into_FB) = @_;
	
	my %hash_KP_to_FB; # mapping between KP gene identifer and its FlyBase FB id
	my %hash_FlyBase_to_array_of_KP;  #  used to identify multiple KP to FB mappings - only keep 1st mapping - filter others
	
	my $fh_gene_id_conversion_file_input = IO::File->new( "< $gene_id_conversion_file" )
	  or die "ERROR - cannot open file gene_id_conversion_file $gene_id_conversion_file\n";
	my $have_we_seen_entire_table = $NO;

	my $count_dupe_mapping_KP_2_FB = 0;
	my $count_failed_lines = 0;
	my $count_succeeded_lines = 0;
	my $count_num_lines = 0;
	while (<$fh_gene_id_conversion_file_input>) {
		
		#  iterate across entire input file
		
		my $curr_conversion_line = $_;
		
		next if ( m/^\#/);  #  skip over lines starting with symbol # they're comments
		
		chomp;  #  strip off trailing newline char
		
#		print_it "curr_conversion_line $curr_conversion_line\n";
		
		# KP=KP07686-454	FBgn=FBgn0032375	CG=CG14932
		my ( $KP_ID_raw, $FlyBase_ID_raw, $CG_ID_raw ) = split('\t', $curr_conversion_line);

		if (	defined $KP_ID_raw && $KP_ID_raw =~ /KP/ && 
				defined $FlyBase_ID_raw && $FlyBase_ID_raw =~ /FBgn/ &&
				defined $CG_ID_raw && $CG_ID_raw =~ /CG/ 
		) {
			
# print_it "$curr_conversion_line";
# print_it "parsed: KP ->$KP<- "; # seeing KP=KPIsoform039  KP=KPST2
# print_it "FlyBase ->$FlyBase<- "; # seeing FBgn=   FBgn=CG31624-RA:FBgn0051624,CG31988-RA:FBgn0051988

# FBgn=CG3560-RA:FBgn0030733,CG32576-RA:FBgn0052576
#  FBgn=CG4026-RA:FBgn0032147,CG18854-RA:FBgn0042174,CG18854-RC:FBgn0042174,CG18854-RB:FBgn0042174
# FlyBase ->FBgn=CG32823-RB:FBgn0052823,CG18000-RK:FBgn0013761,CG33497-RA:FBgn0053497,CG33499-RA:FBgn0053499,CG18000-RA:FBgn0013761,CG18000-RB:FBgn0013761,CG18000-RC:FBgn0013761,CG18000-RD:FBgn0013761,CG18000-RE:FBgn0013761,CG18000-RF:FBgn0013761,CG18000-RG:FBgn0013761,CG18000-RH:FBgn0013761,CG18000-RI:FBgn0013761,CG18000-RJ:FBgn0013761,CG9580-RA:FBgn0025801,CG18000-RL:FBgn0013761,CG18000-RM:FBgn0013761<-
			# print_it "CG ->$CG<-";    # seeing CG=   CG=CG1034,CG14578,CG10837
			# print_it "\n";
			
			
			#            'KP12754-750' => 'CG10246-RA:FBgn0013771,CG10247-RA:FBgn0033981',
			
			my ( $stub_kp, $KP_id) = split(/=/, $KP_ID_raw);
			my ( $stub_fb, $FlyBase_id) = split(/=/, $FlyBase_ID_raw);
			
#			print_it "KP_id $KP_id   FlyBase_id $FlyBase_id\n";
			
			
			#    11 FB id  
			
			if ( $KP_id =~ /Isoform/ && $FlyBase_id eq '' ) {
				
#				print_it "WARNING KP_id $KP_id contain Isoform so no FB ... skipping $curr_conversion_line";
#				print_it "Isoform     KP_id $KP_id\n";
#				print_it "Isoform     FlyBase_id $FlyBase_id\n\n";
				
				$count_failed_lines++;
				
			} elsif ( $FlyBase_id eq '' ) {
				
#				print_it "WARNING - seeing blank FlyBase_id on KP_id $KP_id ... skipping as no FB $curr_conversion_line";
				$count_failed_lines++;
				
			} elsif ( exists $hash_KP_to_FB{$KP_id}) {
				
				die "ERROR - seeing KP_id $KP_id previously exists in hash_KP_to_FB ... skipping as no KP id $curr_conversion_line";
				
			} else {
				
				if ( length($FlyBase_id) > $length_FlyBase_id ) {
					
					print_it "WARNING - extra long FB with KP_id $KP_id  FlyBase_id $FlyBase_id\n";
					$count_failed_lines++;
										
				} else {
				
					if (exists $hash_FlyBase_to_array_of_KP{$FlyBase_id} &&
						$flag_filter_dupe_mappings_KP_into_FB eq $YES) {
							
						print "NOTICE - filter on dupe mapping KP_id $KP_id into  FlyBase_id $FlyBase_id\n";
						$count_dupe_mapping_KP_2_FB++;
						
					} else {
						
						$hash_KP_to_FB{$KP_id} = $FlyBase_id;
						
						push @{$hash_FlyBase_to_array_of_KP{$FlyBase_id}}, $KP_id;
						
	#					print_it "just saved hash_KP_to_FB with KP_id $KP_id    FlyBase_id $FlyBase_id\n";
						
						$count_succeeded_lines++;
					}
				}
				
			}
			
		} else {
			
#			print_it "nonstandard: $curr_conversion_line";
			$count_failed_lines++;
		}
		$count_num_lines++;
	}
	print_it "end of parse_conversion_file of $count_num_lines lines == " .
		"$count_succeeded_lines and $count_failed_lines and " .
		"count_dupe_mapping_KP_2_FB $count_dupe_mapping_KP_2_FB\n";

	return ( \%hash_KP_to_FB, \%hash_FlyBase_to_array_of_KP );
	
}  #  end of parse_conversion_file


sub parse_ncbi_geo {
	
	my ( $ncbi_dataset, $num_data_lines, $num_timepoints,
			$flag_have_we_already_done_KP_validation, 
			$ref_hash_validated_KP_genes_with_sequence, $length_KP_gene_id, $mask_value,
			$value_weight_per_this_timepoint ) = @_;

	my ( @array_time_series, @array_mask, @array_weight_per_timepoint, @array_KP_gene_identifers );


	my @array_treatments;
	my $table_line_limit = $num_data_lines
	  ;    #  if not 0 then limit number of rows read from input file table
	my $curr_table_line_count = 0;

	my $artificially_limit_num_timepoints =
	  $num_timepoints;  #  if not 0 then limit number of columns from input file
	  
	print_it "artificially_limit_num_timepoints $artificially_limit_num_timepoints\n";
	
	my $curr_artificially_limit_num_timepoints;
	my $fh_ncbi_dataset_input = IO::File->new( "< $ncbi_dataset" )
	  or die "ERROR - cannot open file ncbi_dataset $ncbi_dataset\n";
	my $have_we_seen_entire_table = $NO;
	my $doing_table_mode          = $NO;

	#  my $line_count = 0;
	while (<$fh_ncbi_dataset_input>) {
		
		#  iterate across NCBI input file

		# $line_count++;
		# print_it "curr line is  $_\n";
		
		# $curr_artificially_limit_num_timepoints = 0;
		$_ =~ s/"//g;    #  remove the double quote as in "
		
		chomp;  #  remove trailing newline off end of line
		
		if (/series_matrix_table_begin/)
		{
			print_it "seeing series_matrix_table_begin\n";
			$doing_table_mode = $YES;
			next;
		}
		elsif (/series_matrix_table_end/)
		{
			print_it "seeing series_matrix_table_end $curr_table_line_count $curr_table_line_count\n";
			$doing_table_mode = $NO;
			$have_we_seen_entire_table =
			  $YES;      #  indicates we viewed all rows of input file
		}
		elsif ( $doing_table_mode eq $YES )
		{

# my ( $first_element_of_row, @array_elements ) = split;  #  skip first element by putting into $ignore_first_element
			my ( $first_element_of_row, @array_elements ) = split('\t')
			  ;    #  skip first element by putting into $ignore_first_element
			
			if (/ID_REF/) {

				# print_it "seeing table header row $_\n";
				@array_treatments = @array_elements;
				
			} else {
				
				# -------------- now parse table -------------- #

				my ( $stub, $gene_identifer_KP_format ) = split(/=/, $first_element_of_row);
				
#				print_it "first_element_of_row $first_element_of_row ";
#				print_it "KP_format $gene_identifer_KP_format " . length($gene_identifer_KP_format) . "\n";
				if ( $length_KP_gene_id != length($gene_identifer_KP_format)) {
					
					die "ERROR - invalid length of KP formatted gene identifier - it must be length $length_KP_gene_id\n";
				}
			
				if ($flag_have_we_already_done_KP_validation eq $YES) {
					
					if (exists $ref_hash_validated_KP_genes_with_sequence->{$gene_identifer_KP_format}) {
						
						# OK validation was done - and this KP gene id is OK
						
					} else {
						#  stens todo - look at this  print_it "NOTICE this KP is not valid $gene_identifer_KP_format ... skipping \n";
						next;
					}
				}
			
				$curr_table_line_count++;
	# print_it "$curr_table_line_count line to parse table assume full row here\n";
				if ( $table_line_limit > 0 )
				{    #  do we care about limiting number of rows to visit
					if ( $curr_table_line_count > $table_line_limit )
					{
						print_it "stopping while loop since saw limit of " .
							"table_line_limit $table_line_limit with " .
							"curr_table_line_count $curr_table_line_count\n";
						$have_we_seen_entire_table =
						  $YES;  #  satisfies flag on seeing table since we have
						last;
					}
				}

				# my %hash_one_row;
				my ( @array_one_row_data, @array_one_row_mask );
				my $index_element = 0;

	   # print_it "array_elements is ",  scalar(@array_elements), " and no more\n";
				if ( scalar(@array_elements) > 0 )
				{
					$curr_artificially_limit_num_timepoints = 0;
				}

				# print_it "array_elements @array_elements\n";
				my $previous_time_series_data_point = '';
				foreach my $curr_time_series_data_point (@array_elements)
				{
					if ( $curr_time_series_data_point eq '' )
					{
						print_it "Seeing empty element on curr_table_line_count $curr_table_line_count\n";
						$curr_time_series_data_point = $previous_time_series_data_point;
						if ( $previous_time_series_data_point eq '' )
						{
							die "ERROR bad input data seeing blank previous_element on row $_\n";
						}
					}
					$curr_artificially_limit_num_timepoints++;
					if ( $artificially_limit_num_timepoints > 0 )
					{
						if ( $curr_artificially_limit_num_timepoints >
							 $artificially_limit_num_timepoints )
						{
							last
							  ; #  have reached limit of num columns to pluck from curr row
						}
					}

	   # print_it"$array_treatments[$index_element] curr_element  $curr_element\n";
	   # $hash_one_row{$array_treatments[$index_element]} = $curr_element;
					push @array_one_row_data, $curr_time_series_data_point;
					push @array_one_row_mask, $mask_value;
					$previous_time_series_data_point = $curr_time_series_data_point
					  ; # will be used as fill-in value for missing values in row
					$index_element++;
				}

				# print_it "curr_artificially_limit_num_timepoints $curr_artificially_limit_num_timepoints\n";
				push @array_time_series, \@array_one_row_data;
				push @array_mask,        \@array_one_row_mask;
				push @array_KP_gene_identifers, $gene_identifer_KP_format;

				# --------------
			}
		}
	}
	if ( $have_we_seen_entire_table eq $NO )
	{
		die "ERROR - failed to see all rows of table\n";
	}

# print_Dumper ( \@array_time_series );
print_it "curr_artificially_limit_num_timepoints $curr_artificially_limit_num_timepoints\n";
	my $index = 1;
	my $pop_this_num_timepoint_weights =
	  $curr_artificially_limit_num_timepoints;
	if ( $artificially_limit_num_timepoints == 0 )
	{
		$pop_this_num_timepoint_weights++;
	}

	# print_it "pop_this_num_timepoint_weights $pop_this_num_timepoint_weights\n";
	for ( $index = 1 ; $index < $pop_this_num_timepoint_weights ; $index++ )
	{
		push @array_weight_per_timepoint, $value_weight_per_this_timepoint;
	}

	# print_it "array_weight_per_timepoint @array_weight_per_timepoint\n";
	return ( \@array_time_series, \@array_mask,
			 \@array_weight_per_timepoint,
			 $curr_artificially_limit_num_timepoints,
			 \@array_KP_gene_identifers );
			 
}  #  end of parse_ncbi_geo


1;  #   loaded OK

