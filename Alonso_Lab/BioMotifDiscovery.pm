
package BioMotifDiscovery;

use warnings;
use Data::Dumper;
use strict;
use Exporter  qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = ();
our @EXPORT_OK   = qw(de_novo_motif_discovery   firstName lastName print  );
our %EXPORT_TAGS = ( DEFAULT => [qw(&de_novo_motif_discovery)],
                  );
                  
my $YES = 'YES';
my $NO  = 'NO';

my $fh_standard_out;  #  standard out is send into this supplied file handle
my $send_print_to_file = $NO;  #  flag to determine whether to send standard out to output file

#   $| = 1;  #  turn buffering off

    use Carp;
    my $Debugging = 0;  #  http://www.perl.com/doc/manual/html/pod/perltoot.html

    sub new {
    	
    	my $proto = shift;
        my $class = ref($proto) || $proto;
    	
    	my $self = {
    		
#    		_fh_standard_out => undef,
    		
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
    

sub de_novo_motif_discovery {
	
	my ( $self,
		$gene_pairs_counts_across_multi_runs, 
		$ref_hash_gene_KB_to_3_prime_UTR_s, 
		$ref_hash_UTR_segment_freq_success, $ref_hash_UTR_segment_freq_mismatch) = @_;
	
	print "inside de_novo_motif_discovery\n";
	
#	print_Dumper("ref_hash_gene_KB_to_3_prime_UTR_s", $ref_hash_gene_KB_to_3_prime_UTR_s);

	my $fh_gene_pairs_counts_across_multi_runs = 
		IO::File->new( "< $gene_pairs_counts_across_multi_runs" )
		or die "ERROR - cannot open gene_pairs_counts_across_multi_runs " .
				"$gene_pairs_counts_across_multi_runs: $!\n";
	
	#  count of pairs -> count of children > parent > array children
	my %hash_pair_count_2_count_num_children_of_parent;  #  a post production of below preliminary hash
	
	my %hash_pair_count_2_parent_2_children;  #  easy to retrieve best matching genes in sorted order of pairing
	
	while (<$fh_gene_pairs_counts_across_multi_runs>) {
		
		chomp;
		
		my $one_line = $_;
		
		next if (! $one_line =~ m/^SBS/);
		
#		print "$one_line\n";
		
		my ( $stuba, $FlyBase_id_parent, $FlyBase_id_child, $pairing_count) = split(/ /, $one_line);
		
#		print "$FlyBase_id_parent    $FlyBase_id_child   $pairing_count\n";
		
		my $ref_hash_parents_2_array_of_children_at_this_count;
		if (exists $hash_pair_count_2_parent_2_children{$pairing_count}) {
			
			$ref_hash_parents_2_array_of_children_at_this_count = $hash_pair_count_2_parent_2_children{$pairing_count};
		}
		my $ref_array_children_genes;
		if (exists $ref_hash_parents_2_array_of_children_at_this_count->{$FlyBase_id_parent}){
			
			$ref_array_children_genes = $ref_hash_parents_2_array_of_children_at_this_count->{$FlyBase_id_parent};
		}
		push @{$ref_array_children_genes}, $FlyBase_id_child;
		$ref_hash_parents_2_array_of_children_at_this_count->{$FlyBase_id_parent} 
				= $ref_array_children_genes;
		$hash_pair_count_2_parent_2_children{$pairing_count} 
				= $ref_hash_parents_2_array_of_children_at_this_count;
	}
#	print_Dumper("hash_pair_count_2_parent_2_children", \%hash_pair_count_2_parent_2_children);
	
	BioUtility::reorganize_data_structure(\%hash_pair_count_2_parent_2_children, 
										\%hash_pair_count_2_count_num_children_of_parent);
	
#	print_Dumper("hash_pair_count_2_count_num_children_of_parent", \%hash_pair_count_2_count_num_children_of_parent);
	
	seek_patterns( \%hash_pair_count_2_count_num_children_of_parent, 
				   $ref_hash_gene_KB_to_3_prime_UTR_s,
				   $ref_hash_UTR_segment_freq_success, $ref_hash_UTR_segment_freq_mismatch);
				   
#	print_Dumper("about to leave de_novo_motif_discovery hash_UTR_segment_freq_success", 
#					$ref_hash_UTR_segment_freq_success);
				   
	
}  #  de_novo_motif_discovery

sub seek_patterns {
	
	my ( $ref_hash_pair_count_2_count_num_children_of_parent, 
		 $ref_hash_gene_KB_to_3_prime_UTR_s, 
		 $ref_hash_UTR_segment_freq_success, $ref_hash_UTR_segment_freq_mismatch) = @_;
		 
#	my $min_gene_pairing_count_to_process_outer = 9;  #  filter todo dynamically calc this
	my $min_gene_pairing_count_to_process_outer = 3;  #  filter todo dynamically calc this

#	my $min_gene_pairing_count_to_process_inner = 4;  #  filter
	my $min_gene_pairing_count_to_process_inner = 3;  #  filter
	
#	print_Dumper("seek PP patterns  ref_hash_pair_count_2_count_num_children_of_parent", $ref_hash_pair_count_2_count_num_children_of_parent);
#	print_Dumper("seek patterns  ref_hash_gene_KB_to_3_prime_UTR_s", $ref_hash_gene_KB_to_3_prime_UTR_s);

	my $count_num_keys = scalar( keys %{$ref_hash_pair_count_2_count_num_children_of_parent});
	my $curr_num_key = 1;
	foreach my $pairing_count (reverse sort { $a <=> $b } keys %{$ref_hash_pair_count_2_count_num_children_of_parent}) {
		
		return if ( $pairing_count < $min_gene_pairing_count_to_process_outer);  #  filter
		
		print "count_num_keys $count_num_keys   pairing_count $pairing_count\n";
		
#		ssssssssss
				
		my $ref_hash_count_of_children = $ref_hash_pair_count_2_count_num_children_of_parent->{$pairing_count};
					
		foreach my $num_children_this_parent_this_pairing_count (reverse sort { $a <=> $b } keys %{$ref_hash_count_of_children}) {
			
			next if ($num_children_this_parent_this_pairing_count < $min_gene_pairing_count_to_process_inner);
			
			my $ref_hash_parent = $ref_hash_count_of_children->{$num_children_this_parent_this_pairing_count};

			foreach my $FlyBase_id_parent ( sort keys %{$ref_hash_parent}) {
				
				foreach my $FlyBase_id_child ( sort @{$ref_hash_parent->{$FlyBase_id_parent}}) {
					
#					print "seek_patterns_per_prioritized_gene_pair $pairing_count " .
#						"$num_children_this_parent_this_pairing_count $FlyBase_id_parent $FlyBase_id_child\n";
					
					seek_patterns_per_prioritized_gene_pair($FlyBase_id_parent, $FlyBase_id_child,
						$ref_hash_gene_KB_to_3_prime_UTR_s->{$FlyBase_id_parent},
						$ref_hash_gene_KB_to_3_prime_UTR_s->{$FlyBase_id_child},
						$ref_hash_UTR_segment_freq_success,
						$ref_hash_UTR_segment_freq_mismatch, $curr_num_key, $count_num_keys);
				}	
			}	
		}
		$curr_num_key++;
	}
}

sub seek_patterns_per_prioritized_gene_pair {
	
	my ($FlyBase_id_parent, $FlyBase_id_child,
		$ref_hash_gene_KB_to_3_prime_UTR_s_FlyBase_id_parent,
		$ref_hash_gene_KB_to_3_prime_UTR_s_FlyBase_id_child,
		$ref_hash_UTR_segment_freq_success,
		$ref_hash_UTR_segment_freq_mismatch, $curr_num_key, $count_num_keys) = @_;
		
	my $previous_parent_UTR_seq = '';
	foreach my $curr_parent_3_prime_UTR_seq (sort @{$ref_hash_gene_KB_to_3_prime_UTR_s_FlyBase_id_parent}) {
		
		next if($curr_parent_3_prime_UTR_seq eq  $previous_parent_UTR_seq);  #  skip dupe UTR
		
		my $previous_child_UTR_seq = '';
		foreach my $curr_child_3_prime_UTR_seq (sort @{$ref_hash_gene_KB_to_3_prime_UTR_s_FlyBase_id_child}) {
			
			next if ($curr_child_3_prime_UTR_seq eq $previous_child_UTR_seq);  #  skip dupe UTR
			
			print "$curr_num_key of $count_num_keys $FlyBase_id_parent $curr_parent_3_prime_UTR_seq child $FlyBase_id_child $curr_child_3_prime_UTR_seq\n";
			
			seek_patterns_this_pair_of_UTR($FlyBase_id_parent, $curr_parent_3_prime_UTR_seq,
										   $FlyBase_id_child,  $curr_child_3_prime_UTR_seq,
										   $ref_hash_UTR_segment_freq_success,
										   $ref_hash_UTR_segment_freq_mismatch);
										   
			$previous_child_UTR_seq = $curr_child_3_prime_UTR_seq;
		}
		$previous_parent_UTR_seq = $curr_parent_3_prime_UTR_seq;
	}
}

sub seek_patterns_this_pair_of_UTR {
	
	my ($FlyBase_id_parent, $curr_parent_3_prime_UTR_seq,
		$FlyBase_id_child,  $curr_child_3_prime_UTR_seq,
		$ref_hash_UTR_segment_freq_success,
		$ref_hash_UTR_segment_freq_mismatch) = @_;
										   
	my $min_length_of_motif = 8;
	my $max_length_of_motif = 8;
	
	my $parent_UTR_length = length($curr_parent_3_prime_UTR_seq);
	my $child_UTR_length  = length($curr_child_3_prime_UTR_seq);
	
#	print "curr_parent_3_prime_UTR_seq $curr_parent_3_prime_UTR_seq\n";
		
	my $curr_motif_length = $max_length_of_motif;
	while ($curr_motif_length >= $min_length_of_motif) {
		
		my $curr_parent_ptr_left = 0;
		my $curr_parent_ptr_right = $curr_parent_ptr_left + $curr_motif_length;
		
		while ($curr_parent_ptr_right <= $parent_UTR_length) {
			
#			print "$curr_motif_length $curr_parent_ptr_left $curr_parent_ptr_right ";
			
			my $curr_parent_UTR_segment = substr $curr_parent_3_prime_UTR_seq, $curr_parent_ptr_left, $curr_motif_length;
			
#			print "$curr_parent_UTR_segment  <-->  ";
			
			my $curr_child_ptr_left = 0;
			my $curr_child_ptr_right = $curr_child_ptr_left + $curr_motif_length;
			
			while ($curr_child_ptr_right <= $child_UTR_length) {
				
#				print "$curr_motif_length $curr_parent_ptr_left $curr_parent_ptr_right ";
#				print "$curr_parent_UTR_segment  <-->  ";

				
#				print "$curr_motif_length $curr_child_ptr_left $curr_child_ptr_right ";
				
				my $curr_child_UTR_segment = substr $curr_child_3_prime_UTR_seq, $curr_child_ptr_left, $curr_motif_length;
				
#				print "$curr_child_UTR_segment\n";
				
				if ( $curr_parent_UTR_segment eq $curr_child_UTR_segment) {
					
#					print "SUCCESS found match $curr_parent_UTR_segment ===> ";
#					
#					print "$curr_motif_length $curr_parent_ptr_left $curr_parent_ptr_right ";
#					print "$curr_parent_UTR_segment  <-->  ";
#					print "$curr_motif_length $curr_child_ptr_left $curr_child_ptr_right ";
#					print "$curr_child_UTR_segment\n";

					$ref_hash_UTR_segment_freq_success->{$curr_parent_UTR_segment}++;
					$ref_hash_UTR_segment_freq_success->{'total_number_seq_looked_at'}++;
					
				} else {
					
#					print "seeing mismatch\n";
					
					$ref_hash_UTR_segment_freq_mismatch->{$curr_parent_UTR_segment}++;
					$ref_hash_UTR_segment_freq_mismatch->{$curr_child_UTR_segment}++;
					$ref_hash_UTR_segment_freq_mismatch->{'total_number_seq_looked_at'}++;
					$ref_hash_UTR_segment_freq_mismatch->{'total_number_seq_looked_at'}++;
				}
				$curr_child_ptr_left++;
				$curr_child_ptr_right++;
			}
			$curr_parent_ptr_left++;
			$curr_parent_ptr_right++;
		}
		$curr_motif_length--;
	}
}

sub sum_total_motifs_found {
	
	my ( $ref_hash_UTR_motif_count_2_ref_array_motifs_success) = @_;
	
	my $total_motifs = 0;
	foreach my $curr_motif_count ( keys %{$ref_hash_UTR_motif_count_2_ref_array_motifs_success}) {
		
		my $ref_array_curr_motif = $ref_hash_UTR_motif_count_2_ref_array_motifs_success->{$curr_motif_count};
		
#		print "\ncurr_motif_count $curr_motif_count @{$ref_array_curr_motif}\n";
		
		my $count_num_motif_this_count = scalar(@{$ref_array_curr_motif});
			
		$total_motifs += $curr_motif_count * $count_num_motif_this_count;
		
#		print "total_motifs $total_motifs\n";
	}
	return $total_motifs;
}


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