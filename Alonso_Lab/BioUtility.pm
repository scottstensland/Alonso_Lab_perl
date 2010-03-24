
package BioUtility;

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

my $fh_standard_out;     #  standard out is send into this supplied file handle
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
	
	
	

sub reorganize_data_structure {
	
	#  morph input hash into output hash returning identical data restructured
	
	my ($ref_hash_pair_count_2_parent_2_children, $ref_hash_pair_count_2_count_num_children_of_parent) = @_;
	
	foreach my $pairing_count (keys %{$ref_hash_pair_count_2_parent_2_children}) {
		
		my $ref_hash_parents_2_array_of_children_at_this_count = $ref_hash_pair_count_2_parent_2_children->{$pairing_count};
		
		foreach my $FlyBase_id_parent (keys %{$ref_hash_parents_2_array_of_children_at_this_count}) {
			
			my $ref_array_children_genes = $ref_hash_parents_2_array_of_children_at_this_count->{$FlyBase_id_parent};
			
			my $num_children_this_parent_this_pairing_count = scalar(@{$ref_array_children_genes});
			
			# --- done retrieving original datastructure --- now need to populate newly synthesized datastructure
			
#			print "$pairing_count  $FlyBase_id_parent $num_children_this_parent_this_pairing_count @{$ref_array_children_genes}\n";

			my $ref_hash_count_of_children;
			if (exists $ref_hash_pair_count_2_count_num_children_of_parent->{$pairing_count}) {
				
				$ref_hash_count_of_children = $ref_hash_pair_count_2_count_num_children_of_parent->{$pairing_count};
			}
			my $ref_hash_parent;
			if (exists $ref_hash_count_of_children->{$num_children_this_parent_this_pairing_count}) {
				
				$ref_hash_parent = $ref_hash_count_of_children->{$num_children_this_parent_this_pairing_count};
			}
			$ref_hash_parent->{$FlyBase_id_parent} = $ref_array_children_genes;
			$ref_hash_count_of_children->{$num_children_this_parent_this_pairing_count} = $ref_hash_parent;
			$ref_hash_pair_count_2_count_num_children_of_parent->{$pairing_count} = $ref_hash_count_of_children;
		}
	}
	
#	print_Dumper ( "ref_hash_pair_count_2_count_num_children_of_parent", $ref_hash_pair_count_2_count_num_children_of_parent);

#  this routine restructures input datastructure which looks like :	
#
#          '10' => {
#                    'FBgn0040368' => [
#                                       'FBgn0029588',
#                                       'FBgn0029589',
#                                       'FBgn0040038',
#                                       'FBgn0040892'
#                                     ],
#
#    into same data but more useful structure of
#
#          '10' => {
#                    '4' => {
#                             'FBgn0040368' => [
#                                                'FBgn0029588',
#                                                'FBgn0029589',
#                                                'FBgn0040038',
#                                                'FBgn0040892'
#                                              ],
#
#    so we can more easily retrieve top scoring gene pairs
#    For example:  in above score of 10 means parent FBgn0040368 was found to pair with each of its children
#                  across 10 runs of script.  4 indicates number of children at this score.
#                  This permits a sorted ordering of highest scores 10, followed by highest number of children 4

}  #  reorganize_data_structure

sub flip_hash_into_hash_values_of_ref_arrays_of_keys {
	
	my ( $self, $ref_given_hash, $ref_synth_hash_key_2_array_of_values) = @_;
	
	#  transform same data from input into output hash
	#  original set of (keys/values) --> becomes (values/ref_arrays_of_keys)
		
	while (my ($key, $value) = each %{$ref_given_hash}) {
	
#		print "$key $value\n";
	
		my $ref_array_keys_this_value;
		if (exists $ref_synth_hash_key_2_array_of_values->{$value}) {
		
			$ref_array_keys_this_value = $ref_synth_hash_key_2_array_of_values->{$value};
		}
		push @{$ref_array_keys_this_value}, $key;
		$ref_synth_hash_key_2_array_of_values->{$value} = $ref_array_keys_this_value;
	}
	
#	print_Dumper("ref_synth_hash_key_2_array_of_values", $ref_synth_hash_key_2_array_of_values);
}
	
#	sub firstName {
#	    my ( $self, $firstName ) = @_;
#	    $self->{_firstName} = $firstName if defined($firstName);
#	    return $self->{_firstName};
#	}
#
#
#	sub lastName {
#	    my ( $self, $lastName ) = @_;
#	    $self->{_lastName} = $lastName if defined($lastName);
#	    return $self->{_lastName};
#	}
#
#	sub print {
#	    my ($self) = @_;
#	
#	    #print Person info
#	    printf( "Name:%s %s\n\n", $self->firstName, $self->lastName );
#	}

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