#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;

my $help = 0;
my $man = 0;
my $maker_file = "";
my $web_file = "";
my $output = "";

GetOptions('help|?' => \$help, 
		man => \$man,
		"maker_gff=s" => \$maker_file,
		"web_gff=s" => \$web_file,
		"out_gff=s" => \$output
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(1) unless $maker_file;
pod2usage(1) unless $web_file;
pod2usage(-exitval => 0, -verbose => 2) if $man;

#__END__

=head1 NAME

reannotation_parser.pl


=head1 SYNOPSIS

 perl reannotation_parser.pl [options] --maker_gff <file> --web_gff <file> 
	Options:
	--help		shows help message

needs Bio::FeatureIO

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message

=item B<--maker_gff>

Annotation file (gff3) from maker

=item B<--web_gff>

Reannotation file (gff3) from webapollo

=item B<--out_gff>

Name for outputfile

=back

=head1 DESCRIPTION

The parser is taking two gff3 files as input. One is exported from webapollo and contains reannotated genes. The other gff3 file contains the original maker annotation for the same accession. 
The script tests the webapollo gff3 for correct use of names, substructurs, and tag=value relationships.
If genes do not follow the conventions in the SOP, it reports those genes to STDERR.
The remaining, correct genes are compared to the maker gff3 file, and genes that were reannotated replace the original maker genes. The final gff3 structure is finally printed to the outputfile.

=cut

######
######
my $maker_gff = import_gff($maker_file);
my $web_gff = import_gff($web_file);



###
### Parser for maker gff3 file
###

my $contig="";
my %maker_genes=();
my $maker_gene_curr="";
my %contig_matches=();

print STDERR "read maker file: started\n";
while (my $feature = $maker_gff->next_feature){
	if ($feature->primary_tag eq 'contig'){
		$contig=$feature->{_gsf_seq_id};
		push (@{$contig_matches{$contig}},$feature);
	}
	elsif ($feature->primary_tag eq 'gene'){
		$maker_gene_curr = $feature->{_gsf_tag_hash}{ID}[0];
		push (@{$maker_genes{$maker_gene_curr}},$feature);
	} 
	elsif ($feature->source_tag eq 'maker'){
		push (@{$maker_genes{$maker_gene_curr}},$feature);
	}
	else {
		push (@{$contig_matches{$contig}},$feature);
	}
}
print STDERR "read maker file: done\n";

###
### Parser for webapollo gff3 file: Parses lines into genes and tests for correctness of given names and attributes
###

my %web_genes=();
my $gene_curr="";
my %mRNAs=();
my %genes=();

print STDERR "read webapollo file: started\n";
while (my $feature = $web_gff->next_feature) {
##A few renamings to match maker file structure in the end
	$feature->{_source_tag}="webapollo";
	$feature->{_gsf_tag_hash}{Alias}[0]=$feature->{_gsf_tag_hash}{ID}[0];
	$feature->{_gsf_tag_hash}{ID}[0]=$feature->{_gsf_tag_hash}{Name}[0];
##gene level
	if ($feature->primary_tag eq 'gene'){
		$gene_curr = get_name($feature);
		if (exists($genes{$gene_curr})){
			print STDERR "Duplicated gene '$gene_curr'. Correct and rerun script.\n"} else {
			$genes{$gene_curr}=1;
		}		
		%mRNAs=();
	}
	if ($feature->primary_tag=~/pseudogene|tRNA|snRNA|snoRNA|ncRNA|rRNA|miRNA|repeat_region|transposable_element/){
	my $tag=$feature->primary_tag;	
	my $name=get_name($feature);
	print STDERR "Did not expect to have feature type '$tag' (in gene '$name'). Will result in several follow-up error messages. Rename with 'gene' and rerun script\n";
	}
##mRNA level
	if ($feature->primary_tag eq 'mRNA'){
		my $transcript=get_name($feature); 
		test_names($gene_curr,$transcript,'','');
		if (exists($mRNAs{$transcript})){
			print STDERR "The transcript '$transcript' already exists. Correct and rerun script.\n"} else {
			$mRNAs{$transcript}=1;
		}
	}
##create gff file structure again
	push (@{$web_genes{$gene_curr}},$feature);	
}


print STDERR "read webapollo file: done\n";
##Check correctness of tag-value pairs. If a pair/putpair tag is detected, load also partner gene for testing.
print STDERR "test webapollo file: started\n";
foreach my $key (keys %genes){
	my $paired_gene="";
	if ($web_genes{$key}[0]->{_gsf_tag_hash}{putpair}[0]){
		$paired_gene=$web_genes{$key}[0]->{_gsf_tag_hash}{putpair}[0];
		test_tag_values($web_genes{$key}[0],$web_genes{$paired_gene}[0]);
	}
	elsif ($web_genes{$key}[0]->{_gsf_tag_hash}{pair}[0]){
		$paired_gene=$web_genes{$key}[0]->{_gsf_tag_hash}{pair}[0];
		test_tag_values($web_genes{$key}[0],$web_genes{$paired_gene}[0]);
	}	
	else {
		test_tag_values($web_genes{$key}[0]);
	}
}
print STDERR "test webapollo file: done\n";

###
###Replace maker genes with webapollo genes
###

print STDERR "create final gff file: started\n";
my %final_genes=replace(\%maker_genes,\%web_genes);
###continue here to write output gff3 file
my $gffout = Bio::Tools::GFF->new(-file => ">$output", -gff_version => 3); 
foreach my $contig (keys %contig_matches) {
	$gffout->write_feature($contig_matches{$contig}[0]);
	foreach my $gene (keys %final_genes){
		if ($final_genes{$gene}[0]->{_gsf_seq_id} eq $contig){
			foreach (@{$final_genes{$gene}}) {
				$gffout->write_feature($_)
			}
		}
	}
	shift $contig_matches{$contig};
	foreach (@{$contig_matches{$contig}}){
		$gffout->write_feature($_);
	}
}
print STDERR "create final gff file: done\n";
print STDERR "Script finished\n";









###
###subroutines
###

=head2 IMPORT_GFF

this subroutine imports a gff file and returns its content

=cut

## import a gff file
sub import_gff {
	my $file = $_[0];
	my $in  = Bio::Tools::GFF->new(-file => $file , -gff_version => 3);
#	$in->ignore_sequence(1);
##remove previous line and add next line if you want to use fasta sequences that might be attached at the end of the gff3 file
	$in->features_attached_to_seqs(1);

	return $in;
}




=head2 GET_NAME

Gets the human readable name of a feature

=cut

sub get_name {
	my $feature =$_[0];
	my $name = $feature->{_gsf_tag_hash}{Name}[0];
	return $name;
}




### tests

=head2 TEST_NAMES

If called while reading a mRNA feature: This subroutine reads gene name and corresponding transcript name. It tests if both follow the naming convention from the SOPs and correctly correspond to each other.
If called while reading the value for a 'merged' or 'paired' tag: The subroutine tests if the value follows the naming convention from the SOPs.

=cut

sub test_names {
        my $gene_name = shift @_;
        my $trans_name = shift @_;
        my $value = shift @_;
        my $tag = shift @_;
        if ($tag eq 'merged'){
                if ($gene_name){
                        unless ($value=~/^\d+\|G\d+$/ || $value=~/^\d+\|G\d+\.N\d+$/ || $value=~/^\d+\|G\d+\.\d+$/) {
                                print STDERR "In the %merged% tag for gene '$gene_name', the name '$value' does not follow SOPs. Correct and rerun script\n";}
                } else {print STDERR "Did not find gene name while parsing a %merged% tag. Cannot give more information here. Error will be reported in the transcript feature.\n"}
        } elsif ($tag eq 'paired') {
        	if ($gene_name){
                       unless ($value=~/^\d+\|G\d+$/ || $value=~/^\d+\|G\d+\.N\d+$/ || $value=~/^\d+\|G\d+\.\d+$/) {
                                print STDERR "In the %pair%/%putpair% tag for gene '$gene_name', the name '$value' does not follow SOPs. Correct and rerun script\n";}
                } else {print STDERR "Did not find gene name while parsing a %pair%/%putpair% tag. Cannot give more information. Error will be reported in the transcript feature.\n"}        	
        } else {
                if ($gene_name){
                        unless ($gene_name=~/^\d+\|G\d+$/ || $gene_name=~/^\d+\|G\d+\.N\d+$/ || $gene_name=~/^\d+\|G\d+\.\d+$/) {
                                print STDERR "Gene name '$gene_name' does not follow SOPs. Correct and rerun script\n";}
                } else {print STDERR "Gene name not defined for transcript '$trans_name'. Probably follow up error. Check previous error messages. Correct and rerun.\n"}

                if ($gene_name && $trans_name){
                        unless ($trans_name=~/^\d+\|T\d+-R\d+$/ || $trans_name=~/^\d+\|T\d+\.N\d+-R\d+$/ || $trans_name=~/^\d+\|T\d+\.\d+-R\d+$/) {
                                print STDERR "Transcript name '$trans_name' does not follow SOPs. Correct and rerun script\n";}
                        } else {print STDERR "Gene and/or transcript name not defined. Probably follow up error. Check previous error messages.\n"}

                if ($gene_name){
                        $gene_name=~/(^\d+\|)G(.+)/;
                        unless ($trans_name=~/${1}T${2}/){
                        print STDERR "Transcript '$trans_name' does not correctly correspond to parental gene name '$gene_name'. Correct and rerun script\n";}
                } 
        }
}

=head2 TEST_TAG_VALUES

Ensures that correct tag-value pairs are given

=cut

sub test_tag_values {
	my $feature = shift @_;
	my $paired_gene = shift @_;
	my $fusion = $feature->{_gsf_tag_hash}{fusion}[0];
	my $putpair = $feature->{_gsf_tag_hash}{putpair}[0];
	my $pair = $feature->{_gsf_tag_hash}{pair}[0];
	my $Note = $feature->{_gsf_tag_hash}{Note}[0];
	my $reinspection = $feature->{_gsf_tag_hash}{reinspection}[0];
	my $truncated = $feature->{_gsf_tag_hash}{truncated}[0];
	my $pseudogene = $feature->{_gsf_tag_hash}{pseudogene}[0];
	my $noevidence = $feature->{_gsf_tag_hash}{noevidence}[0];
	my $merged = $feature->{_gsf_tag_hash}{merged}[0];
	my $corbound = $feature->{_gsf_tag_hash}{corbound}[0];
	my $cortrans = $feature->{_gsf_tag_hash}{cortrans}[0];
	my $misassembly = $feature->{_gsf_tag_hash}{misassembly}[0];
	my $delete = $feature->{_gsf_tag_hash}{delete}[0];
	my $date_creation = $feature->{_gsf_tag_hash}{date_creation}[0];
	my $owner = $feature->{_gsf_tag_hash}{owner}[0];
	my $ID = $feature->{_gsf_tag_hash}{ID}[0];
	my $Alias = $feature->{_gsf_tag_hash}{Alias}[0];
	my $Name = $feature->{_gsf_tag_hash}{Name}[0];
	my $date_last_modified = $feature->{_gsf_tag_hash}{date_last_modified}[0];

### test for tags that are not allowed (will mostly be typos)
	my %tags=('fusion'=>1, 'putpair'=>1, 'pair'=>1, 'Note'=>1, 'reinspection'=>1, 'truncated'=>1, 'pseudogene'=>1, 
		'noevidence'=>1, 'merged'=>1, 'corbound'=>1, 'cortrans'=>1, 'misassembly'=>1, 'delete'=>1, 'date_creation'=>1, 
		'owner'=>1, 'ID'=>1, 'Alias'=>1, 'Name'=>1, 'date_last_modified'=>1);
	foreach my $key (keys %{$feature->{_gsf_tag_hash}}) {
		unless (exists $tags{$key}){print STDERR "The tag '$key'  in gene '$Name' is not defined. Allowed tags are:", Dumper(\%tags),"\n.Correct and rerun script.\n"}
		}
### test for correct tag-value relationships. 
## fusions
	if ($Name=~/^\d+\|G\d+\.\d+$/){unless ($fusion){print STDERR "Tag fusion=1 expected for gene '$Name', but not found. Correct and rerun.\n"}}
	if ($fusion){
		unless ($Name=~/^\d+\|G\d+\.\d+$/) {
			print STDERR "The gene '$Name' is tagged as fusion, but its name does not follow the SOP for fusions(^\\d+\\|G\\d+\\.\\d+$). Please check and if necessary correct and rerun script.\n"
		}
	}

## reinspection
	if ($reinspection){print STDERR "Gene '$Name' is tagged with reinspection. The comment is '$Note'. If no comment is given, you will get an 'uninitialized value $Note' error. You might want to finalize the gene model and rerun.\n"}

## pseudogene 
	if ($pseudogene){unless ($pseudogene=~/AT\dG\d+/){print STDERR "Gene '$Name' is tagged as 'pseudogene'. Value '$pseudogene' does not follow the expected regex (AT\\dG\\d+) for an Arabidopsis locus name (e.g. AT5G46470). Correct and rerun script.\n"}}

##merged genes
	if ($merged){
		my @merged_genes=split(' ',$merged);
		my $size = @ merged_genes;
		if ($size <= 1) {print STDERR "Gene '$Name' is tagged with 'merged'. At least two genes are expected as value for 'merged', but found is only ", Dumper(\@merged_genes),"\n"}
		my $gene_present='';
		foreach (@merged_genes){
			test_names($Name,'',$_,'merged');
			if ($Name eq $_){$gene_present=1};
		}
		unless ($gene_present){print STDERR "Gene '$Name' is not present in values for its own 'merged' tag. Correct and rerun\n"}
	}

##putpair and pair
	if ($pair){
		test_names($Name,'',$pair,'paired');
		unless ($pair eq $paired_gene->{_gsf_tag_hash}{Name}[0]){
			my $paired_pair = $paired_gene->{_gsf_tag_hash}{pair}[0];
			my $paired_name = $paired_gene->{_gsf_tag_hash}{Name}[0];
			print STDERR "Gene '$Name' and '$pair' are tagged as paired, but reciprocal tag in '$paired_name' is '$paired_pair'. Correct and rerun.\n"
		}
	}
	if ($putpair){
		test_names($Name,'',$putpair,'paired');	
               unless ($putpair eq $paired_gene->{_gsf_tag_hash}{Name}[0]){
                        my $paired_pair = $paired_gene->{_gsf_tag_hash}{putpair}[0];
                        my $paired_name = $paired_gene->{_gsf_tag_hash}{Name}[0];
                        print STDERR "Gene '$Name' and '$putpair' are tagged as paired, but reciprocal tag in '$paired_name' is '$paired_pair'. Correct and rerun.\n"
                }

	}

##check values for all tags that have an expected value of '1'.
	my %presence_tags=(truncated=>$truncated,
		noevidence=>$noevidence,
		fusion=>$fusion,
		corbound=>$corbound,
		cortrans=>$cortrans,
		misassembly=>$misassembly,
		delete=>$delete);
	foreach my $key (keys %presence_tags){
		if (defined $presence_tags{$key}){unless ($presence_tags{$key}==1){print STDERR "Expected value for '$key' is '1'. Correct this in gene '$Name' and rerun.\n"}}
	}
##remove empty tags again
	foreach my $key (keys %{$feature->{_gsf_tag_hash}}){
		$feature = remove_empty($feature,$key);	
	}
}

=head2 REMOVE_EMPTY

Removes empty tags

=cut

sub remove_empty {
	my $feature_curr=shift @_;
	my $tag=shift @_;

	unless ($feature_curr->{_gsf_tag_hash}{$tag}[0]){
		delete $feature_curr->{_gsf_tag_hash}{$tag};
	}
	return $feature_curr;
}


=head2 REPLACE

Compare and replace maker genes with corresponding reannotated genes from webapollo and report a final list of genes.

=cut


sub replace {
	my $maker_genes=shift @_;
	my %maker_genes=%$maker_genes;
	my $web_genes=shift @_;
	my %web_genes=%$web_genes;
	my %replaced_genes=();
	my %deleted=();
		foreach my $reanno_gene (keys %web_genes){
##gene neither merged nor fused		
			if ($web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{Name}[0]=~/^\d+\|G\d+$/ &&! $web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{merged}[0]){
				delete $maker_genes{$reanno_gene};
				$maker_genes{$reanno_gene}=$web_genes{$reanno_gene};	
			} 
##merged genes
			elsif ($web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{merged}[0]) {
				my @merged=split(' ',$web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{merged}[0]);
				foreach (@merged) {
					$_=~/^(\d+\|G\d+)/;
					unless ($deleted{$1}){
						delete $maker_genes{$1} or die "Trying to replace some maker genes with the merged gene '$web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{Name}[0]'. Tag 'merged' contains '$_', but '$1' is not in maker gff file. Correct and replace. \n";
						$deleted{$1}=1;
					}
				}
				$maker_genes{$reanno_gene}=$web_genes{$reanno_gene};
						
			}
##fusion genes
			elsif ($web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{fusion}[0]){
				my $Name=$web_genes{$reanno_gene}[0]->{_gsf_tag_hash}{Name}[0];
				$Name=~/^(\d+\|G\d+)\.\d+/;
				unless ($deleted{$1}){
					delete $maker_genes{$1};
					$deleted{$1}=1;
				}
				$maker_genes{$reanno_gene}=$web_genes{$reanno_gene};
			}
			else {
				print STDERR "Gene 'reanno_gene' does not fit into any replacement sceme. Contact developer.\n"
			}
##remove genes that are flagged with 'delete'
			if (exists $maker_genes{$reanno_gene}[0]->{_gsf_tag_hash}{delete}){
				delete $maker_genes{$reanno_gene};
				next;
			}
##remove empty tags
        		foreach my $key (keys $maker_genes{$reanno_gene}[0]->{_gsf_tag_hash}){
				$maker_genes{$reanno_gene}[0] = remove_empty($maker_genes{$reanno_gene}[0],$key);
        		}

		}			
	return %maker_genes;	
}
