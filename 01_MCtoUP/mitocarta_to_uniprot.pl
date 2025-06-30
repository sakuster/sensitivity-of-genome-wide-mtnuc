#!/usr/bin/env perl

use strict;
use warnings;
use sloan;
use Bio::SearchIO;

#usage statement
my $usage = "\nUSAGE: perl $0 mitocarta_fasta uniprot_fasta blast_file > output_file\n\n";

#read file names from command line
my $mitocarta_file = shift or die ($usage);
my $uniprot_file = shift or die ($usage);
my $blast_file = shift or die ($usage);

#read in fasta files into hash data structures
my %mitocarta = fasta2hash($mitocarta_file);
my %uniprot = fasta2hash($uniprot_file);

#loop over uniprot fasta headers and store a count of each gene name.
#Name is assume to be after the second pipe (|) and followed by an underscore.
my %up_name_count_hash;
foreach (keys %uniprot){
	my @sh1 = split (/\|/, $_);
	my @sh2 = split (/\_/, $sh1[2]);
	++$up_name_count_hash{uc($sh2[0])};
}

#Initialize data structure that will store main results.
#hash of array.
#Key = mitocarta unique ID
#array element 0 = mitocarta gene name
#array element 1 = count of identical gene names in uniprot
#array element 2 = uniprot ID of top blast hit
#array element 3 = uniprot gene name of top blast hit
#array element 4 = blast % ID
#array element 5 = blast % ID without gaps
#array element 6 = blast % coverage of query sequence
my %mc_id_name_HoA; 

#loop over mitocarta fasta headers.
#Store the relationship between unique ID (1st element) and gene name (3rd element).
#Store count of uniprot genes with identical names (case insensitive)
#store lists of unique mitocarta names and sequences
my %unique_mc_name_hash;
my %unique_mc_seq_hash;
foreach (sort keys %mitocarta){
	my @sh = split (/\t/, $_);
	$mc_id_name_HoA{$sh[0]}[0] = $sh[2];
	if (exists $up_name_count_hash{$sh[2]}){
		$mc_id_name_HoA{$sh[0]}[1] = $up_name_count_hash{$sh[2]};
	}else{
		$mc_id_name_HoA{$sh[0]}[1] = 0;		
	}
	$unique_mc_name_hash{$sh[2]} = 1;
	$unique_mc_seq_hash{$mitocarta{$_}} = 1;
}

#print unique name and seq counts for mc fasta
my @unique_names = keys %unique_mc_name_hash;
my @unique_seqs = keys %unique_mc_seq_hash;

print STDERR "\nNumber of unique gene names in $mitocarta_file: ", scalar (@unique_names), "\n";
print STDERR "Number of unique sequences in $mitocarta_file: ", scalar (@unique_seqs), "\n\n";


#loop over blast output. For each mitocarta query, store uniprot ID and gene name of top blast hit, percent ID, and percent coverage.
#check for and report any queries that have the same mitocarta gene name but have different top blast hits
my %hit_name_hash; #key is mitocarta gene name; value is uniprot ID

my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => "$blast_file");
while( my $result_obj = $SearchIO_obj->next_result ) {
	my $query_name = $result_obj->query_name;
	unless (exists ($mc_id_name_HoA{$query_name})){
		print STDERR "\nWARNING: $query_name found in $blast_file but not in $mitocarta_file\n"
	}
	my $query_mc_gene = $mc_id_name_HoA{$query_name}[0];

	my $query_length = $result_obj->query_length;
	if ( my $hit_obj = $result_obj->next_hit ) {
		my $hit_name = $hit_obj->name;
		my @sn = split (/\|/, $hit_name);
		$mc_id_name_HoA{$query_name}[2] = $sn[1];
		my @sn2 = split (/\_/, $sn[2]);
		$mc_id_name_HoA{$query_name}[3] = $sn2[0];
		my $hsp_obj = $hit_obj->next_hsp;
		$mc_id_name_HoA{$query_name}[4] = $hsp_obj->frac_identical;
		$mc_id_name_HoA{$query_name}[5] = $hsp_obj->num_identical / ($hsp_obj->length('total') - $hsp_obj->gaps);
		$mc_id_name_HoA{$query_name}[6] = $hsp_obj->length('query') / $query_length;
		if (exists($hit_name_hash{$query_mc_gene})){
			$hit_name_hash{$query_mc_gene} eq $sn[1] or print STDERR "WARNING: $query_name has mitocarta gene name $query_mc_gene with top blast hit to $sn[1] but previous queries with that gene name hit $hit_name_hash{$query_mc_gene}.\n";
		}else{
			$hit_name_hash{$query_mc_gene} = $sn[1];
		}
	}else{
		$mc_id_name_HoA{$query_name}[2] = "NO_BLAST_HIT";
		$mc_id_name_HoA{$query_name}[3] = "NO_BLAST_HIT";
		$mc_id_name_HoA{$query_name}[4] = "NO_BLAST_HIT";
		$mc_id_name_HoA{$query_name}[5] = "NO_BLAST_HIT";
		$mc_id_name_HoA{$query_name}[6] = "NO_BLAST_HIT";
		if (exists($hit_name_hash{$query_mc_gene})){
			$hit_name_hash{$query_mc_gene} eq "NO_BLAST_HIT" or print STDERR "WARNING: $query_name has mitocarta gene name $query_mc_gene and returns no blast hit but previous queries with that gene name hit $hit_name_hash{$query_mc_gene}.\n";
		}else{
			$hit_name_hash{$query_mc_gene} = "NO_BLAST_HIT";
		}
	}
}

print "MitocartaID\tMitocartaGene\tUniprotNameMatches\tUniprotBlastID\tUniprotBlastName\tBlastIdentity\tBlastIdentityNoGaps\tBlastCoverage\n";

#loop over main data structure and print output
foreach (sort keys %mc_id_name_HoA){
	my @output_line = @{$mc_id_name_HoA{$_}};
	print "$_\t", join("\t", @output_line), "\n";
}
