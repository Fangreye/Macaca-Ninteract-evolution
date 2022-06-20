#!/usr/bin/env perl
use strict;
use warnings;


# This program will read in two files. The first contains the coordinates of
# all N-mt interact genes, their acronyms, and whether or not (1 or 0) they interact
# directly with mt genes.

# The other file is a file with pi in windows with coordinates.  First
# the mean pi of Ninteract (1) and non-interact (0) windows will be calculated.
# Then permutations will be performed where the difference between these categories is 
# recalculated after the interaction is permuted n times. This will allow a p value of the
# observed pi to be estimated.

# this program deliverately ignores all genes in chrX.  A separate script will be generated
# that is only for chrX

# execute like this:
# ./All_N_mt_allinteract_pi_permutation.pl FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez_orientation.txt fst_bor/nem_mau_windowstats.concat

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $outputfile = $ARGV[2];

my @windowsites;
my @pi_values;
my $sumsites=0;
my $counter=0;
my @temp;
my $y;
my $x;
my %OXPHOS;

# first open up the OXPHOS gene info (# FINAL_OXPHOS_ARP2_MRP_MTREPLICATION_allinteractexceptC2_andallother_genez_orientation.txt)
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(($temp[0] ne 'gene')&&($temp[0] ne 'chrX')){ # deliberately ignores chrX
		if($temp[6] eq '+'){ # the gene is in the forward orientation
			$OXPHOS{$temp[0]."_".$temp[3]."_".$temp[4]}{"gene"} = $temp[9];
			$OXPHOS{$temp[0]."_".$temp[3]."_".$temp[4]}{"complex"} = $temp[9];
			$OXPHOS{$temp[0]."_".$temp[3]."_".$temp[4]}{"mt_interact"} = $temp[10];
		}
		elsif($temp[6] eq '-'){ # the gene is in the reverse orientation
			$OXPHOS{$temp[0]."_".$temp[4]."_".$temp[3]}{"gene"} = $temp[9];
			$OXPHOS{$temp[0]."_".$temp[4]."_".$temp[3]}{"complex"} = $temp[9];
			$OXPHOS{$temp[0]."_".$temp[4]."_".$temp[3]}{"mt_interact"} = $temp[10];
		}
		else{
			print "something wrong with gene orientation $line\n";
		}
	}	
}	
close DATAINPUT;

# this will print out a file that I may use for plotting later

#my $outputfile = $inputfile2."_pi__density.txt"; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

# now open up the pi data
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}

my @temp1;
my $N_interact_window=0;
my $gene_containing_window=0;
my $number_of_genes_in_this_window=0;
my $number_of_Ninteract_genes_in_this_window=0;
my $Ninteract_acronym="";
my $number_of_Ninteract_genes_spanning_a_window=0;

print OUTFILE "chr\tpos\tpi\tcontainsgenes\tcontainsNinteractgenez\tnumber_of_genes\tnumber_of_Ninteractgenez\tNinteract_acronym\n";

while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split(',',$line);
		$N_interact_window=0;
		$gene_containing_window=0;
		$number_of_genes_in_this_window=0;
		$number_of_Ninteract_genes_in_this_window=0;
		$Ninteract_acronym="-";
		# cycle through each gene
		foreach my $key (keys %OXPHOS){
			@temp1=split('_',$key);
			# check if this window contains one or more N_mt genes
			if(($temp1[0] eq $temp[0])&&($temp1[1] >= $temp[1])&&($temp1[1] <= $temp[2])){
					$gene_containing_window=1; # only count the start of genes in windows
					$number_of_genes_in_this_window+=1; # only count the start of genes in windows
					$OXPHOS{$key}{"start_pi"} = $temp[5];
					#$OXPHOS{$key}{"start_pi_sites"} = $temp[4];
					if($OXPHOS{$key}{"mt_interact"} == 1){
						$N_interact_window=1;
						$number_of_Ninteract_genes_in_this_window+=1; # only count the start of genes in windows
						if($number_of_Ninteract_genes_in_this_window == 1){
							$Ninteract_acronym=$OXPHOS{$key}{"gene"};
						}	
						else{
							$Ninteract_acronym=$Ninteract_acronym.",".$OXPHOS{$key}{"gene"};
						}	
					}	
			} # start is in block
		}
		if($temp[0] ne 'scaffold'){
			print OUTFILE $temp[0],"\t",$temp[1],"\t",$temp[5],"\t",$gene_containing_window,"\t",$N_interact_window,
			"\t",$number_of_genes_in_this_window,"\t",$number_of_Ninteract_genes_in_this_window,"\t",$Ninteract_acronym,"\n";
		}	
}

close OUTFILE;
close DATAINPUT2;

# Now the OXPHOS hash has pi and coordinates of all blocks that have genes
my @pi_for_perms; # this has only the mtinteractors and all genes
my $pi_associated=0; # all OXPHOS,MRP, ARP2 genes
my $pi_non_associated=0;
my $n_pi_associated=0; # this also works for only_N_mt
my $n_pi_non_associated=0; # this includes all non associated genes, including those anywhere 
my $N_interact_counter=0;

# now calculate the average pi for associated and non-associated OXPHOS genes
foreach my $key (keys %OXPHOS){
	if((exists($OXPHOS{$key}{"start_pi"})) &&
		($OXPHOS{$key}{"start_pi"} ne 'nan')){

		if($OXPHOS{$key}{"mt_interact"} == 1){
			$N_interact_counter+=1;
			$pi_associated += $OXPHOS{$key}{"start_pi"};
			$n_pi_associated += 1;	
		}
		elsif($OXPHOS{$key}{"mt_interact"} == 0){
			$pi_non_associated += $OXPHOS{$key}{"start_pi"};
			$n_pi_non_associated += 1;
		}
		push(@pi_for_perms,$OXPHOS{$key}{"start_pi"});
	}
}


# now report values
print "Number of associated blocks ",$n_pi_associated,"\n";
print "Number of N_interact genes ",$N_interact_counter,"\n";
print "Mean pi all associated N_mt genes ",$pi_associated/$n_pi_associated,"\n";
print "Mean pi non-associated all genes ",$pi_non_associated/$n_pi_non_associated,"\n";
#print "Mean pi non-associated only N_mt ",$pi_non_associated_only_N_mt/$n_pi_non_associated_only_N_mt,"\n";


#################
# ALL COMPLEXES perms including all genes
#################

# calculate test statistic for all genes
my $test_stat = ($pi_associated/$n_pi_associated) - ($pi_non_associated/$n_pi_non_associated);

# do permutations for all_complex_mt_interact vs all_other_genes
# first make an array that will be shuffled with the same number of 1s and 0s as the $OXPHOS{$key}{"complex"} variable
my @associated_or_not_array = (('1') x $n_pi_associated, ('0') x $n_pi_non_associated);

my $perms=1000;
my $pi_associated_perm=0;
my $pi_not_associated_perm=0;
my $n_pi_associated_perm=0;
my $n_pi_not_associated_perm=0;

# first check if the length of the permutation array is the same as the pi array
# for this analysis
if($#associated_or_not_array ne $#pi_for_perms){
	print "Problem: length associated_or_not_array ",$#associated_or_not_array,
	" length pi array ",$#pi_for_perms,"\n";
}

my @perm_diffs;
@perm_diffs={};
for ($y = 0 ; $y < $perms; $y++ ) {
	fisher_yates_shuffle( \@associated_or_not_array );    # permutes @array in place
	$pi_associated_perm=0;
	$pi_not_associated_perm=0;
	$n_pi_associated_perm=0;
	$n_pi_not_associated_perm=0;
	for ($x = 0 ; $x <= $#associated_or_not_array; $x++ ) {
		if($associated_or_not_array[$x] == 1){
			$pi_associated_perm+= $pi_for_perms[$x];
			$n_pi_associated_perm +=1;
		}
		elsif($associated_or_not_array[$x] == 0){	
			$pi_not_associated_perm+= $pi_for_perms[$x];
			$n_pi_not_associated_perm +=1;
		}
	}
	push(@perm_diffs,($pi_associated_perm/$n_pi_associated_perm) - ($pi_not_associated_perm/$n_pi_not_associated_perm));	
}

my @perm_diffs_sorted = sort { $a <=> $b } @perm_diffs;
my $switch=0;
my $pval=0;
# now figure out where the test stat is
for ($y = 0 ; $y <= $#perm_diffs_sorted; $y++ ) {
	if(($test_stat <= $perm_diffs_sorted[$y])&&($switch==0)){
		$pval=$counter;
		$switch = 1;
	}
	$counter+=1;
}	

#print "@perm_diffs_sorted\n";
print "Test stat for test including all genes (if positive, then not significant):",$test_stat,"\n";
print "P = ",($pval/$perms),"\n";




# fisher_yates_shuffle( \@array ) : 
    # generate a random permutation of @array in place
    sub fisher_yates_shuffle {
        my $array = shift;
        my $i;
        for ($i = @$array; --$i; ) {
            my $j = int rand ($i+1);
            next if $i == $j;
            @$array[$i,$j] = @$array[$j,$i];
        }
    }


