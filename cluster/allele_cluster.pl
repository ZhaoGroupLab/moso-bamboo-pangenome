use strict;
use warnings;

# Create a hash to store the genes and their corresponding sets
my %gene_sets;
my $list=$ARGV[0];
my $out=$ARGV[1];

my $log=$ARGV[2];

open LI,"$list"or die "cannot open the file $list\n";

# Loop through each input file
my $n=1;
open LOG,">$log";
print LOG "start reading list $list\n";
close LOG;
while( my $file= <LI>) {
    chomp $file;
    # Open the input file
    open(my $fh, '<', $file) or die "Can't open $file: $!";
    open LOG,">>$log";
    print LOG (localtime) . " open the $n file:$file\n";
    close LOG;
    # Loop through each line of the input file
    while (my $line = <$fh>) {
        # Remove newline character
        chomp($line);
	next if($line=~/^#/);
        # Split the line into two genes
        my ($gene1, $gene2) = (split/\t/, $line)[0,1];
        # Check if either gene is already in a gene set
        my $gene1_set = $gene_sets{$gene1};
        my $gene2_set = $gene_sets{$gene2};

        # If both genes are already in a set, merge the sets
        if ($gene1_set and $gene2_set) {
            # If the sets are the same, do nothing
            if ($gene1_set eq $gene2_set) {
                next;
            }
            # Otherwise, merge the sets by adding all genes from one set to the other
            else {
                my @genes_to_move = grep { $_ ne $gene1_set } keys %gene_sets;
                foreach my $gene (@genes_to_move) {
                    if ($gene_sets{$gene} eq $gene1_set) {
                        $gene_sets{$gene} = $gene2_set;
                    }
                }
            }
        }
        # If one gene is already in a set, add the other gene to the set
        elsif ($gene1_set) {
            $gene_sets{$gene2} = $gene1_set;
        }
        elsif ($gene2_set) {
            $gene_sets{$gene1} = $gene2_set;
        }
        # If neither gene is in a set, create a new set and add both genes to it
        else {
            my $set_id = scalar(keys %gene_sets) + 1;
            $gene_sets{$gene1} = $set_id;
            $gene_sets{$gene2} = $set_id;
        }
    }

    # Close the input file
    close($fh);
    $n++;
}

close LI;
# Create a hash to store the gene sets and their corresponding genes
my %gene_set_genes;

# Loop through each gene and add it to its corresponding set
foreach my $gene (keys %gene_sets) {
    my $gene_set = $gene_sets{$gene};
    push(@{$gene_set_genes{$gene_set}}, $gene);
}

# Print the output file

open OUT,">$out";
foreach my $gene_set (sort {$a <=> $b} keys %gene_set_genes) {
    my $genes = join("\t", @{$gene_set_genes{$gene_set}});
    printf OUT ("OG%08d:\t%s\n",$gene_set,$genes);

}
close OUT;

