use strict;
use warnings;

my $txt=$ARGV[0];
my $all=$ARGV[1];
my $out=$ARGV[2];
my %gene;
my @private;
my $ID=1;
open OUT,">$out" or die "cannot open the file $out\n";
open TXT,"$txt"or die "cannot open the file $txt\n";
while(my $line=<TXT>){
	chomp $line;
	my @genes = (split/\t/,$line);
	my @newgene = @genes[1..$#genes];
	my $gene_s=join("\t", @newgene);
	printf OUT ("OG%06d:\t%s\n",$ID,$gene_s);
	$ID++;
	foreach(@genes){
		$gene{$_}=1;
	}
}
close TXT;

open ALL,"$all" or die "cannot open the file $all\n";
while(my $line=<ALL>){
        chomp $line;
	if(!exists $gene{$line}){
		printf OUT ("OG%06d:\t%s\n",$ID,$line);
		$ID++;
	}
}
close ALL;

close OUT;
