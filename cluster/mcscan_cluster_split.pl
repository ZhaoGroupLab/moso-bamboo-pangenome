use strict;
use warnings;

# 读取A文件，构建基因对哈希表
our %gene_pairs;
my $anchors=$ARGV[0];
my $mcscan=$ARGV[1];
my $out=$ARGV[2];
my @GS;
open(my $file_A, "<", "$anchors") or die "Cannot open A.txt: $!";
our $count = 0;  # 计数器，用于生成编号
while (my $line = <$file_A>) {
    chomp $line;
    my ($gene1, $gene2, $similarity, $distance) = split("\t", $line);
    my $gene_pair;
    if($gene1 gt $gene2){
        $gene_pair="$gene1-$gene2";
    }else{
        $gene_pair="$gene2-$gene1";
    }
    $gene_pairs{$gene_pair} = [$similarity, $distance];
}
close($file_A);

# 读取B文件，处理每一行基因集合

open(my $file_B, "<", "$mcscan") or die "Cannot open B.txt: $!";
open OUT,">$out";
my %gene_sets;
while (my $line = <$file_B>) {
    chomp $line;
    my @genes = split("\t", $line);
    shift @genes;  # 忽略第一列
    my $genes_num = scalar @genes;
    if($genes_num == 1){
        $count++;
        my $genes_list = join("", @genes);
        printf OUT ("GS%07d:\t%s\n", $count, $genes_list);
    }else{
    my %gene_sets_part=split_set(@genes);
    foreach my $id(keys %gene_sets_part){
        my @sorted_genes = (sort {$a cmp $b} @{$gene_sets_part{$id}});
        my $genes_list = join("\t", @sorted_genes);
        if($genes_list){
            $count++;
            printf OUT ("GS%07d:\t%s\n", $count, $genes_list);
        }
    }
    }
}

close OUT;
close($file_B);
   




sub split_set{
    my @genes=@_;

    my %gene_sets;
    my %gene_set_genes;
    my %gene_pairs_has;
    my $in_count=0;
    foreach my $gene1 (@genes){
        foreach my $gene2 (@genes){
            my $gene_pair="$gene1-$gene2";
            if(exists $gene_pairs{$gene_pair}){
                $gene_pairs_has{$gene_pair}=$gene_pairs{$gene_pair};
            }
        }
    }
    my @sorted_pairs = sort {
        $gene_pairs{$a}[1] <=> $gene_pairs{$b}[1] || # 按 distance 正序排序
        $gene_pairs{$b}[0] <=> $gene_pairs{$a}[0]    # 如果 distance 相同，则按 similarity 倒序排序
    } keys %gene_pairs_has;

foreach my $pair (@sorted_pairs) {

    my ($gene1,$gene2) = (split /-/, $pair)[0, 1];
    my ($similarity, $distance) = @{$gene_pairs_has{$pair}};
    #print "$gene1, $gene2,$similarity, $distance\n";
    
    my $gene1_set = $gene_sets{$gene1};
    my $gene2_set = $gene_sets{$gene2};

    if($gene1_set && $gene2_set){



        #print "1:$gene1, $gene2,$similarity, $distance,$gene_set_genes{$gene1_set},$gene_set_genes{$gene2_set}\n";
        if($gene1_set eq $gene2_set){
            next;
        }else{
            #print "---$gene1_set---$gene2_set---\n";
            my $gene1_set_list=join("\t",@{$gene_set_genes{$gene1_set}});
            my $gene2_set_list=join("\t",@{$gene_set_genes{$gene2_set}});
            #print "$gene1_set_list   vs    $gene2_set_list\n";
            # my $gene1_set_list=join("\t",@{$gene_set_genes{$gene1_set}});
            # my $gene2_set_list=join("\t",@{$gene_set_genes{$gene2_set}});
            # print "$gene1_set_list---vs---$gene2_set_list\n";
            #     print "ref(\$gene_set_genes{\$gene1_set}): " . ref($gene_set_genes{$gene1_set}) . "\n";
            #     print "ref(\$gene1_set): " . ref($gene1_set) . "\n";
            #     print "ref(\$gene_set_genes{\$gene2_set}): " . ref($gene_set_genes{$gene2_set}) . "\n";
            #     print "ref(\$gene2_set): " . ref($gene2_set) . "\n";

                my @genes1=@{$gene_set_genes{$gene1_set}};
                my @genes2=@{$gene_set_genes{$gene2_set}};


                my $genes1_num = scalar @genes1;
                my $genes2_num = scalar @genes2;
                 #print "$genes1_num\t$genes2_num\n";

            unless(same_groups(\@{$gene_set_genes{$gene1_set}}, \@{$gene_set_genes{$gene2_set}})){
                    #print "before:$gene_sets{$gene1}\t$gene_sets{$gene2};\n";
                    my $new_gene_set="$gene1_set-$gene2_set";
                    foreach my $genein1 (@{$gene_set_genes{$gene1_set}}){
                        $gene_sets{$genein1} = $new_gene_set;
                    }
                    foreach my $genein2 (@{$gene_set_genes{$gene2_set}}){
                        $gene_sets{$genein2} = $new_gene_set;
                    }
                    push @{$gene_set_genes{$new_gene_set}}, @{$gene_set_genes{$gene1_set}};
                    push @{$gene_set_genes{$new_gene_set}}, @{$gene_set_genes{$gene2_set}};
                    delete $gene_set_genes{$gene1_set};
                    delete $gene_set_genes{$gene2_set};
                    #print "after:$gene_sets{$gene1}\t$gene_sets{$gene2};\n";

            }
        }
    }elsif($gene1_set && !$gene2_set) {
        #print "2:$gene1, $gene2,$similarity, $distance,$gene_set_genes{$gene1_set},$gene_set_genes{$gene2_set}\n";
        unless (have_same_group($gene2, @{$gene_set_genes{$gene1_set}})) {
            $gene_sets{$gene2} = $gene1_set;
            push(@{$gene_set_genes{$gene1_set}}, $gene2);
        }else{
            $in_count++;
            my $set_id="$gene2-$in_count";
            $gene_sets{$gene2} = $set_id;
            push(@{$gene_set_genes{$set_id}}, $gene2);
        }
    }
    elsif (!$gene1_set && $gene2_set ) {
        #print "3:$gene1, $gene2,$similarity, $distance,$gene_set_genes{$gene1_set},$gene_set_genes{$gene2_set}\n";
        unless (have_same_group($gene1, @{$gene_set_genes{$gene2_set}})) {
            $gene_sets{$gene1} = $gene2_set;
            push(@{$gene_set_genes{$gene2_set}}, $gene1);
        }else{
            $in_count++;
            my $set_id="$gene1-$in_count";
            $gene_sets{$gene1} = $set_id;
            push(@{$gene_set_genes{$set_id}}, $gene1);
        }
    }
    # If neither gene is in a set, create a new set and add both genes to it
    elsif (!$gene2_set && !$gene1_set) {
        #print "4:$gene1, $gene2,$similarity, $distance,$gene_set_genes{$gene1_set},$gene_set_genes{$gene2_set}\n";
        $in_count++;
        my $set_id="$pair-$in_count";
        $gene_sets{$gene1} = $set_id;
        push(@{$gene_set_genes{$set_id}}, $gene1);
        $gene_sets{$gene2} = $set_id;
        push(@{$gene_set_genes{$set_id}}, $gene2);
    }
    
}

    return %gene_set_genes;
}


sub have_same_group {
    my ($gene, @genes) = @_;
    
    my $jug=0;
    my $gene_group = $gene =~ /^(.+?)Gene\d+$/ ? $1 : "";
    #print "$gene\t$gene_group\n";
    foreach my $existing_gene (@genes) {
        my $group = $existing_gene =~ /^(.+?)Gene\d+$/ ? $1 : "";
        #print "have_same_group:$gene\t$gene_group\t$existing_gene\t$group\n";

        if ($gene_group eq $group) { 
            $jug++;
        }
    }
    #print "x,$jug\n";
    return $jug;
}

sub have_same_grouptry {
    my ($gene, @genes) = @_;
    
    my $jug=0;
    my $gene_group = $gene =~ /^(.+?)Gene\d+$/ ? $1 : "";
    #print "$gene\t$gene_group\n";
    foreach my $existing_gene (@genes) {
        my $group = $existing_gene =~ /^(.+?)Gene\d+$/ ? $1 : "";
        #print "have_same_group:$gene\t$gene_group\t$existing_gene\t$group\n";

        if ($gene_group eq $group) { 
            $jug++;
        }
    }
    #print "x,$jug\n";
    return $jug;
}


sub same_groups {
    my ($genes1, $genes2) = @_;
    my $jug=0;
    my $genes1_num = scalar @$genes1;
    my $genes2_num = scalar @$genes2;
    #print "$genes1_num\t$genes2_num\n";
    foreach my $gene1 (@$genes1) {

        $jug+=(have_same_grouptry($gene1,@$genes2));

    }
    #print "y,$jug\n";
    return $jug;
}
sub split_gene_sets {
    my ($gene_set, $gene_pairs) = @_;
    my @allele_sets;

    my @sorted_genes = sort {
        $gene_pairs->{$a}->{$b}->[1] <=> $gene_pairs->{$b}->{$a}->[1]
            || $gene_pairs->{$a}->{$b}->[0] <=> $gene_pairs->{$b}->{$a}->[0]
            || $a cmp $b
    } keys %$gene_pairs;

    foreach my $gene_pair (@sorted_genes) {
        my ($gene1, $gene2) = split(':', $gene_pair);

        my $added = 0;
        foreach my $allele_set (@allele_sets) {
            if (gene_exists_in_set($gene1, $allele_set) && !gene_exists_in_set($gene2, $allele_set)) {
                if (!have_same_group($gene2, @$allele_set)) {
                    push @$allele_set, $gene2;
                    $added = 1;
                    last;
                }
            } elsif (gene_exists_in_set($gene2, $allele_set) && !gene_exists_in_set($gene1, $allele_set)) {
                if (!have_same_group($gene1, @$allele_set)) {
                    push @$allele_set, $gene1;
                    $added = 1;
                    last;
                }
            }
        }

        unless ($added) {
            push @allele_sets, [$gene1, $gene2];
        }
    }
    
    return @allele_sets;
}

sub gene_exists_in_set {
    my ($gene, $allele_set) = @_;

    foreach my $existing_gene (@$allele_set) {
        if ($gene eq $existing_gene) {
            return 1;
        }
    }

    return 0;
}


