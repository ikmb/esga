#!/usr/bin/env perl

use strict;
use warnings;

# accept SNAP ("ZFF") format such as
# >chr1
# Einit   5343    5444    +       8.956   0       0       2       chr1-snap.1
# Exon    6401    6498    +       -0.203  0       2       1       chr1-snap.1
# Eterm   6727    6781    +       7.370   1       0       1       chr1-snap.1
# output standard GFF

my $current_id = '';
my $current_chr = '';
my $real_start = 0;
my $real_end = 0;
my @arr = ();

while(my $line = <>) {
    chomp($line);
    if(substr($line, 0, 1) eq ">") {
        $line =~ s/>//;
        $current_chr = $line;
    } else {
        my ($type, $start, $end, $strand, $score, $blank1, $blank2, $blank3, $id) = split('\t', $line);
        if($type eq "Einit") {
            if($current_id ne "") {
                print "$current_chr\tSNAP\tmRNA\t$real_start\t$real_end\t$score\t$strand\t$blank3\tID=$current_id\n";
            }
            foreach (@arr) {
                print "$_\n";
            }
            @arr = ();
            $current_id = $id;
            $real_start = $start;
        }
        push(@arr, "$current_chr\tSNAP\texon\t$start\t$end\t$score\t$strand\t$blank3\tParent=$current_id");
        $real_end = $end;
    }
}
