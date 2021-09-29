#!/usr/bin/env perl
use strict;
use Data::Dumper;
 
my @exons;
my $gene_count = 0;
my $current_seq = "";
while(my $line = <STDIN>)
{
  if($line =~ m/^###/)
  {
    flush(\@exons);
    next;
  }
  my @fields = split(/\t/, $line);
  if($fields[2] eq "mRNA")
  {
    flush(\@exons);
  }
  elsif($fields[2] eq "exon")
  {
    if ($fields[0] ne $current_seq)
    {
      $current_seq = $fields[0];
      printf(">%s\n", $current_seq);
    }
    push(@exons, \@fields);
  }
}
flush();
 
sub flush
{
  my $num_exons = scalar(@exons);
  return if($num_exons == 0);
   
  my $group = sprintf("%s.%d", $exons[0]->[0], $gene_count);
  $gene_count++;
   
  if($num_exons == 1)
  {
    my($start, $end) = ($exons[0]->[3], $exons[0]->[4]);
    if($exons[0]->[6] eq "-")
    {
      ($start, $end) = ($exons[0]->[4], $exons[0]->[3]);
    }
    printf("Esngl\t%lu\t%lu\t%s\n", $start, $end, $group);
  }
  else  {
    my @fexons = ();
    if($exons[0][6] eq "-") {
	    @fexons = @exons;
            #printf STDERR "Reversing group ${group}\n";
	    #printf STDERR Dumper(\@fexons) . "\n";
    } else {
	@fexons = @exons;
    } 
    for(my $i = 0; $i < $num_exons; $i++) {
      my $exon_type = "Exon";
      if($i == 0) {
        $exon_type = "Einit";
      }
      elsif($i == $num_exons - 1) {
        $exon_type = "Eterm";
      }
       
      my($start, $end) = ($fexons[$i][3], $fexons[$i][4]);

      if($fexons[0][6] eq "-") {
        ($start, $end) = ($fexons[$i][4], $fexons[$i][3]);
      }
       
      printf("%s\t%lu\t%lu\t%s\n", $exon_type, $start, $end, $group);
    }
  }
  @exons = ();
}
