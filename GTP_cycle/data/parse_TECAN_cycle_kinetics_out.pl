#!/usr/bin/perl -w
use strict;

# prefixes (experiments so far)

my ($arg) = @ARGV;
my ($InFile, $IndexFile, $OutFile, $prefix);
if ($arg =~ /^[0-9]+$/) {   ### this is if $arg is just the date e.g. 20160721
    $InFile = $arg."_TP_Cycle_kinetics.asc";
    $IndexFile = $arg."_TP_Cycle_kinetics_INDEX.txt";
    $OutFile = $arg."_TP_Cycle_kinetics_parsed.txt";
    $prefix = $arg; 
} else {
    $InFile = $arg.".asc";
    $IndexFile = $arg."_INDEX.txt";
    $OutFile = $arg."_parsed.txt";
    $prefix = substr($arg, 0, 8);
}
open(IND, $IndexFile) || die "Cannot open $IndexFile: $!\n";
my %index;
my $head = <IND>;
while (my $line = <IND>) {
    chomp $line;
    my ($well, $protein, $protein_conc, $GAP_conc, $GEF_conc, $MOG1_conc) = split("\t", $line);
    if ($well =~ /([A-Z])([0-9]+)/) {
        my $row = $1; my $col = $2;
        $index{$row}->{$col}->{'protein'} = $protein;
        $index{$row}->{$col}->{'protein_conc'} = $protein_conc;
        $index{$row}->{$col}->{'GAP_conc'} = $GAP_conc;
        $index{$row}->{$col}->{'GEF_conc'} = $GEF_conc;
        $index{$row}->{$col}->{'MOG1_conc'} = $MOG1_conc;
    }
}
open(IN, $InFile) || die "Cannot open $InFile: $\n";
my ($time, %kinetic_data);
while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /([0-9]+) s/) {
        $time = $1;
    } elsif ($time && $line =~ /^([A-Z])\t(.+)/) {
        my $row = $1; my $values_line = $2;
        if (exists $index{$row}) {  
            my @values = split("\t", $values_line);
            for (my $col = 1; $col < scalar @values; $col++) {
                if ($values[$col -1]) {
                    my $val = $values[$col - 1];
                    if ($val eq "Overflow") {$val = "NA"};
                    $kinetic_data{$time}->{$row}->{$col} = $val;
                }
            }
        }
    }
}

open(OUT, ">$OutFile") || die "Cannot open $OutFile: $!\n";
#print OUT "time\trow\tcol\tprotein\tprotein_conc\tGAP_concentration\tGEF_concnetration\tMOG1_concentration\tfluorescence\n";
foreach my $time (sort {$a <=> $b} keys %kinetic_data) {
    foreach my $row (sort keys %{$kinetic_data{$time}}) {
        foreach my $col (sort {$a <=> $b} keys %{$kinetic_data{$time}->{$row}}) {
            print OUT $prefix, "\t", $time, "\t", $row, "\t", $col, "\t", $index{$row}->{$col}->{'protein'}, "\t",
            $index{$row}->{$col}->{'protein_conc'}, "\t", $index{$row}->{$col}->{'GAP_conc'}, "\t",
            $index{$row}->{$col}->{'GEF_conc'}, "\t", $index{$row}->{$col}->{'MOG1_conc'}, "\t", $kinetic_data{$time}->{$row}->{$col}, "\n";
        }
    }
}

system("cat *parsed.txt >cycle_kinetics.txt")

