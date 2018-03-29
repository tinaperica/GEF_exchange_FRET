#!/usr/bin/perl -w
use strict;

# prefixes (experiments so far)

my ($arg, $control) = @ARGV;
my ($InFile, $IndexFile, $OutFile, $prefix);
if ($arg =~ /^[0-9]+$/) {   ### this is if $arg is just the date e.g. 20160721
    $InFile = $arg."_TP_FRET_kinetics.asc";
    $IndexFile = $arg."_TP_FRET_kinetics_INDEX.txt";
    $OutFile = $arg."_TP_FRET_kinetics_parsed.txt";
    $prefix = $arg; 
} else {
    $InFile = $arg.".asc";
    $IndexFile = $arg."_INDEX.txt";
    $OutFile = $arg."_parsed.txt";
    $prefix = substr($arg, 0, 8);
}
if ($control) {
    print STDERR $control;
    if ($arg =~ /^[0-9]+$/) {   ### this is if $arg is just the date e.g. 20160721
        $InFile = $arg."_TP_FRET_control_kinetics.asc";
        $IndexFile = $arg."_TP_FRET_control_kinetics_INDEX.txt";
        $OutFile = $arg."_TP_FRET_control_kinetics_parsed.txt";
        $prefix = $arg; 
    } else {
        $InFile = $arg.".asc";
        $IndexFile = $arg."_INDEX.txt";
        $OutFile = $arg."_parsed.txt";
        $prefix = substr($arg, 0, 8);
    }
}
open(IND, $IndexFile) || die "Cannot open $IndexFile: $!\n";
my %index;
while (my $line = <IND>) {
    chomp $line;
    if (! $control) {
        my ($well, $protein, $concentration, $GEF_conc, $time_cutoff) = split(",", $line);
        if ($well =~ /([A-Z])([0-9]+)/) {
            my $row = $1; my $col = $2;
            $index{$row}->{$col}->{'protein'} = $protein;
            $index{$row}->{$col}->{'conc'} = $concentration;
            $index{$row}->{$col}->{'GEF'} = $GEF_conc;
            $index{$row}->{$col}->{'time'} = $time_cutoff;
        }
    } elsif ($control) {
        my ($well, $protein, $concentration, $nucleotide, $GEF_conc) = split(",", $line);
        if ($well =~ /([A-Z])([0-9]+)/) {
            my $row = $1; my $col = $2;
            $index{$row}->{$col}->{'protein'} = $protein;
            $index{$row}->{$col}->{'conc'} = $concentration;
            $index{$row}->{$col}->{'nucleotide'} = $nucleotide;
            $index{$row}->{$col}->{'GEF'} = $GEF_conc;
        }
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
#print OUT "time\trow\tcol\tprotein\tGsp1_concentration\tGEF_concentration\tfluorescence\n";
foreach my $time (sort {$a <=> $b} keys %kinetic_data) {
    foreach my $row (sort keys %{$kinetic_data{$time}}) {
        foreach my $col (sort {$a <=> $b} keys %{$kinetic_data{$time}->{$row}}) {
            print OUT $prefix, "\t", $time, "\t", $row, "\t", $col, "\t", $index{$row}->{$col}->{'protein'}, "\t", $index{$row}->{$col}->{'conc'}, "\t",
            $index{$row}->{$col}->{'GEF'}, "\t", $kinetic_data{$time}->{$row}->{$col};
            if ($control) {
                print OUT "\t", $index{$row}->{$col}->{'nucleotide'}, "\n";
            } else {
                print OUT "\t", $index{$row}->{$col}->{'time'}, "\n";
            }
        }
    }
}

system("cat *_parsed.txt > GEF_FRET_kinetics.txt");

