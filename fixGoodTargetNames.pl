#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Bio::SeqIO;

my $help = 0;
my $in;
my $out;


GetOptions ("in=s"     => \$in,
            "out=s"         => \$out,
            "help|man"      => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";


if (!$in or !$out or $help) {
    die "Must supply --in and --out.\n";    
}

# After the first round of ARC, the target names in the FASTA file look like this:
# >CTSandBTS5iter_:_contig00003|E19A4|OPA_:_Contig002_contig00003|E19A4|OPA

# We want to revise those to look something like this:
# >contig00003--E19A4--OPA---CTSandBTS5iter
#
# So we'll apped the names with ---CTSandBTS5iterGoodTarget, and replace the vertical bars with --.
# First need to harvest the target name between the _:_(.*)_:_

my $seqIn = Bio::SeqIO->new(-file => $in,
                            -format => 'fasta');

my $seqOut = Bio::SeqIO->new(-file => ">$out",
                             -format => 'fasta');

while (my $seq = $seqIn->next_seq()) {
    if ($seq->display_id() =~ /\_\:\_(.*)\_\:\_/) {
        # $1 is the target name
        my $newID = $1 . "---CTSandBTS5iterGoodTarget";
        $newID =~ s/\|/\-\-/g;
        $seq->display_id($newID);
    } else {
        die "Target didn't have the _:_(.*)_:_ pattern in its name!\n";
    }
    $seqOut->write_seq($seq);
}
