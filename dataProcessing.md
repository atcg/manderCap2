**Mander Data Processing:**
===========================

Prerequisites:
--------------
  * Bioperl
  * spades
  * bowtie2
  * ARC
  * Standalone blast+
  * The Parallel::ForkManager module for Perl



Step 1:
-------
Filter adapter contamination and quality trim reads:<br>
`perl scriptBelow.pl --out qc_to_fastq_join --reads concatenatedReads/ --adatpers HSEM020_adapters/ --log runQC_to_fastq-join.log --threads 12 > runQC_to_fastq-join.stdoutanderr 2>&1`
```perl
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use Parallel::ForkManager;

my $help = 0;
my $outDir;
my $logFile;
my $readsDir;
my $adaptersDir;
my $threadsMax = 1;

GetOptions  ("out=s"    => \$outDir,
             "reads=s"  => \$readsDir,
             "adapters=s" => \$adaptersDir,
             "log=s"    => \$logFile,
             "threads=i"=> \$threadsMax,
             "help|man" => \$help) || pod2usage(2);

if (!$outDir or !$logFile or !$readsDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}

my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}


open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";

# First gather up all the file names of the files in the specified reads directory
my %readFilesHash;
my %sampleNamesHash;
chdir($readsDir);
opendir my $readsDirFH, "./";
my @readsFiles = readdir $readsDirFH or die "Couldn't readdir $readsDirFH: $!\n";
closedir $readsDirFH;
chdir $startingDir;
my $trimmomaticDir = $outDir . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}

# Now do sequence QC using Trimmomatic
{
my @trimmomaticCommands;
print $logFH "Generating trimmomatic commands on all read files.\n";
foreach my $file (@readsFiles) {
    if ($file =~ /(.*)_(.*)_L001_R1_001.fastq/) { # Only looking for the R1s here, because we only want one trimmomatic run per library (or in the case of too many reads for a single fastq.gz file, do a trimmomatic run for each pair of read files)
	
	# Example file names:
	#02F_0806_BTS_TACGGTTG-AACACGCT_L001_R1_001.fastq
	#02F_0806_BTS.adapters
        # $1 = sample name, something like 02F_0806_BTS
        # $2 = index sequence, something like TACGGTTG-AACACGCT
        print $logFH "--------------------------------------------------\n";
        print $logFH "Sequence file found: $1\_$2\_L001_R1\_001.fastq\n";
        my $R1File = $readsDir . "$1\_$2\_L001_R1\_001.fastq";
        print $logFH "R1 file: $readsDir" . "$1\_$2\_L001_R1_001.fastq\n";
        my $R2File = $readsDir . "$1\_$2\_L001_R2\_001.fastq";
        print $logFH "R2 file: $readsDir" . "$1\_$2\_L001_R2_001.fastq\n";
        my $readGroupName = $1;
        
        my $adaptersFile = $adaptersDir . $1 . ".adapters";
        print $logFH "Adapters file to be used for the sequence group: " . $adaptersFile . "\n";
        
        my $R1OutFilePaired = "$trimmomaticDir/$1\_R1_paired_trimmed.fastq";
        print $logFH "R1 paired trimmomatic output file: $trimmomaticDir/$1\_R1_paired_trimmed.fastq\n";
        my $R1OutFileSingles = "$trimmomaticDir/$1\_R1_singles_trimmed.fastq";
        print $logFH "R1 singles trimmomatic output file: $trimmomaticDir/$1\_R1_singles_trimmed.fastq\n";
        my $R2OutFilePaired = "$trimmomaticDir/$1\_R2_paired_trimmed.fastq";
        print $logFH "R2 paired trimmomatic output file: $trimmomaticDir/$1\_R2_paired_trimmed.fastq\n";
        my $R2OutFileSingles = "$trimmomaticDir/$1\_R2_singles_trimmed.fastq";
        print $logFH "R2 singles trimmomatic output file: $trimmomaticDir/$1\_R2_singles_trimmed.fastq\n";
        $sampleNamesHash{$readGroupName}{'R1_paired_trimmed'} = $R1OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R1_singles_trimmed'} = $R1OutFileSingles;
        $sampleNamesHash{$readGroupName}{'R2_paired_trimmed'} = $R2OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R2_singles_trimmed'} = $R2OutFileSingles;        
        push (@trimmomaticCommands, "java -Xmx8G -jar ~/bin/trimmomatic/trimmomatic-0.32.jar PE -threads 2 -phred33 $R1File $R2File $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adaptersFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");
    }
}
print $logFH "--------------------------------------------------\n";
print $logFH "Finished generating trimmomatic commands on all read files.\n\n\n";

print $logFH "Running all trimmomatic commands\n";
my $counter = 0;
my $forkManager = new Parallel::ForkManager($threadsMax);
foreach my $trimCommand (@trimmomaticCommands) {
    $counter++;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Trimmomatic command $counter: \n\t"; # Indent the next line to make it easier to find the commands in the text
    print $logFH $trimCommand . "\n";
    sleep 10;
    print "\n";
    $forkManager->start and next;
    print "\n";
    system("$trimCommand");
    print "Finished running the following:\n\t$trimCommand\n\n";
    $forkManager->finish;
}
$forkManager->wait_all_children;
print $logFH "--------------------------------------------------\n";
print $logFH "Finished running all trimmomatic commands\n";
print $logFH "--------------------------------------------------\n\n";
}

# Now we can iterate through the %sampleNamesHash that we populated earlier to do the following:
#  1. Find read groups that have both 001 and 002, and cat the 002s onto the end of the 001s
#  1b. Delete the 002 files
#  2. Run fastq-join on the output R1 and R2 paired files
#  3. Merge the R1 and R2 singles files together into one singles file, and add the fastq-join joined reads with the "true" single reads
my @filesToZip;
my @assemblyFiles;

my $forkManager2 = new Parallel::ForkManager(2); # We want to be able to run in parallel two biological samples for the gunzip, fastq-join, gzip, and assembly stages

foreach my $readGroup (sort keys %sampleNamesHash) {
    sleep 10;
    print $logFH "\n\nStarting to process $readGroup through the fastq-join and assembly stages\n\n";
    $forkManager2->start and next;
        print $logFH "\n--------------------------------------------------\n";

        my $fastqJoinOutDir = $outDir . "/fastq-join/";
        unless (-d $fastqJoinOutDir) {
            mkdir $fastqJoinOutDir;
        }
        chdir $fastqJoinOutDir;
        print $logFH "\n--------------------------------------------------\n";
        print $logFH "\nRunning fastq-join on $readGroup\n";
        print $logFH "\nfastq-join command used:\n\t";
        my $fastqJoinOutputName = $readGroup . ".%.fastq";
        print $logFH "\nfastq-join -v ' ' -m 8 ../trimmomatic/$readGroup\_R1_paired_trimmed.fastq ../trimmomatic/$readGroup\_R2_paired_trimmed.fastq -o $fastqJoinOutputName";
        system("fastq-join -v ' ' -m 8 ../trimmomatic/$readGroup\_R1_paired_trimmed.fastq ../trimmomatic/$readGroup\_R2_paired_trimmed.fastq -o $fastqJoinOutputName");
        print $logFH "\nFinished running fastq-join on $1\n";
        print $logFH "\n--------------------------------------------------\n\n";
        
        my $joinedFile = $readGroup . ".join.fastq";
        my $newR1file = $readGroup . ".un1.fastq";
        my $newR2file = $readGroup . ".un2.fastq";

        my $singlesFile1 = "../trimmomatic/$readGroup\_R1_singles_trimmed.fastq";
        my $singlesFile2 = "../trimmomatic/$readGroup\_R2_singles_trimmed.fastq";
        my $joinedAndSinglesFile = $readGroup . "_joined_and_both_singles.fastq";
        
        print $logFH "\n--------------------------------------------------\n";
        print $logFH "\nCombining fastq-joined reads with the singletons Trimmomatic made from R1 and R2.\n";
        print $logFH "\nConcatenation command:\n\t";
        print $logFH "\ncat $joinedFile $singlesFile1 $singlesFile2 > $joinedAndSinglesFile\n";
        system("cat $joinedFile $singlesFile1 $singlesFile2 > $joinedAndSinglesFile");
        print $logFH "\nRemoving the joined and singles files\n";
        print $logFH "\nFile removal command:\n\t";
        print $logFH "\nunlink($joinedFile)\n";
        unlink($joinedFile);
        print $logFH "\nFinished combining joined reads with singletons and removing the consitutuent files.\n";
        print $logFH "\n--------------------------------------------------\n\n";
        chdir $startingDir;

    $forkManager2->finish;

    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running all fastq-join and assembly commands\n";
    print $logFH "--------------------------------------------------\n\n";
}
$forkManager2->wait_all_children;


#Documentation
__END__

=head1 NAME

reads_to_assemblies.pl ##CHANGE

=head1 SYNOPSIS 

perl reads_to_assemblies.pl --out <directory_name> --log <file> --threads <integer> --reads <directory> --kmer <integer>

 Options:
   -out=s           Name of directory where all output will be saved
   -reads=s         Name of directory where all raw reads are
   -log=s           Log filename
   -threads=i       Maximum number of threads (default 1)
   -kmer=i          kmer value to use for Abyss assemblies
   -help|man        Prints out documentation

=head1 DESCRIPTION

This program takes raw read files in fastq.gz format and does sequence quality
trimming and adapter contamination removal. It then merges overlapping paired
end reads, and combines the single orphaned reads from the qc process with the
joined reads to create the "single-end" library for each sample.


=cut
```

Step 2:
-------
Combine all the reads for the CTS and BTS libraries into CTS and BTS samples:

```bash
#!/bin/bash
cat 01A_0801_CTS.un1.fastq 01B_0802_CTS.un1.fastq 01C_0803_CTS.un1.fastq 01D_0804_CTS.un1.fastq 01E_0805_CTS.un1.fastq 01F_0806_CTS.un1.fastq 01G_0911_DSN_CTS.un1.fastq 02B_0906_CTS.un1.fastq 02C_0907_CTS.un1.fastq 02D_0908_CTS.un1.fastq 02E_0909_CTS.un1.fastq 02F_0910_CTS.un1.fastq 02G_0911_CTS.un1.fastq 03A_0905_CTS.un1.fastq 03D_0812_CTS.un1.fastq 03E_0901_CTS.un1.fastq 03F_0902_CTS.un1.fastq 03G_0903_CTS.un1.fastq 03H_0904_CTS.un1.fastq > AllCTS.un1.fastq
cat 01A_0905_BTS.un1.fastq 01B_0906_BTS.un1.fastq 01C_0907_BTS.un1.fastq 01D_0908_BTS.un1.fastq 01E_0901_BTS.un1.fastq 01E_0909_BTS.un1.fastq 01F_0902_BTS.un1.fastq 01F_0910_DSN_BTS.un1.fastq 01G_0903_BTS.un1.fastq 01H_0904_BTS.un1.fastq 02A_0905_BTS.un1.fastq 02F_0806_BTS.un1.fastq 02G_0807_BTS.un1.fastq 02H_0808_BTS.un1.fastq 03A_0809_BTS.un1.fastq 03B_0810_BTS.un1.fastq 03C_0811_BTS.un1.fastq 03H_0912_BTS.un1.fastq > AllBTS.un1.fastq
cat 01A_0801_CTS.un2.fastq 01B_0802_CTS.un2.fastq 01C_0803_CTS.un2.fastq 01D_0804_CTS.un2.fastq 01E_0805_CTS.un2.fastq 01F_0806_CTS.un2.fastq 01G_0911_DSN_CTS.un2.fastq 02B_0906_CTS.un2.fastq 02C_0907_CTS.un2.fastq 02D_0908_CTS.un2.fastq 02E_0909_CTS.un2.fastq 02F_0910_CTS.un2.fastq 02G_0911_CTS.un2.fastq 03A_0905_CTS.un2.fastq 03D_0812_CTS.un2.fastq 03E_0901_CTS.un2.fastq 03F_0902_CTS.un2.fastq 03G_0903_CTS.un2.fastq 03H_0904_CTS.un2.fastq > AllCTS.un2.fastq
cat 01A_0905_BTS.un2.fastq 01B_0906_BTS.un2.fastq 01C_0907_BTS.un2.fastq 01D_0908_BTS.un2.fastq 01E_0901_BTS.un2.fastq 01E_0909_BTS.un2.fastq 01F_0902_BTS.un2.fastq 01F_0910_DSN_BTS.un2.fastq 01G_0903_BTS.un2.fastq 01H_0904_BTS.un2.fastq 02A_0905_BTS.un2.fastq 02F_0806_BTS.un2.fastq 02G_0807_BTS.un2.fastq 02H_0808_BTS.un2.fastq 03A_0809_BTS.un2.fastq 03B_0810_BTS.un2.fastq 03C_0811_BTS.un2.fastq 03H_0912_BTS.un2.fastq > AllBTS.un2.fastq
cat 01A_0801_CTS_joined_and_both_singles.fastq 01B_0802_CTS_joined_and_both_singles.fastq 01C_0803_CTS_joined_and_both_singles.fastq 01D_0804_CTS_joined_and_both_singles.fastq 01E_0805_CTS_joined_and_both_singles.fastq 01F_0806_CTS_joined_and_both_singles.fastq 01G_0911_DSN_CTS_joined_and_both_singles.fastq 02B_0906_CTS_joined_and_both_singles.fastq 02C_0907_CTS_joined_and_both_singles.fastq 02D_0908_CTS_joined_and_both_singles.fastq 02E_0909_CTS_joined_and_both_singles.fastq 02F_0910_CTS_joined_and_both_singles.fastq 02G_0911_CTS_joined_and_both_singles.fastq 03A_0905_CTS_joined_and_both_singles.fastq 03D_0812_CTS_joined_and_both_singles.fastq 03E_0901_CTS_joined_and_both_singles.fastq 03F_0902_CTS_joined_and_both_singles.fastq 03G_0903_CTS_joined_and_both_singles.fastq 03H_0904_CTS_joined_and_both_singles.fastq > AllCTS_joined_and_both_singles.fastq
cat 01A_0905_BTS_joined_and_both_singles.fastq 01B_0906_BTS_joined_and_both_singles.fastq 01C_0907_BTS_joined_and_both_singles.fastq 01D_0908_BTS_joined_and_both_singles.fastq 01E_0901_BTS_joined_and_both_singles.fastq 01E_0909_BTS_joined_and_both_singles.fastq 01F_0902_BTS_joined_and_both_singles.fastq 01F_0910_DSN_BTS_joined_and_both_singles.fastq 01G_0903_BTS_joined_and_both_singles.fastq 01H_0904_BTS_joined_and_both_singles.fastq 02A_0905_BTS_joined_and_both_singles.fastq 02F_0806_BTS_joined_and_both_singles.fastq 02G_0807_BTS_joined_and_both_singles.fastq 02H_0808_BTS_joined_and_both_singles.fastq 03A_0809_BTS_joined_and_both_singles.fastq 03B_0810_BTS_joined_and_both_singles.fastq 03C_0811_BTS_joined_and_both_singles.fastq 03H_0912_BTS_joined_and_both_singles.fastq > AllBTS_joined_and_both_singles.fastq
```

Step 3:
-------
Run the ARC assembly pipeline for 5 iterations using spades as the assembler and bowtie2 as the mapper. Here is the ARC_config.txt file:
```
## Name=value pairs:		
## reference: contains reference sequences in fasta format		
## numcycles: maximum number of times to try remapping		
## mapper: the mapper to use (blat/bowtie2)		
## assembler: the assembler to use (newbler/spades)		
## nprocs: number of cores to use		
## format: fasta or fasta, all must be the same		
## verbose: control mapping/assembly log generation (True/False)		
## urt: For Newbler, enable use read tips mode (True/False)		
## map_against_reads: On iteration 1, skip assembly, map against mapped reads (True/False)		
## assemblytimeout: kill assemblies and discard targets if they take longer than N minutes		
##		
## Columns:		
## Sample_ID:Sample_ID		
## FileName: path for fasta/fasta file		
## FileType: PE1, PE2, or SE		
## FileFormat: fasta or fasta		
# reference=targets.fasta		
# numcycles=5		
# mapper=bowtie2		
# assembler=spades		
# nprocs=31		
# format=fastq		
# verbose=True		
# urt=True		
# map_against_reads=False		
# assemblytimeout=300		
# bowtie2_k=5		
# rip=True		
# cdna=False		
# subsample=1		
# maskrepeats=True
Sample_ID	FileName	FileType
CTSandBTS5iter	CTS_and_BTS.un1.fastq	PE1
CTSandBTS5iter	CTS_and_BTS.un2.fastq	PE2
CTSandBTS5iter	CTS_and_BTS_joined_and_both_singles.fastq	SE
```

Step 4:
-------
Harvest the reciprocal best blast hits (RBBHs) from the ARC assembled contigs (contigs.fasta):<br>
`perl belowScript.pl --assembly contigs.fasta --targets ../targets.fasta --out RBBHs.fasta`

```perl
#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

my $help = 0;
my $assembly;
my $targets;
my $outFile;

GetOptions ("assembly=s"    => \$assembly,
            "targets=s"     => \$targets,
            "out=s"         => \$outFile,
            "help|man"      => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

#GetOptions  ("assemblybl2targ=s"      => \$assemblybl2targ,
#             "targbl2assembly=s"      => \$targbl2assembly,
#             "assemblyseqs=s"         => \$assemblySeqs,
#             "out=s"                  => \$outFile,
#             "help|man"               => \$help)           || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$targets or !$assembly or !$outFile or $help) {
    die "Must supply --assembly and --targets and --out.\nBlast commands should look like\
    blastn -db assembly -query targets.fasta -outfmt 6 -max_target_seqs 1 -out targets_bl2_assembly.blast and \
    blastn -db targets -query assembly.fasta -outfmt 6 -max_target_seqs 1 -out assembly_bl2_targets.blast";
}

system("makeblastdb -in $assembly -dbtype nucl -out findRBBHassembly");
system("makeblastdb -in $targets -dbtype nucl -out findRBBHtargets");
system("blastn -db findRBBHassembly -query $targets -outfmt 6 -max_target_seqs 1 -out findRBBH_targets_bl2_assembly.blast");
system("blastn -db findRBBHtargets -query $assembly -outfmt 6 -max_target_seqs 1 -out findRBBH_assembly_bl2_targets.blast");


# Make an index of all the sequences in the assembly
my %assemblyHash;
my $assemblyIn = Bio::SeqIO->new(-file => $assembly,
                                -format => 'fasta');
while (my $seq = $assemblyIn->next_seq()) {
    $assemblyHash{$seq->display_id()} = $seq;
}


# Find the reciprocal best blast hits
my $assembly2targetsBlast = Bio::SearchIO->new(-file => "findRBBH_assembly_bl2_targets.blast",
                                            -format => 'blasttable');

my %assem2targBlastResults;
while (my $result = $assembly2targetsBlast->next_result()) {
    $assem2targBlastResults{$result->query_name()} = ($result->next_hit())->name(); # This finds the name of the first hit and makes that the value for the key that is the name of the query
}

my $targ2assembBlast = Bio::SearchIO->new(-file => "findRBBH_targets_bl2_assembly.blast",
                                        -format => 'blasttable');
my %targ2assemBlastResults;
while (my $result = $targ2assembBlast->next_result()) {
    # print $result->query_name() . "\n";
    $targ2assemBlastResults{$result->query_name()} = ($result->next_hit())->name();  # This finds the name of the first hit and makes that the value for the key that is the name of the query
}

my $seqOut = Bio::SeqIO->new(-file=>">$outFile",
                            -format => 'fasta');

my $counter = 0;
foreach my $seqName (sort keys %targ2assemBlastResults) {
    # $seqName should be something like "Contig63"
    # $targ2assemBlastResults{$seqName} will be something like "6" (just a single number)
    # We want to see if $assem2targBlastResults{[the name of the contig in the assembly that was best match for target]} is the same as the name of the contig in the targets that matched the assembly
    if ($assem2targBlastResults{$targ2assemBlastResults{$seqName}} eq $seqName) {
        $counter++;

        my $seqDisplayName = $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id();
        my $newSeqName = $seqDisplayName . "_$seqName";
        $assemblyHash{$targ2assemBlastResults{$seqName}}->display_id($newSeqName);
        $seqOut->write_seq($assemblyHash{$targ2assemBlastResults{$seqName}});
    }
}
print "$counter total RBB hits found for assembly\n";

unlink("findRBBH_targets_bl2_assembly.blast", "findRBBH_targets_bl2_assembly.blast");
```

Step 5:
-------
Mask chimeric sequence and repetitive sequence as best we can. We do this by performing a self blast of the reciprocal best blast hits, then masking (replacing with Ns) all the HSP regions of the RBBHs that match another target:<br>
`makeblastdb -in RBBHs.fasta -dbtype nucl`<br>
`blastn -db RBBHs.fasta -query RBBHs.fasta -outfmt 6 -evalue 1e-20 -out RBBHs_bl2_self.blast`<br>
`perl belowScript.pl --in RBBHs_bl2_self.blast --sequences RBBHs.fasta --out RBBHs5iter_chimeramasked.fasta`

```perl
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;

my $help = 0;
my $inFile;
my $sequences;
my $outFile;

GetOptions  ("in=s"      => \$inFile,
             "sequences=s" => \$sequences,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outFile or $help) {
    die "Must supply --in and --out.\n";
}


#### Example output that we're parsing: ####
# All_:_contig140316|NEMF|11_:_Unfinished001_contig140316|NEMF|11 All_:_contig140316|NEMF|11_:_Unfinished001_contig140316|NEMF|11 100.00  1451    0       0       1       1451    1       1451    0.0     2680
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   100.00  2045    0       0       1       2045    1       2045    0.0     3777
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig162323|MLLT4|8_:_Unfinished001_contig162323|MLLT4|8 84.21   190     7       9       1       167     111     300     1e-39    163
# All_:_contig140356|USP3|3_:_Unfinished001_contig140356|USP3|3   All_:_contig71398|6-Mar|1_:_Unfinished005_contig71398|6-Mar|1   83.78   185     7       5       1       162     33      217     7e-37    154

my $searchIO = Bio::SearchIO->new(-file=>$inFile,
                                  -format=>'blasttable');
my $seqIn = Bio::SeqIO->new(-file=>$sequences,
                            -format=>'fasta');

my %seqHash;
while (my $seq = $seqIn->next_seq()) {
    
    $seqHash{$seq->display_id()} = $seq->seq();
}
my %maskingHash;
while (my $result = $searchIO->next_result()) {
    my $hitCounter = 0;
    while (my $hit = $result->next_hit()) {
        # Look for multiple hits. First one should be 100% and be to itself
        # For all other hits, we want to mask the regions that are blasting to those regions
        if ($hitCounter == 0) {
            $hitCounter++;
            next;
        }
        while (my $hsp = $hit->next_hsp()) {
            # Translate the query and subject names into the actual target names
            my $subjectName = $hit->name();
            my $queryName = $result->query_name();

            # Now find the offending sections:
            my $queryStart = $hsp->start('query');
            my $queryEnd = $hsp->end('query');
            my $queryString = $queryStart . ":" . $queryEnd;
            my $subjectStart = $hsp->start('subject');
            my $subjectEnd = $hsp->end('subject');
            my $subjectString = $subjectStart . ":" . $subjectEnd;
            # And store them:
            push(@{$maskingHash{$queryName}}, $queryString);
            push(@{$maskingHash{$subjectName}}, $subjectString);
        }
        # For all of these remaining hits, we want to eventually mask the overlapping sequence in both the query and the subject contigs
        $hitCounter++;   
    }
}


foreach my $target (sort keys %maskingHash) {
    # Let's individually mask every one of these regions...
    
    foreach my $badRegion (@{$maskingHash{$target}}) {
        my @fields = split(/:/, $badRegion);
        unless ($fields[0] < $fields[1]) {
            die '$fields[0] is less than $fields[1] for target ' . $target . ". We are assuming they're sorted...\n";
        }
        #print "Sequence: " . $seqHash{$target};
        #print "Badbeginning: $fields[0]\n";
        #print "Badend: fields[1]\n";
        my $badLength = $fields[1] - $fields[0] + 1; # The +1 is to include both ends
        #print "Badlength = $badLength\n";
        
        my $replacement = "N" x $badLength; # This will give, for instance, NNNN if $badLength = 4
        substr($seqHash{$target}, $fields[0], $badLength, $replacement);
        #print $seqHash{$target} . "\n";
    }
}


# print Dumper(\%maskingHash);
# print Dumper(\%seqHash);


open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";

foreach my $target (sort keys %seqHash) {
    print $outFH ">$target\n$seqHash{$target}\n";
}
```

Step 6:
-------
Note the portions of the chimera-masked regions that correspond to the actual target sequences:<br>
`makeblastdb -in RBBHs5iter_chimeramasked.fasta -dbtype nucl`<br>
`blastn -db RBBHs5iter_chimeramasked.fasta -query ../targets.fasta -outfmt 6 -out targets_bl2_RBBHs5iter_chimeramasked_e1e-20.blast -evalue 1e-20`

Step 7:
-------
Map the reads for each library to the chimera-masked RBBHs, run Picard MarkDuplicates, and summarize using samtools flagstat:<br>

`cd /dir/with/libraryFastqs/`<br>
`perl belowScript.pl`

```perl
#!/usr/bin/perl

use strict;
use warnings;

my @samples = ("01A_0801_CTS", "01A_0809_F1", "01A_0905_BTS", "01B_0802_CTS", "01B_0810_F1", "01B_0906_BTS", "01C_0803_CTS", "01C_0811_F1", "01C_0907_BTS", "01D_0804_CTS", "01D_0812_F1", "01D_0908_BTS", "01E_0805_CTS", "01E_0901_BTS", "01E_0909_BTS", "01F_0806_CTS", "01F_0902_BTS", "01F_0910_DSN_BTS", "01G_0807_F1", "01G_0903_BTS", "01G_0911_DSN_CTS", "01H_0808_F1", "01H_0904_BTS", "02A_0801_F1", "02A_0905_BTS", "02B_0802_F1", "02B_0906_CTS", "02C_0803_F1", "02C_0907_CTS", "02D_0804_F1", "02D_0908_CTS", "02E_0805_F1", "02E_0909_CTS", "02F_0806_BTS", "02F_0910_CTS", "02G_0807_BTS", "02G_0911_CTS", "02H_0808_BTS", "02H_0912_F1", "03A_0809_BTS", "03A_0905_CTS", "03B_0810_BTS", "03B_0906_F1", "03C_0811_BTS", "03C_0907_F1", "03D_0812_CTS", "03D_0908_F1", "03E_0901_CTS", "03E_0909_F1", "03F_0902_CTS", "03F_0910_F1", "03G_0903_CTS", "03G_0911_F1", "03H_0904_CTS", "03H_0912_BTS");

foreach my $sample (@samples) {
    my $command1 = "bwa mem -t 30 -M /home/evan/manderReads/finished_CTSandBTS5iter/RBBHs5iter_chimeramasked.fasta " . $sample . ".un1.fastq " . $sample . ".un2.fastq | samtools view -@ 30 -bS - | samtools sort -@ 30 -T temp -o " . $sample . ".paired.bam -";
    my $command2 = "bwa mem -t 30 -M /home/evan/manderReads/finished_CTSandBTS5iter/RBBHs5iter_chimeramasked.fasta " . $sample . "_joined_and_both_singles.fastq | samtools view -@ 30 -bS - | samtools sort -@ 30 -T temp -o " . $sample . ".singles.bam -";
    my $command3 = "samtools merge -@ 30 " . $sample . ".merged.bam " . $sample . ".paired.bam " . $sample . ".singles.bam";
    my $command4 = "java -jar ~/bin/picard-tools-1.125/picard.jar MarkDuplicates I=" . $sample . ".merged.bam O=" . $sample . ".mergedMD.bam METRICS_FILE=" . $sample . ".mergedMD.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false";
    my $command5 = "samtools flagstat " . $sample . ".mergedMD.bam";
    system($command1);
    system($command2);
    system($command3);
    system($command4);
    system($command5);
}
```

Step 8:
-------
Create the depth files.

```perl
#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

my @samples = ("01A_0801_CTS", "01A_0809_F1", "01A_0905_BTS", "01B_0802_CTS", "01B_0810_F1", "01B_0906_BTS", "01C_0803_CTS", "01C_0811_F1", "01C_0907_BTS", "01D_0804_CTS", "01D_0812_F1", "01D_0908_BTS", "01E_0805_CTS", "01E_0901_BTS", "01E_0909_BTS", "01F_0806_CTS", "01F_0902_BTS", "01F_0910_DSN_BTS", "01G_0807_F1", "01G_0903_BTS", "01G_0911_DSN_CTS", "01H_0808_F1", "01H_0904_BTS", "02A_0801_F1", "02A_0905_BTS", "02B_0802_F1", "02B_0906_CTS", "02C_0803_F1", "02C_0907_CTS", "02D_0804_F1", "02D_0908_CTS", "02E_0805_F1", "02E_0909_CTS", "02F_0806_BTS", "02F_0910_CTS", "02G_0807_BTS", "02G_0911_CTS", "02H_0808_BTS", "02H_0912_F1", "03A_0809_BTS", "03A_0905_CTS", "03B_0810_BTS", "03B_0906_F1", "03C_0811_BTS", "03C_0907_F1", "03D_0812_CTS", "03D_0908_F1", "03E_0901_CTS", "03E_0909_F1", "03F_0902_CTS", "03F_0910_F1", "03G_0903_CTS", "03G_0911_F1", "03H_0904_CTS", "03H_0912_BTS");

my $forkManager = new Parallel::ForkManager(5);
foreach my $sample (@samples) {
    $forkManager->start and next;
    my $command = "samtools depth -q 20 -Q 20 " . $sample . ".mergedMD.bam > " . $sample . ".depth"; 
    system($command);
    $forkManager->finish;
}
$forkManager->wait_all_children;
```

Step 9:
-------
Harvest the depths across all samples for each locus:

`mkdir CTSandBTS5iterLocusDepthFiles`<br>
`perl belowScript.pl --out CTSandBTS5iterLocusDepthFiles`
```perl
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $help = 0;
my $outDir;

GetOptions  ("out=s"      => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$outDir or $help) {
    die "Must supply --out <outputDirectory>.\n";
}

unless (-d $outDir) {
    mkdir $outDir;
}

my @samples = ( "01A_0801_CTS","01A_0809_F1","01A_0905_BTS","01B_0802_CTS","01B_0810_F1","01B_0906_BTS","01C_0803_CTS","01C_0811_F1","01C_0907_BTS","01D_0804_CTS","01D_0812_F1","01D_0908_BTS","01E_0805_CTS","01E_0901_BTS","01E_0909_BTS","01F_0806_CTS","01F_0902_BTS","01F_0910_DSN_BTS","01G_0807_F1","01G_0903_BTS","01G_0911_DSN_CTS","01H_0808_F1","01H_0904_BTS","02A_0801_F1","02A_0905_BTS","02B_0802_F1","02B_0906_CTS","02C_0803_F1","02C_0907_CTS","02D_0804_F1","02D_0908_CTS","02E_0805_F1","02E_0909_CTS","02F_0806_BTS","02F_0910_CTS","02G_0807_BTS","02G_0911_CTS","02H_0808_BTS","02H_0912_F1","03A_0809_BTS","03A_0905_CTS","03B_0810_BTS","03B_0906_F1","03C_0811_BTS","03C_0907_F1","03D_0812_CTS","03D_0908_F1","03E_0901_CTS","03E_0909_F1","03F_0902_CTS","03F_0910_F1","03G_0903_CTS","03G_0911_F1","03H_0904_CTS","03H_0912_BTS");

foreach my $sample (@samples) {
    my $depthFile = $sample . ".depth";
    open(my $depthFH, "<", $depthFile) or die "Couldn't open $depthFile for reading: $!\n";
    while (my $line = <$depthFH>) {
        my @fields = split(/\t/, $line);
        if ($fields[0] =~ /.*\_\:\_(.*)\_\:\_/) {
            my $locus = $1;
            $locus =~ s/\|/\-\-/g; # Replace all the |'s with --'s
            my $fileName = $outDir . "/" . $locus . ".depth";
            open(my $outFH, ">>", $fileName) or die "Couldn't open $fileName for appending: $!\n";
            print $outFH "$sample\t$line";
        }
    }
}
```

Step 10:
--------
Plot the 

First, harvest the filenames of all the locus depth files: <br>
`cd CTSandBTS5iterLocusDepthFiles;`<br>
`ls -al | grep "\.depth" | awk '{print $(NF)}' > depthFiles.txt`<br>

Then plot all of them into their own png files, drawing a black bar where the actual
target sequence overlaps.


```R
library("data.table") # Might have to install this
depthFiles <- fread("depthFiles.txt", sep="\n", header=FALSE)
rainbowColors <- rainbow(55) # 55 because we have 55 potential libraries to plot for each locus 
blastResults <- fread("targets_bl2_RBBHschimeramasked.blast")

# Now depthFiles$V1 holds all of the individual loci files--we'll read them one-by-one and create the plots
for (i in 1:length(depthFiles$V1)) {
  locusData <- fread(depthFiles$V1[i])
  samples <- unique(locusData$V1) # Not all samples are in all targets. Get # for each target
  imageName <- paste("plots/", depthFiles$V1[i], ".png", sep="")
  
  #depthFiles$V1[i] will be something like gi62424-Wnt1-OPA.depth. We need to translate that to gi62424|Wnt1|OPA, which we can then look up in the blastResults table to show us where to draw our black bar showing the target region
  locusName <- gsub("@", "|", depthFiles$V1[i]) # Change ampersands in filenames to bars
  locusName <- gsub(".depth", "", locusName) # Remove the file extension
  targetStart <- blastResults$V9[blastResults$V1 == locusName] # Harvest the subject HSP start site from blast output
  targetEnd <-blastResults$V10[blastResults$V1 == locusName] # Harvest the subject HSP end site from blast output
  write(paste("Processing ", depthFiles$V1[i], ". locus: ", locusName, ". targetStart: ", targetStart, ". targetEnd: ", targetEnd, "\n"))
  
  png(filename = imageName, width=1000, height=600, units="px", pointsize=12)
  # Make the initial plot
  plot(locusData$V3[locusData$V1 == samples[1]], locusData$V4[locusData$V1 == samples[1]], type="l", xlim=c(0,max(locusData$V3)), ylim=c(0,50), col=rainbowColors[1], lwd=1.2, main=depthFiles$V1[i], xlab="Position in contig", ylab="Read depth")
  
    # Add all the lines for the rest of the samples.
  # Start at #2 because we just did #1 above
  for (i in 2:length(samples)) {
    lines(locusData$V3[locusData$V1 == samples[i]], locusData$V4[locusData$V1 == samples[i]], type="l", col=rainbowColors[i], lwd=1.2)
  }
  
  # Plot black line where the target blasted to the assembled contig that represents that target.
  # This actually works if there are multiple HSPs for that target, because targetStart and
  # targetEnd become arrays, and the segments call parses them properly.
  segments(targetStart, 0, targetEnd, 0, lwd=5)
  dev.off()
  rm(locusData)
  gc() # The depth files are potentially large--do as much as we can to reduce memory usage
}
```

Finishing up:
-------------
Here is what the results look like for a few loci:
