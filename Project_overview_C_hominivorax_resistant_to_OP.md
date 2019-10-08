# 📋 Project overview: C hominivorax resistant to OP


| **Status**      | In progress 🔨  |
| --------------- | --------------- |
| **Team**        | Sophie, Tatiana |
| **Description** | RNA-seq         |
| **Deadline**    | 30/03/2017      |

# Timeline
- What are the key dates?


## Action items
[x] @someone let’s get this done
[ ]  
# Background
## Goal

Where’s this project headed?


## Research
- Add related docs by typing a plus sign (+)


# Project

10.03.2017

## 1 Files

Sequenced at VetmedUni Wien (1st replicate)
a. /Volumes/HD2/seqs/Chom/RNAseq-resistance/CS/

1st_run/

- Chom_illumina01-ACGTT-control.fastq     10.332.646
- Chom_illumina01-TGCAT-resistant.fastq     9.191.047

2nd_run/

- Chom_Illumina02-ACGTT.fastq     2.128
- Chom_Illumina02-TGCAT.fastq     8.533

* Thousands of sequences; quality not guaranteed. They were not used.

3rd_run/
* mixed with Dmelanogaster samples and were sorted by using sort_reads_fastq.pl

ACGTT = control
TGCAT = resistant


    perl /Volumes/HD2/scripts/Illumina/sort_reads_fastq.pl -s /Volumes/HD2/seqs/Chom/RNAseq-resistance/CS/3rd\ run/s_2+3_sequence_read1+2_101011/s_2_1_sequence.txt -b ACGTT -n Chom-control
    
    perl /Volumes/HD2/scripts/Illumina/sort_reads_fastq.pl -s /Volumes/HD2/seqs/Chom/RNAseq-resistance/CS/3rd\ run/s_2+3_sequence_read1+2_101011/s_3_1_sequence.txt -b TGCAT -n Chom-resistant


- Chom-control-ACGTT.fastq        2.696.205
- Chom-resistant-TGCAT.fastq     4.134.932

Sequenced at ESALQ (2nd replicate)
b. /Volumes/HD2/seqs/Calliphoridae-Oestridae/1st-run/Sample_Chom-C/ and /Volumes/HD2/seqs/Calliphoridae-Oestridae/1st-run/Sample_Chom-R/

13.03.2017

1.1 Concatenating files
a. First replicate (/Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/)

- Control-1st-R1.fastq = 13.028.851 sequences
- Resistant-1st-R1.fastq = 13.325.979 sequences

b. Second replicate - paired-end (/Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/)

- Control-2nd-R1.fastq = 26.807.116
- Control-2nd-R2.fastq = 26.807.116
- Resistant-2nd-R1.fastq = 24.558.107
- Resistant-2nd-R2.fastq = 24.558.107


## 2 Quality trimming 

Trimmomatic: A flexible read trimming tool for Illumina NGS data
Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Version 0.36

2.1. First replicate

- Control-1st-R1.fastq = 13.028.851 sequences
    java -jar /Applications/seq_an/Trimmomatic-0.36/trimmomatic-0.36.jar SE /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/Control-1st-R1.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/Control-1st-R1-trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Quality encoding detected as phred64
Input Reads: 13028851 Surviving: 11055926 (84.86%) Dropped: 1972925 (15.14%)


- Resistant-1st-R1.fastq = 13.325.979 sequences
    java -jar /Applications/seq_an/Trimmomatic-0.36/trimmomatic-0.36.jar SE /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/Resistant-1st-R1.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/Resistant-1st-R1-trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

Quality encoding detected as phred64
Input Reads: 13325979 Surviving: 11279633 (84.64%) Dropped: 2046346 (15.36%)

2.2. Second replicate

- Control-2nd-R1.fastq/Control-2nd-R2.fastq = 26.807.116 sequences
    java -jar /Applications/seq_an/Trimmomatic-0.36/trimmomatic-0.36.jar PE /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R1.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R2.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R1-trimmed-paired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R1-trimmed-unpaired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R2-trimmed-paired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Control/Control-2nd-R2-trimmed-unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Quality encoding detected as phred33
Input Read Pairs: 26807116 Both Surviving: 26277392 (98.02%) Forward Only Surviving: 361984 (1.35%) Reverse Only Surviving: 103989 (0.39%) Dropped: 63751 (0.24%)


- Resistant-2nd-R1.fastq/Resistant-2nd-R2.fastq = 24.558.107 sequences
    java -jar /Applications/seq_an/Trimmomatic-0.36/trimmomatic-0.36.jar PE /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R1.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R2.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R1-trimmed-unpaired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R2-trimmed-paired.fastq /Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/2nd-replicate/Resistant-2nd-R2-trimmed-unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

Quality encoding detected as phred33
Input Read Pairs: 24558107 Both Surviving: 24091747 (98.10%) Forward Only Surviving: 314072 (1.28%) Reverse Only Surviving: 94667 (0.39%) Dropped: 57621 (0.23%)

## 3 Collapsing sequences 

Two scripts; for single and paired-end sequences:

a. collapse_reads_fastq.pl (single-end)

    #                                                        
    #  PROGRAM: collapse_reads_fastq.pl                          07.01.2013                   
    #
    #  DESCRIPTION: Script to collapse identical sequences and report unique
    #               reads with the highest average quality in fastq format.  
    #
    #  (header modification)
    #
    #  USAGE: collapse_reads_fastq.pl -s Input fastq sequence file 
    #         -e Quality scores encoding (33 or 64) -n Project name
    #
    #  AUTHORS: Tatiana Torres (tttorres@ib.usp.br)
    #           Sophie Tandonnet
    #
    #  CITATION: Tandonnet, S. and Torres, T.T. Unpublished. 3' RNA-seq
    #            method for non-model species.
    #
    #  LAST MODIFIED: 22.03.2014 
    #
    
    
    
    #!usr/bin/perl -w
    
    use Getopt::Std;
    use warnings;
    use strict;
    
    
    
    ### Declare and initialize variables 
    
    my $input_seq = my $project_name = my $encoding = "";
    
    my $seq_id = my $seq = my $qul_id = my $qul = ""; 
    my $next_seq = my $next_qul = 0;
    
    my %opts = my %summary = ();
    
    my %counts = my % seqID = my %qualID = my %qual = ();
    
    ### Read command line
    
    my ($USAGE) = "\nUSAGE: $0\n".
                          "\t\t-s Input sequence file - sequences in FASTQ format\n".
                          "\t\t-e Quality scores encoding - 33 or 64 [default = 64]\n".            
                  "\t\t-n Project name [default = collapsed]\n\n";
    
    getopts('s:e:n:h', \%opts);
    chomp(%opts);
    
    if ($opts{s}) {
       $input_seq = $opts{s};
       chomp $input_seq;
    }
    
    if ($opts{e}) {
        $encoding = $opts{e};
        chomp $encoding;
    } else { $encoding = 64 }
    
    if ($opts{n}) {
       $project_name = $opts{n};
       chomp $project_name;
    } else { $project_name = "collapsed" }
    
    
    if (!($input_seq) || ($opts{h})) {
            print $USAGE;
            exit;
    }
    
    ###Open input file
    open(SEQ, "<$input_seq")  || die "Could not open file $input_seq.\n";
    
    
    ### Create new directory to store new files
    
    mkdir $project_name || die "Could not create folder $project_name\n";
    
    ### Output files
    
    my $out = $project_name."/". $project_name . ".fastq"; #collapsed sequences
    open(OUT, ">$out") || die "Could not open file $out.\n";
    
    my $table = $project_name."/". $project_name . "stats.txt"; #duplicated reads counts
    open(TAB, ">$table") || die "Could not open file $out.\n";
            
    
    ### Sorting sequences 
    while(<SEQ>) {
            
            my $seqID = $_;        # @HWI-ST365_0150:7:1101:1441:2340#CAAGTA/1
        my $seq = <SEQ>;       # TATTAACTGTTTCATAAAAATTTGTGATAGAATTTTATTTTTTAGTC
        my $qualID = <SEQ>;    # +HWI-ST365_0150:7:1101:1441:2340#CAAGTA/1
        my $qual = <SEQ>;      # daddaeefffffddfc`c\Zdddd`fdfcdfdfffffffffffcfee
        
        if ($seqID{$seq}){
            ++$counts{$seq};
            
            my $totalQual = 0;
            
            while ($qual =~ /./g) {
                            my $qbase = ord($&) - $encoding;
                $totalQual += $qbase;
                    }
                    
            my $recQual = $qual{$seq};  # quality of recorded sequence
            
            my $recTotalQual = 0;
            
            while ($recQual =~ /./g) {
                            my $qbase = ord($&) - $encoding;
                $recTotalQual += $qbase;
                    }
                    
            if ($recTotalQual < $totalQual){
                $seqID{$seq}  = $seqID;
                $qualID{$seq} = $qualID;
                $qual{$seq}   = $qual;
            }
            
        } else {
            $counts{$seq} = 1;        # number of occurrences of the sequence
            $seqID{$seq}  = $seqID;   # sequence ID
            $qualID{$seq} = $qualID;  # quality ID
            $qual{$seq}   = $qual;    # quality
        }
        
    }
        
    foreach my $seq(keys(%seqID)){
        print OUT $seqID{$seq}. $seq . $qualID{$seq}. $qual{$seq};
    }
    
    foreach $seq (sort { $counts{$b} <=> $counts{$a} } (keys(%counts))) {
        my $counts = $counts{$seq};
        chomp($seq);
        print TAB "$seq\t$counts\n";
    }
    
    exit;

b. collapse_reads_fastq_PE.pl (paired-end reads)

    #                                                        
    #  PROGRAM: collapse_reads_fastq_PE.pl                          27.05.2013                   
    #
    #  DESCRIPTION: Script to collapse identical sequences and report unique
    #               reads with the highest average quality in fastq format.  
    #
    #  USAGE: collapse_reads_fastq.pl -s Input fastq sequence file 
    #         -e Quality scores encoding (33 or 64) -n Project name
    #
    #  AUTHORS: Tatiana Torres (tttorres@ib.usp.br)
    #           Sophie Tandonnet
    #
    #  CITATION: Tandonnet, S. and Torres, T.T. Unpublished. 3' RNA-seq
    #            method for non-model species.
    #
    #  LAST MODIFIED: 22.03.2014 
    #
    
    
    
    #!usr/bin/perl -w
    
    use Getopt::Std;
    use warnings;
    use strict;
    
    
    
    ### Declare and initialize variables 
    
    my $inputLeft = my $inputRight = "";
    
    my $tableLeft = my $tableRright = "";
    
    my $projectName = my $encoding = "";
    
    my $seq_id = my $seq = my $qul_id = my $qul = ""; 
    my $next_seq = my $next_qul = 0;
    
    my %opts = my %summary = ();
    
    my %counts = my % seqID = my %qualID = my %qual = ();
    
    my $counter = 0;
    
    
    ### Read command line
    
    my ($USAGE) = "\nUSAGE: $0\n".
                          "\t\t-l Input sequence file - left reads in FASTQ format\n".
                          "\t\t-r Input sequence file - right reads in FASTQ format\n".
                          "\t\t-e Quality scores encoding - 33 or 64 [default = 33]\n".            
                  "\t\t-n Project name [default = collapsed]\n\n";
    
    getopts('l:r:e:n:h', \%opts);
    chomp(%opts);
    
    if ($opts{l}) {
       $inputLeft = $opts{l};
       chomp $inputLeft;
    }
    
    if ($opts{r}) {
        $inputRight = $opts{r};
        chomp $inputRight;
    }
    
    if ($opts{e}) {
        $encoding = $opts{e};
        chomp $encoding;
    } else { $encoding = 33 }
    
    if ($opts{n}) {
       $projectName = $opts{n};
       chomp $projectName;
    } else { $projectName = "collapsed" }
    
    
    if (!($inputLeft && $inputRight ) || ($opts{h})) {
            print $USAGE;
            exit;
    }
    
    ###Open input file
    open(LEFT,  "<$inputLeft")  || die "Could not open file $inputLeft.\n";
    open(RIGHT, "<$inputRight") || die "Could not open file $inputRight.\n";
    
    ### Create new directory to store new files
    
    mkdir $projectName || die "Could not create folder $projectName\n";
    
    ### Output files
    
    my $outLeft  = $projectName."/". $projectName . "_R1.fastq"; #collapsed sequences
    open(OUTL, ">$outLeft")  || die "Could not open file $outLeft.\n";
    
    my $outRight = $projectName."/". $projectName . "_R2.fastq"; #collapsed sequences
    open(OUTR, ">$outRight") || die "Could not open file $outRight.\n";
    
    my $table = $projectName."/". $projectName . "_stats.txt"; #duplicated reads counts
    open(TAB, ">$table") || die "Could not open file $table.\n";
    
            
    
    ### Sorting sequences 
    while(<LEFT>) {
        
        ++$counter;
            
        # left read
            my $seqID_L  = $_;        # @HWI-1KL182:77:C26REACXX:4:1101:7500:2222 1:N:0:CAGATC
        my $seq_L    = <LEFT>;    # TATTAACTGTTTCATAAAAATTTGTGATAGAATTTTATTTTTTAGTC
        my $qualID_L = <LEFT>;    # +HWI-ST365_0150:7:1101:1441:2340#CAAGTA/1
        my $qual_L   = <LEFT>;    # daddaeefffffddfc`c\Zdddd`fdfcdfdfffffffffffcfee
        
        # right read
            my $seqID_R  = <RIGHT>;   # @HWI-1KL182:77:C26REACXX:4:1101:7500:2222 2:N:0:CAGATC
        my $seq_R    = <RIGHT>;   # TATTAACTGTTTCATAAAAATTTGTGATAGAATTTTATTTTTTAGTC
        my $qualID_R = <RIGHT>;   # +HWI-ST365_0150:7:1101:1441:2340#CAAGTA/2
        my $qual_R   = <RIGHT>;   # daddaeefffffddfc`c\Zdddd`fdfcdfdfffffffffffcfee
        
        # concatenated sequence: left+right
            my $seqID  = $seqID_L  ."~". $seqID_R;
        my $seq    = $seq_L    ."~". $seq_R;
        my $qualID = $qualID_L ."~". $qualID_R;
        my $qual   = $qual_L   ."~". $qual_R;
        
        if ($seqID{$seq}){
            
            ++$counts{$seq};
            
            my $totalQual = 0;
            
            while ($qual =~ /./g) {
                            my $qbase = ord($&) - $encoding;
                $totalQual += $qbase;
                    }
                    
            my $recQual = $qual{$seq};  # quality of recorded sequence
            
            my $recTotalQual = 0;
            
            while ($recQual =~ /./g) {
                            my $qbase = ord($&) - $encoding;
                $recTotalQual += $qbase;
                    }
            
                    
            if ($recTotalQual < $totalQual){
                $seqID{$seq}  = $seqID;
                $qualID{$seq} = $qualID;
                $qual{$seq}   = $qual;
            }
            
        } else {
            
            $counts{$seq} = 1;        # number of occurrences of the sequence
            $seqID{$seq}  = $seqID;   # sequence ID (hash with observed sequences)
            $qualID{$seq} = $qualID;  # quality ID
            $qual{$seq}   = $qual;    # quality
        }
        
    }
        
    foreach my $seq(keys(%seqID)){
        
        (my $seqID_L, my $seqID_R)   = split("~", $seqID{$seq});
        (my $seq_L, my $seq_R)       = split("~", $seq);
        (my $qualID_L, my $qualID_R) = split("~", $qualID{$seq});
        (my $qual_L, my $qual_R)     = split("~", $qual{$seq});
        
        print OUTL $seqID_L . $seq_L . $qualID_L . $qual_L;
        print OUTR $seqID_R . $seq_R . $qualID_R . $qual_R;
    }
    
    foreach $seq (sort { $counts{$b} <=> $counts{$a} } (keys(%counts))) {
        my $counts = $counts{$seq};
        $seq =~ s/\n~/\./g;
        chomp($seq);
        print TAB "$seq\t$counts\n";
    }
    
    exit;


3.1. Single end files (1st replicate)

- Control-1st-R1-trimmed.fastq = 11.055.926 sequences
    perl /Users/Tatiana/scripts/Illumina/collapse_reads_fastq.pl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/1st-replicate/Control-1st-R1-trimmed.fastq -n Control-1st-R1-collapsed

Surviving: 4.806.021


- Resistant-1st-R1-trimmed.fastq = 11.279.633 sequences
    perl /Users/Tatiana/scripts/Illumina/collapse_reads_fastq.pl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/1st-replicate/Resistant-1st-R1-trimmed.fastq -n Resistant-1st-R1-collapsed

Surviving: 5.845.907

3.2. Paired-end files (2nd replicate)

- Control-2nd-R1-trimmed-paired.fastq/Control-2nd-R2-trimmed-paired.fastq = 26.277.392 sequences
    perl /Users/Tatiana/scripts/Illumina/collapse_reads_fastq_PE.pl -l /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/2nd-replicate/Control-2nd-R1-trimmed-paired.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/2nd-replicate/Control-2nd-R2-trimmed-paired.fastq -n Control-2nd-collapsed

Surviving: 20.453.006


- Resistant-2nd-R1-trimmed-paired.fastq/Resistant-2nd-R2-trimmed-paired.fastq = 24.091.747 sequences
    perl /Users/Tatiana/scripts/Illumina/collapse_reads_fastq_PE.pl -l /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/1-trimming/2nd-replicate/Resistant-2nd-R2-trimmed-paired.fastq -n Resistant-2nd-collapsed

Surviving: 19.130.494


14.03.2017

## 4 Assembly 

RNA-Seq De novo Assembly Using Trinity
Trinity version: trinityrnaseq_r20140717

All sequences
If you have RNA-Seq data from multiple libraries and you want to run them all through Trinity in a single pass, simply combine all your left.fq files into one left.fq file, and combine all right.fq files into one right.fq file. Then run Trinity using these separately concatenated left and right input files.

Files:

- Control-collapsed_R1.fastq (Control, 1st and 2nd replicate)
- Control-2nd-collapsed_R2.fastq (Control, 2nd replicate)
- Resistant-collapsed_R1.fastq (Resistant, 1st and 2nd replicate)
- Resistant-2nd-collapsed_R2.fastq (Resistant, 2nd replicate)

Running Trinity:

    Trinity --JM 20G --seqType fq --left Control-collapsed_R1.fastq,Resistant-collapsed_R1.fastq --right Control-2nd-collapsed_R2.fastq,Resistant-2nd-collapsed_R2.fastq --CPU 20 --inchworm_cpu 20 --bflyCPU 20 --normalize_reads --full_cleanup --out /mnt/HD4/Chom-resistance/Trinity

Start: Wednesday, March 15, 2017: 19:29:56
Finish: Wednesday, March 15, 2017: 20:28:12

Number of contigs: 40277

Average contig length: 902.187

    perl -ne 'print "$1\n" if /\>.*len\=(\d*)/g' headers | awk '{s+=$1; n++} END {print s/n}' 



14.03.2017

## 5 Assembly quality

BUSCO 
Assessing genome assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs (BUSCO)

Using Diptera database: http://busco.ezlab.org/datasets/diptera_odb9.tar.gz

    python /Applications/seq_an/busco/BUSCO.py -i /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta -o Chom-busco -m tran -l /Applications/seq_an/busco/diptera_odb9 -c 20 -e 0.0001 --long 

INFO        Results:
INFO        C:79.1%[S:60.5%,D:18.6%],F:8.2%,M:12.7%,n:2799
INFO        2215 Complete BUSCOs (C)
INFO        1693 Complete and single-copy BUSCOs (S)
INFO        522 Complete and duplicated BUSCOs (D)
INFO        230 Fragmented BUSCOs (F)
INFO        354 Missing BUSCOs (M)
INFO        2799 Total BUSCO groups searched


## 6 Read mapping

6.1. Discarding isoforms, retaining the longest
Perl script: 

    #!usr/bin/perl -w
    
    ##############################################################################
    #  PROGRAM: isoforms.pl                                          25.05.2015  #
    #                                                                            #
    #  DESCRIPTION: Picks the longest isoform from Trinity assembly.             #
    #                                                                            #
    #  AUTHOR: Raquel D. Monfardini and Gisele A. Cardoso                        #
    #                                                                            #
    #  LAST MODIFIED: 27.05.2015                                                 #
    ##############################################################################
    
    use warnings;
    use strict;
    use Getopt::Std;
    
    ### Declare and initialize variables
    
    my %opts;
    my $name = "x";
    my $length = "";
    my $head = "";
    my $i = 0;
    my $seq = "";
    my $input_file;
    my $output_file;
    
    ### Input and output arguments
    
    my ($USAGE) = "\nUSAGE: $0\n".
    "\t\t-i Input file [Trinity sequences in fasta]\n".
    "\t\t-o Output name [default = isoforms.txt]\n\n";
    
    getopts('i:o:h', \%opts);
    chomp(%opts);
    
    if ($opts{i}) {
        $input_file = $opts{i};
        chomp $input_file;
    }
    if ($opts{o}) {
        $output_file = $opts{o};
        chomp $output_file;
    } else {$output_file = "isoforms.txt";}
    
    if (!$opts{i} || ($opts{h})) {
        print $USAGE;
        exit;
    }
    
    ### Open input files
    open (IN, "<$input_file") || die "Could not open file $input_file\n";
    
    ### Open output files
    open (OUT, ">$output_file") || die "Could not open file $output_file\n";
    
    ##Getting the longest isoform
    while (<IN>) {
        if ($_ =~ />(.+)i\d+\ len=(\d+)/) {
            if ($1 eq $name && $2 > $length) {
                $head = $_;
                $name = $1;
                $length = $2;
                $seq = "";
                $i = 1;
            } elsif ($1 eq $name && $2 < $length) { $i = 0 }
            else {
                print OUT $head.$seq;
                $head = $_;
                $name = $1;
                $length = $2;
                $seq = "";
                $i = 1;
            }
        }
        else {
            if ($i == 1) {
                $seq = $seq.$_;
            }
        }    
    }
    
    print OUT $head.$seq."\n";
    
    exit;

Running script: 

    perl /Users/Tatiana/Dropbox/Scripts/isoforms.pl -i /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.longest.fasta

Contigs retained: 33.397

6.2. Reducing redundancy by clustering
CD-HIT
http://weizhongli-lab.org/cd-hit/

a. all contigs

    ./cd-hit-est -i /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.clusters.fasta -c 0.95 -M 0 -T 0

Contigs retained: 37.270

b. longest isoforms

    ./cd-hit-est -i /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.longest.fasta -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.longest.clusters.fasta -c 0.95 -M 0 -T 0 

Contigs retained: 33.038

**Trinity.longest.clusters.fasta was used as database for mapping**

6.3. Building reference for alignments
bowtie2-build 
Bowtie 2 version 2.2.3 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)

    bowtie2-build /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant

6.4. Alignments
bowtie2-build 
Bowtie 2 version 2.2.3 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea) 

a. First replicate (/Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/)

Options not default (for the first replicate):
Input:
  --phred64          qualities are Phred+64
  
 Presets:            
  For --local:
   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

 Alignment:
  --local            local alignment; ends might be soft clipped (off)

 Output:
  --met-stderr       send metrics to stderr (off)

 Performance:
  --threads  20


- Control-1st-R1-trimmed.fastq = 11.055.926 sequences
    bowtie2 --phred64 --local --very-sensitive-local --met-stderr --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -U /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Control-1st-R1-trimmed.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-1st-bowtie2.sam

11055926 reads; of these:

    11055926 (100.00%) were unpaired; of these:
        731674 (6.62%) aligned 0 times
        7683346 (69.50%) aligned exactly 1 time
        2640906 (23.89%) aligned >1 times

93.38% overall alignment rate



- Resistant-1st-R1-trimmed.fastq = 11.279.633 sequences
    bowtie2 --phred64 --local --very-sensitive-local --met-stderr --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -U /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Resistant-1st-R1-trimmed.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-1st-bowtie2.sam

11279633 reads; of these:

    11279633 (100.00%) were unpaired; of these:
        763190 (6.77%) aligned 0 times
        7889860 (69.95%) aligned exactly 1 time
        2626583 (23.29%) aligned >1 times

93.23% overall alignment rate


b. Second replicate - paired-end (/Volumes/HD3/analyses/Tatiana/Chom-resistance/0-raw_sequences/1st-replicate/)

Options not default (for the first replicate):
 Presets:            
  For --local:
   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

Alignment:
  --local            local alignment; ends might be soft clipped (off)

 Output:
  --met-stderr       send metrics to stderr (off)

Paired-end:
  --maxins 1000  maximum fragment length 
  --no-mixed         suppress unpaired alignments for paired reads
  --no-discordant    suppress discordant alignments for paired reads


 Performance:
  --threads  20


- Control-2nd-R1-trimmed-paired.fastq / Control-2nd-R2-trimmed-paired.fastq = 26.277.392
    bowtie2 --local --very-sensitive-local --met-stderr --maxins 1000 --no-mixed --no-discordant --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -1 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R1-trimmed-paired.fastq -2 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R2-trimmed-paired.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-2nd-bowtie2.sam

26277392 reads; of these:

    26277392 (100.00%) were paired; of these:
        1397364 (5.32%) aligned concordantly 0 times
        18985486 (72.25%) aligned concordantly exactly 1 time
        5894542 (22.43%) aligned concordantly >1 times

94.68% overall alignment rate


- Resistant-2nd-R1-trimmed-paired.fastq / Resistant-2nd-R2-trimmed-paired.fastq= 24.091.747
    bowtie2 --local --very-sensitive-local --met-stderr --maxins 1000 --no-mixed --no-discordant --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -1 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq -2 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R2-trimmed-paired.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-2nd-bowtie2.sam

24091747 reads; of these:

    24091747 (100.00%) were paired; of these:
        1281992 (5.32%) aligned concordantly 0 times
        17760654 (73.72%) aligned concordantly exactly 1 time
        5049101 (20.96%) aligned concordantly >1 times

94.68% overall alignment rate

17.03.2017

## 7 Counting reads 

a. eXpress v1.5.1
Streaming quantification for high-throughput sequencing
http://bio.math.berkeley.edu/eXpress/tutorial.html


- Control-1st
    ./express /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-1st-bowtie2.sam -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Control-1st

Start: 2017-Mar-17 17:11:32 
Finish: 2017-Mar-17 17:19:26


- Resistant-1st
    ./express /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-1st-bowtie2.sam -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Resistant-1st

Start: 2017-Mar-17 17:13:26 
Finish: 2017-Mar-17 17:21:21


- Control-2nd
    ./express /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-2nd-bowtie2.sam -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Control-2nd

Start: 2017-Mar-17 17:15:21 
Finish: 2017-Mar-17 17:37:12


- Resistant-2nd
    ./express /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-2nd-bowtie2.sam -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Resistant-2nd

Start: 2017-Mar-17 17:15:29 
Finish: 2017-Mar-17 17:38:00

b. RSEM 1.2.25
RSEM: accurate quantification of gene and isoform expression from RNA-Seq data http://deweylab.biostat.wisc.edu/rsem/

Building RSEM reference

    ./rsem-prepare-reference /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta Chom-resistant

Checking if the SAM files satisfy the requirements


- Control-1st
    ./rsem-sam-validator /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-1st-bowtie2.sam 

[samopen] SAM header is present: 33038 sequences.
............
The input file is valid!


- Resistant-1st
    ./rsem-sam-validator /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-1st-bowtie2.sam

[samopen] SAM header is present: 33038 sequences.
............
The input file is valid!


- Control-2nd
    ./rsem-sam-validator /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Control-2nd-bowtie2.sam

[samopen] SAM header is present: 33038 sequences.
...........................


- Resistant-2nd
    ./rsem-sam-validator /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/Resistant-2nd-bowtie2.sam 

[samopen] SAM header is present: 33038 sequences.
.........................
The input file is valid!

**RSEM does not support gapped alignments!**

Since currently RSEM does not handle indel, local and discordant alignments, the Bowtie2 parameters are set in a way to avoid those alignments. In particular, we use options '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1' by default. The last parameter of '--score-min', '-0.1', is the negative of maximum mismatch rate. This rate can be set by option '--bowtie2-mismatch-rate'. If reads are paired-end, we additionally use options '--no-mixed' and  '--no-discordant'.

Running Bowtie2 again without gaps:


- Control-1st-R1-trimmed.fastq = 11.055.926 sequences
    bowtie2 --phred64 --sensitive --met-stderr --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -U /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Control-1st-R1-trimmed.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Control-1st-bowtie2-ungapped.sam

11055926 reads; of these:

    11055926 (100.00%) were unpaired; of these:
    1814940 (16.42%) aligned 0 times
    8511452 (76.99%) aligned exactly 1 time
    729534 (6.60%) aligned >1 times

83.58% overall alignment rate


- Resistant-1st-R1-trimmed.fastq = 11.279.633 sequences
    bowtie2 --phred64 --sensitive --met-stderr --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 --threads 20 -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -U /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Resistant-1st-R1-trimmed.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Resistant-1st-bowtie2-ungapped.sam

11279633 reads; of these:

    11279633 (100.00%) were unpaired; of these:
        2294719 (20.34%) aligned 0 times
        8265998 (73.28%) aligned exactly 1 time
        718916 (6.37%) aligned >1 times

79.66% overall alignment rate


- Control-2nd-R1-trimmed-paired.fastq / Control-2nd-R2-trimmed-paired.fastq = 26.277.392
    bowtie2 --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 --threads 20 --no-mixed --no-discordant -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -1 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R1-trimmed-paired.fastq -2 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R2-trimmed-paired.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Control-2nd-bowtie2-ungapped.sam 

26277392 reads; of these:

    26277392 (100.00%) were paired; of these:
        6137670 (23.36%) aligned concordantly 0 times
        19048225 (72.49%) aligned concordantly exactly 1 time
        1091497 (4.15%) aligned concordantly >1 times

76.64% overall alignment rate


- Resistant-2nd-R1-trimmed-paired.fastq / Resistant-2nd-R2-trimmed-paired.fastq= 24.091.747
    bowtie2 --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 --threads 20 --no-mixed --no-discordant -x /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Chom-resistant -1 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq -2 /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R2-trimmed-paired.fastq -S /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Resistant-2nd-bowtie2-ungapped.sam


Running RSEM with untapped alignments


- Control-1st
    ./rsem-calculate-expression --sam -p 20 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Control-1st-bowtie2-ungapped.sam Chom-resistant Control-1st 


- Resistant-1st
    ./rsem-calculate-expression --sam -p 20 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Resistant-1st-bowtie2-ungapped.sam Chom-resistant Resistant-1st 


- Control-2nd
    ./rsem-calculate-expression --paired-end --sam -p 20 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Control-2nd-bowtie2-ungapped.sam Chom-resistant Control-2nd 


- Resistant-2nd
    ./rsem-calculate-expression --paired-end --sam -p 20 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/wo_gaps_RSEM/Resistant-2nd-bowtie2-ungapped.sam Chom-resistant Resistant-2nd 

21.03.2017

## 8 Matrices 

*Paste* command to join the rsem.genes.results files side-by-side, then *cut* to select the columns containing the expected_count information, and place them into a final output file. Tail -n+2 is used to print from the second line (without header).

This one-line command assumes the genes (and transcripts) in each files are in the same order. If they are not, the files have to be sorted before joining them together.

**eXpress (estimated and raw counts)**


- Control-1st
    cut -f2,7 /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Control-1st/results.xprs >Control1st.estim.eXpress.txt


- Resistant-1st
    cut -f2,7 /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Resistant-1st/results.xprs >Resistant1st.estim.eXpress.txt


- Control-2nd
    cut -f2,7 /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Control-2nd/results.xprs >Control2nd.estim.eXpress.txt


- Resistant-2nd
    cut -f2,7 /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/eXpress/Resistant-2nd/results.xprs >Resistant2nd.estim.eXpress.txt

Sorting files

    sort /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control1st.genes.eXpress.txt >/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control1st.genes.eXpress.sorted.txt
    sort /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control2nd.genes.eXpress.txt >/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control2nd.genes.eXpress.sorted.txt
    sort /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant1st.genes.eXpress.txt >/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant1st.genes.eXpress.sorted.txt
    sort /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant2nd.genes.eXpress.txt >/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant2nd.genes.eXpress.sorted.txt

Joining files

    join -t $'\t' /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control1st.genes.eXpress.sorted.txt /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Control2nd.genes.eXpress.sorted.txt >temp1.txt
    join -t $'\t' /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant1st.genes.eXpress.sorted.txt /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/estimated/Resistant2nd.genes.eXpress.sorted.txt >temp2.txt
    join -t $'\t' /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/temp1.txt /Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/eXpress-counts/temp2.txt >genes.eXpress.txt

**RSEM**

    paste /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/RSEM/Control-1st/Control-1st.genes.results /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/RSEM/Control-2nd/Control-2nd.genes.results /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/RSEM/Resistant-1st/Resistant-1st.genes.results /Volumes/HD3/analyses/Tatiana/Chom-resistance/7-read_counting/RSEM/Resistant-2nd/Resistant-2nd.genes.results | tail -n+2 | cut -f1,5,12,19,26 > genes.rsem.txt

22.03.2017

## 9 Testing 

9.1 EdgeR 
Version 3.16.5
Empirical Analysis of Digital Gene Expression Data in R

Robinson MD, McCarthy DJ and Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26, pp. -1.



    library(edgeR)
    
    # Reading in the data
    raw.data <- read.table(file = "/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/genes.eXpress.estim.txt", header = FALSE)
    head(raw.data)
    
    # edgeR requires the the dataset to contain only the counts with the row names as the gene ids and the column names as the sample ids
    counts <- raw.data[ , 2:5]
    rownames(counts) <- raw.data[ , 1] # contig names
    colnames(counts) <- c("Control_1st", "Control_2nd","Resistant_1st","Resistant_2nd") # sample names
    
    
    # Summaries
    dim(counts)
    colSums(counts) # Library Sizes
    colSums(counts) / 1e06 # Library Sizes in millions of reads
    table(rowSums(counts))[1:30] # Number of genes with low counts
    
    
    # Building the edgeR Object
    group <- c(rep("C", 2) , rep("R", 2)) #C=control; R=resistant 
    dge <- DGEList(counts, group = group)
    names(dge)
    
    head(dge$counts) # original count matrix
    dge$samples # summary of samples
    
    sum(dge$all.zeros) # How many genes have 0 counts across all samples
    
    # Filtering out low count reads
    # keeping only those genes that have at least 1 read per million in at least 2 samples
    dim(dge)
    keep <- rowSums(cpm(dge)>1) >= 2
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    dim(dge)
    
    # Calculating the normalization factors which correct for the different compositions of the samples
    dge <- calcNormFactors(dge)
    plotSmear(dge)
    dge$samples
    
    # Effective library sizes
    dge$samples$lib.size * dge$samples$norm.factors
    
    # Multi-Dimensional Scaling Plot
    plotMDS(dge, main = "MDS Plot for Count Data", labels = colnames(dge$counts))
    
    
    # Estimating Dispersions
    
    dge <- estimateDisp(dge)
    dge <- estimateCommonDisp(dge)
    dge <- estimateTagwiseDisp(dge)
    
    # Mean-Variance Plot
    
    meanVarPlot <- plotMeanVar( dge , show.raw.vars=TRUE ,
                                show.tagwise.vars=TRUE ,
                                show.binned.common.disp.vars=FALSE ,
                                show.ave.raw.vars=FALSE ,
                                dispersion.method = "qcml" , NBline = TRUE ,
                                nbins = 100 ,
                                pch = 16 ,
                                xlab ="Mean Expression (Log10 Scale)" ,
                                ylab = "Variance (Log10 Scale)" ,
                                main = "Mean-Variance Plot" )
    
    # Testing for DE genes
    results.et <- exactTest(dge)
    
    # summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05
    summary(decideTestsDGE(results.et, p=0.05))
    
    # plotting the log-fold changes of all the genes, and the highlight those that are differentially expressed.
    deGenes <- decideTestsDGE(results.et, p=0.05)
    deGenes <- rownames(results.et)[as.logical(deGenes)]
    plotSmear(results.et, de.tags=deGenes)
    abline(h=c(-1, 1), col=2)
    
    # Writing the results
    
    all.contigs <- topTags(results.et, n=Inf, adjust.method="BH", sort.by="p.value")
    write.table(all.contigs, file="/Volumes/HD3/analyses/Tatiana/Chom-resistance/9-testing/1-EdgeR/allcontigs-eXpress_EdgeR.txt", sep="\t")
    
    options( digits = 4 ) # print only 4 digits
    topTags( results.et , n = 20 , sort.by = "p.value" ) # top 20 DE genes
    
    
    significant.contigs <- topTags( results.et , n = Inf , sort.by = "p.value", p.value = 0.05) 
    
    write.table(significant.contigs, file="/Volumes/HD3/analyses/Tatiana/Chom-resistance/9-testing/1-EdgeR/significant-eXpress_EdgeR.txt", sep="\t")

a. EdgeR + eXpress

Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05:

| Down-regulated | 88    |
| -------------- | ----- |
| No change      | 10389 |
| Up-regulated   | 96    |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.01:

| Down-regulated | 55    |
| -------------- | ----- |
| No change      | 10454 |
| Up-regulated   | 64    |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.001:

| Down-regulated | 23    |
| -------------- | ----- |
| No change      | 10503 |
| Up-regulated   | 47    |


b. EdgeR + RSEM

Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05:

| Down-regulated | 94    |
| -------------- | ----- |
| No change      | 10193 |
| Up-regulated   | 108   |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.01:

| Down-regulated | 58    |
| -------------- | ----- |
| No change      | 10263 |
| Up-regulated   | 74    |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.001:

| Down-regulated | 30    |
| -------------- | ----- |
| No change      | 10314 |
| Up-regulated   | 51    |

  
**Comparison of eXpress and RSEM (both tested in EdgeR, cut-off of 0.05)**

![](https://d2mxuefqeaa7sj.cloudfront.net/s_91E831173CA8F2B570F2F4E249AED6D8451BBCCEF29ED0F5A6798C0E748040C1_1490387244849_venn_result31761.png)



| Names        | Total | Elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------ | ----- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| RSEM eXpress | 178   | c13617_g1_i1 c11928_g1_i2 c12476_g1_i1 c10442_g2_i1 c8704_g1_i2 c13474_g1_i1 c13792_g2_i2 c12989_g2_i3 c9537_g1_i1 c14100_g1_i3 c9103_g1_i1 c8720_g1_i1 c12466_g1_i1 c13581_g1_i2 c9586_g1_i1 c14018_g1_i1 c12978_g1_i1 c12577_g1_i1 c14018_g3_i1 c12174_g1_i1 c29634_g1_i1 c12908_g1_i1 c5372_g1_i1 c7764_g1_i1 c8636_g1_i1 c9771_g1_i2 c14453_g1_i2 c10153_g1_i2 c4804_g1_i1 c13922_g1_i2 c13626_g1_i1 c11717_g1_i1 c11502_g1_i1 c12989_g3_i1 c11961_g1_i1 c13135_g5_i1 c12074_g1_i3 c13331_g1_i1 c11409_g4_i1 c10656_g1_i1 c12192_g3_i4 c10435_g1_i4 c9311_g1_i1 c3260_g1_i1 c12723_g3_i1 c6799_g1_i1 c14544_g1_i4 c9041_g1_i2 c13475_g1_i1 c13820_g1_i7 c13628_g1_i1 c12571_g1_i1 c13959_g1_i1 c13773_g1_i1 c11884_g1_i1 c13135_g4_i1 c11879_g1_i1 c14119_g1_i1 c14015_g4_i2 c12991_g1_i1 c13791_g1_i1 c14257_g1_i1 c8766_g1_i1 c10174_g1_i1 c12563_g1_i1 c14048_g1_i2 c13862_g2_i1 c4195_g1_i1 c13722_g1_i1 c12952_g3_i1 c7864_g1_i1 c8735_g1_i1 c8918_g1_i1 c11167_g1_i1 c12817_g1_i1 c14018_g4_i2 c13437_g3_i1 c13224_g2_i1 c8591_g1_i1 c10463_g1_i1 c9883_g1_i1 c6173_g1_i1 c10762_g1_i1 c14015_g4_i1 c11606_g2_i2 c8280_g1_i1 c9812_g1_i1 c8201_g1_i1 c13820_g1_i1 c29827_g1_i1 c13919_g1_i3 c6588_g1_i1 c6562_g1_i1 c8857_g1_i1 c12649_g1_i1 c14178_g2_i1 c9327_g1_i1 c571_g1_i1 c8899_g1_i1 c10801_g1_i1 c14100_g1_i1 c10504_g1_i1 c8176_g1_i1 c12633_g1_i5 c9754_g1_i1 c14018_g2_i2 c14277_g1_i1 c14195_g1_i1 c13759_g1_i1 c11348_g2_i1 c14139_g1_i1 c14015_g1_i1 c10720_g2_i1 c5697_g1_i1 c11470_g1_i1 c11845_g2_i5 c13812_g3_i2 c12872_g1_i1 c11191_g1_i2 c7673_g2_i1 c10582_g1_i1 c11898_g1_i1 c14434_g2_i2 c9502_g1_i1 c13293_g1_i5 c14131_g2_i3 c10448_g1_i2 c13135_g6_i1 c11659_g1_i14 c12058_g1_i1 c12858_g1_i1 c11959_g1_i1 c13339_g1_i1 c11659_g1_i12 c7673_g1_i2 c10571_g1_i1 c954_g1_i1 c10512_g1_i1 c11186_g1_i1 c13892_g1_i2 c14178_g2_i10 c11147_g1_i1 c6842_g2_i1 c8899_g2_i1 c11548_g1_i1 c9018_g1_i1 c13705_g1_i1 c16157_g1_i1 c14100_g2_i5 c12789_g1_i1 c9248_g1_i2 c12633_g1_i1 c9354_g1_i1 c11218_g1_i1 c13993_g1_i2 c11820_g1_i1 c11823_g1_i1 c9928_g1_i1 c7227_g2_i1 c11854_g1_i1 c13705_g3_i3 c13221_g1_i1 c9552_g1_i1 c12989_g1_i1 c12862_g1_i1 c5364_g1_i1 c9990_g1_i1 c8614_g1_i1 c12991_g2_i2 c9595_g1_i1 c10654_g1_i1 c14306_g1_i2 c6832_g1_i1 c14682_g1_i2 c5746_g1_i1 c9212_g1_i1 c14015_g2_i3 c7834_g1_i3 |
| eXpress      | 6     | c12633_g2_i1 c13673_g1_i1 c13673_g3_i1 c14258_g4_i1 c12192_g3_i1 c14434_g2_i4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| RSEM         | 24    | c11606_g1_i1 c11366_g1_i1 c13264_g1_i3 c28010_g1_i1 c11290_g1_i1 c12633_g2_i3 c12390_g1_i1 c8708_g1_i1 c12626_g1_i1 c3327_g1_i1 c8570_g1_i1 c14015_g5_i1 c7238_g1_i1 c9745_g1_i1 c13597_g1_i1 c6284_g1_i1 c13807_g1_i1 c13085_g1_i1 c11874_g1_i1 c13547_g4_i1 c14015_g4_i3 c12819_g1_i1 c14081_g1_i7 c9044_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |


03.05.2017

9.2 DEseq2 
Version 1.14.1
Differential gene expression analysis based on the negative binomial distribution

Michael I Love, Wolfgang Huber, Simon Anders: Moderated estimation of fold change and disper- sion for RNA-seq data with DESeq2. Genome Biology 2014, 15:550. http://dx.doi.org/10.1186/s13059-014-0550-8 

 Good tutorial with lots of plots: https://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
 

     library('DESeq2')
    
    # Reading in the data
    raw.data <- read.table(file = "/Volumes/HD3/analyses/Tatiana/Chom-resistance/8-matrices/genes.eXpress.estim.txt", header = FALSE,row.names=1)
    head(raw.data)
    
    barplot(colSums(raw.data)*1e-6,names=colnames(raw.data), ylab="Library size (millions)")
    
    # create experiment labels (two conditions)
    colData <- DataFrame(condition=factor(c("C","C","R", "R")))
    
    # create DESeq input matrix
    dds <- DESeqDataSetFromMatrix(raw.data, colData,formula(~ condition))
    
    
    # run DEseq
    dds <- DESeq(dds)
    
    # visualize differentially expressed genes
    plotMA(dds)
    
    # get differentially expressed genes
    res <- results(dds)
    summary(res)
    # summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05
    sum(res$padj < 0.05, na.rm=TRUE )
    
    # order by BH adjusted p-value
    resOrdered <- res[order(res$padj),]
    
    # get differentially expressed gene matrix
    sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.05,]
    
    write.table(sig, file="/Volumes/HD3/analyses/Tatiana/Chom-resistance/9-testing/2-DEseq2/allcontigs-eXpress_DEseq2.txt",sep="\t")

 
a. DEseq2 + eXpress

Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05:

| Down-regulated | 92  |
| -------------- | --- |
| Up-regulated   | 103 |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.01:

| Down-regulated | 62 |
| -------------- | -- |
| Up-regulated   | 69 |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.001:

| Down-regulated | 29 |
| -------------- | -- |
| Up-regulated   | 47 |


b. DEseq2 + RSEM

Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.05:

| Down-regulated | 91  |
| -------------- | --- |
| Up-regulated   | 101 |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.01:

| Down-regulated | 65 |
| -------------- | -- |
| Up-regulated   | 74 |

  
Summary of how many genes are called DE at a false discovery rate (FDR) cut-off of 0.001:

| Down-regulated | 29 |
| -------------- | -- |
| Up-regulated   | 51 |

  
**Comparison of eXpress and RSEM (both tested in DEseq2, cut-off of 0.05)**

![](https://d2mxuefqeaa7sj.cloudfront.net/s_91E831173CA8F2B570F2F4E249AED6D8451BBCCEF29ED0F5A6798C0E748040C1_1493830361055_venn_result20330.png)

| Names        | Total | Elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| ------------ | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| RSEM eXpress | 182   | c13617_g1_i1 c12476_g1_i1 c11928_g1_i2 c10442_g2_i1 c13474_g1_i1 c13792_g2_i2 c12989_g2_i3 c9537_g1_i1 c14100_g1_i3 c13163_g4_i1 c9103_g1_i1 c13581_g1_i2 c12466_g1_i1 c12978_g1_i1 c12577_g1_i1 c14018_g3_i1 c12174_g1_i1 c6842_g1_i1 c12908_g1_i1 c7764_g1_i1 c13806_g4_i1 c8636_g1_i1 c11053_g1_i1 c9771_g1_i2 c14453_g1_i2 c10153_g1_i2 c13922_g1_i2 c13626_g1_i1 c11606_g1_i1 c11717_g1_i1 c11502_g1_i1 c11920_g2_i2 c13135_g5_i1 c11961_g1_i1 c12074_g1_i3 c13331_g1_i1 c10656_g1_i1 c10435_g1_i4 c3260_g1_i1 c14544_g1_i4 c9041_g1_i2 c13475_g1_i1 c13820_g1_i7 c12633_g2_i1 c13628_g1_i1 c13959_g1_i1 c9043_g1_i1 c11884_g1_i1 c13773_g1_i1 c11879_g1_i1 c14119_g1_i1 c14015_g4_i2 c13791_g1_i1 c14257_g1_i1 c12633_g2_i3 c8766_g1_i1 c10174_g1_i1 c12563_g1_i1 c14048_g1_i2 c13513_g1_i1 c4195_g1_i1 c13862_g2_i1 c8998_g1_i1 c13722_g1_i1 c12952_g3_i1 c8735_g1_i1 c8918_g1_i1 c11167_g1_i1 c12920_g1_i7 c12817_g1_i1 c14018_g4_i2 c13224_g2_i1 c13437_g3_i1 c8591_g1_i1 c12390_g1_i1 c10463_g1_i1 c9883_g1_i1 c6173_g1_i1 c10762_g1_i1 c14015_g4_i1 c8708_g1_i1 c8280_g1_i1 c11606_g2_i2 c12626_g1_i1 c8201_g1_i1 c13820_g1_i1 c9899_g1_i1 c9882_g1_i1 c29827_g1_i1 c13919_g1_i3 c13163_g4_i2 c6562_g1_i1 c8857_g1_i1 c12649_g1_i1 c13103_g3_i1 c14178_g2_i1 c9327_g1_i1 c571_g1_i1 c7238_g1_i1 c8899_g1_i1 c12515_g1_i1 c13308_g1_i1 c10801_g1_i1 c14100_g1_i1 c10504_g1_i1 c12633_g1_i5 c12661_g1_i1 c14018_g2_i2 c9754_g1_i1 c14277_g1_i1 c12658_g5_i4 c14195_g1_i1 c13759_g1_i1 c11348_g2_i1 c14139_g1_i1 c14015_g1_i1 c10720_g2_i1 c13812_g3_i2 c11845_g2_i5 c12872_g1_i1 c11191_g1_i2 c7673_g2_i1 c10582_g1_i1 c12770_g1_i1 c13807_g1_i1 c14228_g1_i1 c13293_g1_i5 c9502_g1_i1 c14185_g1_i1 c10448_g1_i2 c14131_g2_i3 c11659_g1_i14 c13019_g1_i3 c13339_g1_i1 c11959_g1_i1 c11659_g1_i12 c7673_g1_i2 c12041_g3_i1 c13547_g4_i1 c954_g1_i1 c10512_g1_i1 c11865_g1_i1 c8610_g1_i1 c14198_g1_i1 c13892_g1_i2 c14178_g2_i10 c11147_g1_i1 c14033_g1_i1 c13934_g1_i1 c6842_g2_i1 c8899_g2_i1 c11548_g1_i1 c9018_g1_i1 c14100_g2_i5 c12789_g1_i1 c9641_g1_i1 c12819_g1_i1 c9248_g1_i2 c12633_g1_i1 c9354_g1_i1 c8345_g1_i2 c11218_g1_i1 c11820_g1_i1 c13993_g1_i2 c11823_g1_i1 c9928_g1_i1 c7227_g2_i1 c11854_g1_i1 c13221_g1_i1 c9552_g1_i1 c12989_g1_i1 c12862_g1_i1 c8614_g1_i1 c12991_g2_i2 c9595_g1_i1 c10654_g1_i1 c6832_g1_i1 c13806_g4_i3 c14682_g1_i2 c9212_g1_i1 c14015_g2_i3 c7834_g1_i3 |
| eXpress      | 13    | c9586_g1_i1 c8877_g1_i1 c12657_g1_i1 c12723_g3_i1 c13673_g1_i1 c13673_g3_i1 c14258_g4_i1 c12361_g2_i3 c13421_g4_i1 c9077_g1_i1 c13182_g1_i1 c13705_g1_i1 c14606_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| RSEM         | 10    | c29634_g1_i1 c13264_g1_i3 c12192_g3_i4 c9311_g1_i1 c14364_g1_i1 c12991_g1_i1 c3327_g1_i1 c9812_g1_i1 c13163_g3_i1 c5697_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |


**Comparison of eXpress-EdgeR, RSEM-EdgeR, eXpress-DEseq2 e RSEM-DEseq2 (cut-off of 0.05)**
http://www.interactivenn.net


![](https://d2mxuefqeaa7sj.cloudfront.net/s_91E831173CA8F2B570F2F4E249AED6D8451BBCCEF29ED0F5A6798C0E748040C1_1493844328091_venn-all.png)



| Names                                               | Total | Elements                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| --------------------------------------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| RSEM-DEseq2 RSEM-EdgeR eXpress-DEseq2 eXpress-EdgeR | 142   | c13617_g1_i1 c12476_g1_i1 c11928_g1_i2 c10442_g2_i1 c13474_g1_i1 c13792_g2_i2 c12989_g2_i3 c9537_g1_i1 c14100_g1_i3 c9103_g1_i1 c13581_g1_i2 c12466_g1_i1 c12978_g1_i1 c12577_g1_i1 c14018_g3_i1 c12174_g1_i1 c12908_g1_i1 c7764_g1_i1 c8636_g1_i1 c9771_g1_i2 c14453_g1_i2 c10153_g1_i2 c13922_g1_i2 c13626_g1_i1 c11717_g1_i1 c11502_g1_i1 c13135_g5_i1 c11961_g1_i1 c12074_g1_i3 c13331_g1_i1 c10656_g1_i1 c10435_g1_i4 c3260_g1_i1 c14544_g1_i4 c9041_g1_i2 c13475_g1_i1 c13820_g1_i7 c13628_g1_i1 c13959_g1_i1 c11884_g1_i1 c13773_g1_i1 c11879_g1_i1 c14119_g1_i1 c14015_g4_i2 c13791_g1_i1 c14257_g1_i1 c8766_g1_i1 c10174_g1_i1 c12563_g1_i1 c14048_g1_i2 c4195_g1_i1 c13862_g2_i1 c13722_g1_i1 c12952_g3_i1 c8735_g1_i1 c8918_g1_i1 c11167_g1_i1 c12817_g1_i1 c14018_g4_i2 c13224_g2_i1 c13437_g3_i1 c8591_g1_i1 c10463_g1_i1 c9883_g1_i1 c6173_g1_i1 c10762_g1_i1 c14015_g4_i1 c8280_g1_i1 c11606_g2_i2 c8201_g1_i1 c13820_g1_i1 c29827_g1_i1 c13919_g1_i3 c6562_g1_i1 c8857_g1_i1 c12649_g1_i1 c14178_g2_i1 c9327_g1_i1 c571_g1_i1 c8899_g1_i1 c10801_g1_i1 c14100_g1_i1 c10504_g1_i1 c12633_g1_i5 c14018_g2_i2 c9754_g1_i1 c14277_g1_i1 c14195_g1_i1 c13759_g1_i1 c11348_g2_i1 c14139_g1_i1 c14015_g1_i1 c10720_g2_i1 c13812_g3_i2 c11845_g2_i5 c12872_g1_i1 c11191_g1_i2 c7673_g2_i1 c10582_g1_i1 c13293_g1_i5 c9502_g1_i1 c10448_g1_i2 c14131_g2_i3 c11659_g1_i14 c13339_g1_i1 c11959_g1_i1 c11659_g1_i12 c7673_g1_i2 c954_g1_i1 c10512_g1_i1 c13892_g1_i2 c14178_g2_i10 c11147_g1_i1 c6842_g2_i1 c8899_g2_i1 c11548_g1_i1 c9018_g1_i1 c14100_g2_i5 c12789_g1_i1 c9248_g1_i2 c12633_g1_i1 c9354_g1_i1 c11218_g1_i1 c11820_g1_i1 c13993_g1_i2 c11823_g1_i1 c9928_g1_i1 c7227_g2_i1 c11854_g1_i1 c13221_g1_i1 c9552_g1_i1 c12989_g1_i1 c12862_g1_i1 c8614_g1_i1 c12991_g2_i2 c9595_g1_i1 c10654_g1_i1 c6832_g1_i1 c14682_g1_i2 c9212_g1_i1 c14015_g2_i3 c7834_g1_i3 |
| RSEM-DEseq2 eXpress-DEseq2                          | 30    | c13163_g4_i1 c6842_g1_i1 c13806_g4_i1 c11053_g1_i1 c11920_g2_i2 c9043_g1_i1 c13513_g1_i1 c8998_g1_i1 c12920_g1_i7 c9899_g1_i1 c9882_g1_i1 c13163_g4_i2 c13103_g3_i1 c12515_g1_i1 c13308_g1_i1 c12661_g1_i1 c12658_g5_i4 c12770_g1_i1 c14228_g1_i1 c14185_g1_i1 c13019_g1_i3 c12041_g3_i1 c11865_g1_i1 c8610_g1_i1 c14198_g1_i1 c14033_g1_i1 c13934_g1_i1 c9641_g1_i1 c8345_g1_i2 c13806_g4_i3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| RSEM-EdgeR eXpress-EdgeR                            | 27    | c8704_g1_i2 c8720_g1_i1 c14018_g1_i1 c5372_g1_i1 c4804_g1_i1 c12989_g3_i1 c11409_g4_i1 c6799_g1_i1 c12571_g1_i1 c13135_g4_i1 c7864_g1_i1 c6588_g1_i1 c8176_g1_i1 c11470_g1_i1 c11898_g1_i1 c14434_g2_i2 c13135_g6_i1 c12058_g1_i1 c12858_g1_i1 c10571_g1_i1 c11186_g1_i1 c16157_g1_i1 c13705_g3_i3 c5364_g1_i1 c9990_g1_i1 c14306_g1_i2 c5746_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| RSEM-DEseq2 eXpress-DEseq2 eXpress-EdgeR            | 1     | c12633_g2_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| RSEM-DEseq2 RSEM-EdgeR eXpress-DEseq2               | 9     | c11606_g1_i1 c12633_g2_i3 c12390_g1_i1 c8708_g1_i1 c12626_g1_i1 c7238_g1_i1 c13807_g1_i1 c13547_g4_i1 c12819_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| RSEM-EdgeR eXpress-DEseq2 eXpress-EdgeR             | 3     | c9586_g1_i1 c12723_g3_i1 c13705_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| RSEM-DEseq2 RSEM-EdgeR eXpress-EdgeR                | 6     | c29634_g1_i1 c12192_g3_i4 c9311_g1_i1 c12991_g1_i1 c9812_g1_i1 c5697_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| eXpress-DEseq2 eXpress-EdgeR                        | 3     | c13673_g1_i1 c13673_g3_i1 c14258_g4_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| RSEM-DEseq2 RSEM-EdgeR                              | 2     | c13264_g1_i3 c3327_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| eXpress-DEseq2                                      | 7     | c8877_g1_i1 c12657_g1_i1 c12361_g2_i3 c13421_g4_i1 c9077_g1_i1 c13182_g1_i1 c14606_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| RSEM-DEseq2                                         | 2     | c14364_g1_i1 c13163_g3_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| eXpress-EdgeR                                       | 2     | c12192_g3_i1 c14434_g2_i4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| RSEM-EdgeR                                          | 13    | c11366_g1_i1 c28010_g1_i1 c11290_g1_i1 c8570_g1_i1 c14015_g5_i1 c9745_g1_i1 c13597_g1_i1 c6284_g1_i1 c13085_g1_i1 c11874_g1_i1 c14015_g4_i3 c14081_g1_i7 c9044_g1_i1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |


**Combining info of the 142 contigs common to all conditions** 


- Sorting files by contain name (and removing quotation marks)
    sort significant-eXpress_EdgeR.txt | tr -d '\"' >significant-eXpress_EdgeR.sorted
    sort significant-RSEM_EdgeR.txt | tr -d '\"' >significant-RSEM_EdgeR.sorted
    
    sort significant-eXpress_DEseq2-0.05.txt | tr -d '\"' >significant-eXpress_DEseq2.sorted
    sort significant-RSEM_DEseq2-0.05.txt | tr -d '\"' >significant-RSEM_DEseq2.sorted


- Selecting common contigs
    join -t $'\t' significant-eXpress_EdgeR.sorted commonAll.sorted >significant-eXpress_EdgeR.common
    join -t $'\t' significant-RSEM_EdgeR.sorted commonAll.sorted >significant-RSEM_EdgeR.common
    
    join -t $'\t' significant-eXpress_DEseq2.sorted commonAll.sorted >significant-eXpress_DEseq2.common
    join -t $'\t' significant-RSEM_DEseq2.sorted commonAll.sorted >significant-RSEM_DEseq2.common


- Selecting columns with fold-change and FDR
    cut -f1,2,5 significant-eXpress_EdgeR.common >significant-eXpress_EdgeR.fc
    cut -f1,2,5 significant-RSEM_EdgeR.common >significant-RSEM_EdgeR.fc
    
    cut -f1,3,7 significant-eXpress_DEseq2.common >significant-eXpress_DEseq2.fc
    cut -f1,3,7 significant-RSEM_DEseq2.common >significant-RSEM_DEseq2.fc


- Joining files
    join -t $'\t' significant-eXpress_EdgeR.fc significant-RSEM_EdgeR.fc >significant-EdgeR.txt
    
    join -t $'\t' significant-eXpress_DEseq2.fc significant-RSEM_DEseq2.fc >significant-DEseq2.txt
    
    join -t $'\t' significant-EdgeR.txt significant-DEseq2.txt >DEcontigsAll.txt


Up-regulated (78 contigs) I think it’s 77 (Sophie)

| Contig        | FC-eXpress_EdgeR | FC-RSEM_EdgeR | FC-eXpress_DEseq2 | FC-RSEM_DEseq2 |
| ------------- | ---------------- | ------------- | ----------------- | -------------- |
| c10656_g1_i1  | 6.842016333      | 6.848053696   | 5.44324037        | 5.496869137    |
| c8280_g1_i1   | 6.354426078      | 6.311289576   | 3.810382862       | 3.890100536    |
| c8614_g1_i1   | 6.041746266      | 5.995270617   | 2.876819818       | 2.856814625    |
| c12789_g1_i1  | 5.27921243       | 5.252913723   | 4.586595475       | 4.566364278    |
| c10504_g1_i1  | 4.58302395       | 4.5792317     | 3.703174111       | 3.758130439    |
| c3260_g1_i1   | 4.407978855      | 4.35521916    | 1.806071446       | 1.887298394    |
| c14015_g4_i2  | 4.398174521      | 4.448322872   | 2.749424976       | 2.784446894    |
| c8918_g1_i1   | 4.368850815      | 4.527812839   | 3.509625151       | 3.627241042    |
| c13135_g5_i1  | 4.260463277      | 4.373193049   | 2.505984046       | 2.564553321    |
| c14015_g4_i1  | 4.25894198       | 4.338213989   | 2.1509051         | 2.255464057    |
| c10654_g1_i1  | 4.232272579      | 4.167843704   | 3.091611149       | 3.058014161    |
| c14018_g4_i2  | 4.12437109       | 4.241213073   | 2.554657444       | 2.674802934    |
| c13791_g1_i1  | 4.088881393      | 4.107063758   | 3.302949565       | 3.326764922    |
| c14015_g2_i3  | 4.024866393      | 3.970260643   | 2.630545814       | 2.703958166    |
| c9041_g1_i2   | 3.991991106      | 4.078177425   | 2.954621751       | 3.052515483    |
| c11961_g1_i1  | 3.860642289      | 3.807248051   | 3.076376421       | 3.038784943    |
| c14018_g3_i1  | 3.799724073      | 3.785773411   | 2.392131613       | 2.400647627    |
| c14018_g2_i2  | 3.725633679      | 3.671140513   | 2.48985681        | 2.54252475     |
| c9212_g1_i1   | 3.522560847      | 3.676081542   | 2.66020444        | 2.746544118    |
| c12174_g1_i1  | 3.436925461      | 3.426080704   | 2.320233804       | 2.328211042    |
| c6173_g1_i1   | 3.353629681      | 3.348723605   | 2.775393989       | 2.785646371    |
| c9771_g1_i2   | 3.325166261      | 3.556876504   | 2.403933258       | 2.577028074    |
| c571_g1_i1    | 3.214137099      | 3.111271329   | 1.591677495       | 1.560973328    |
| c7764_g1_i1   | 3.176851554      | 3.079232303   | 2.278489721       | 2.099453708    |
| c9018_g1_i1   | 3.106146736      | 3.525592915   | 1.87523608        | 1.9689105      |
| c11959_g1_i1  | 3.017193842      | 3.049614585   | 2.298296006       | 2.329795271    |
| c9248_g1_i2   | 3.01638065       | 3.004915832   | 1.896414967       | 1.922355781    |
| c6562_g1_i1   | 2.967961667      | 3.19740059    | 1.815309862       | 1.940493154    |
| c9502_g1_i1   | 2.924031408      | 2.959406765   | 2.076693673       | 2.246010055    |
| c14682_g1_i2  | 2.890681897      | 2.886898078   | 2.533102644       | 2.528294374    |
| c11167_g1_i1  | 2.872975409      | 2.850064027   | 2.351276071       | 2.335726882    |
| c10762_g1_i1  | 2.864630014      | 2.867074199   | 2.230243167       | 2.225388174    |
| c11147_g1_i1  | 2.748165707      | 2.832622649   | 1.648429682       | 1.660943267    |
| c11218_g1_i1  | 2.68678754       | 2.605271699   | 1.804981757       | 1.749099426    |
| c13922_g1_i2  | 2.630422752      | 2.93807576    | 2.092946966       | 2.269956896    |
| c10720_g2_i1  | 2.356789627      | 2.30820505    | 1.563243301       | 1.513698816    |
| c13792_g2_i2  | 2.330690786      | 2.333349024   | 1.801001621       | 1.835065403    |
| c11823_g1_i1  | 2.296449341      | 2.267659348   | 1.870684977       | 1.845094582    |
| c11820_g1_i1  | 2.24105831       | 2.176717079   | 1.552469883       | 1.533272614    |
| c14015_g1_i1  | 2.229219242      | 2.226038548   | 1.983903154       | 1.989077376    |
| c12476_g1_i1  | 2.22896614       | 2.217960017   | 1.745124183       | 1.743632677    |
| c8636_g1_i1   | 2.226170113      | 2.199731482   | 1.530818798       | 1.53308203     |
| c12563_g1_i1  | 2.179513971      | 2.185316893   | 1.79116345        | 1.795536366    |
| c9928_g1_i1   | 2.090368723      | 2.5821064     | 1.396670426       | 1.651869912    |
| c8591_g1_i1   | 2.077216115      | 2.123022943   | 1.402254462       | 1.440956507    |
| c13820_g1_i7  | 2.015995175      | 2.018908625   | 1.838539764       | 1.841297233    |
| c13293_g1_i5  | 2.014615031      | 1.954985814   | 1.790075274       | 1.731645895    |
| c10801_g1_i1  | 1.924606513      | 1.906487507   | 1.618844938       | 1.622778607    |
| c12074_g1_i3  | 1.897240393      | 1.977974383   | 1.543289169       | 1.61504189     |
| c8201_g1_i1   | 1.877232503      | 1.888575936   | 1.713269467       | 1.715846001    |
| c8766_g1_i1   | 1.874270337      | 1.876627435   | 1.649411778       | 1.652357673    |
| c13820_g1_i1  | 1.862134829      | 1.895644851   | 1.662395822       | 1.699309631    |
| c9754_g1_i1   | 1.814779577      | 1.739580692   | 1.424428082       | 1.386576641    |
| c12466_g1_i1  | 1.795824375      | 1.786972269   | 1.620503942       | 1.612543159    |
| c10153_g1_i2  | 1.736383828      | 1.752771023   | 1.443000922       | 1.462380797    |
| c954_g1_i1    | 1.734855672      | 1.72191237    | 1.542394732       | 1.535310342    |
| c11348_g2_i1  | 1.643434687      | 1.649125154   | 1.463932961       | 1.492059359    |
| c13773_g1_i1  | 1.58652933       | 1.589650912   | 1.441608614       | 1.438101983    |
| c11659_g1_i14 | 1.557983561      | 1.55588869    | 1.327807235       | 1.343818437    |
| c10174_g1_i1  | 1.531252284      | 1.542326603   | 1.238089172       | 1.246924921    |
| c9537_g1_i1   | 1.490642238      | 1.487936586   | 1.278601002       | 1.274141678    |
| c7834_g1_i3   | 1.485381396      | 1.462573467   | 1.214097392       | 1.198780747    |
| c13474_g1_i1  | 1.481407942      | 1.465605734   | 1.271665889       | 1.260271014    |
| c11659_g1_i12 | 1.478306192      | 1.441159168   | 1.273688586       | 1.233592034    |
| c9354_g1_i1   | 1.458111522      | 1.435266485   | 1.234712116       | 1.218293481    |
| c9327_g1_i1   | 1.428491727      | 1.693556538   | 1.248835551       | 1.460189751    |
| c13437_g3_i1  | 1.426562223      | 1.442529615   | 1.302302802       | 1.313017606    |
| c13993_g1_i2  | 1.424976157      | 1.418436557   | 1.282307109       | 1.279499676    |
| c13626_g1_i1  | 1.405673311      | 1.434796628   | 1.17297608        | 1.191242533    |
| c9103_g1_i1   | 1.400235617      | 1.375446946   | 1.256201152       | 1.230260807    |
| c11879_g1_i1  | 1.378827501      | 1.383644181   | 1.160333436       | 1.160159244    |
| c13722_g1_i1  | 1.361049356      | 1.401737349   | 1.145461634       | 1.173715045    |
| c8735_g1_i1   | 1.350337518      | 1.319143119   | 1.196604472       | 1.165625994    |
| c13475_g1_i1  | 1.339257168      | 1.427332553   | 1.199714919       | 1.259317908    |
| c13862_g2_i1  | 1.331925624      | 1.35867799    | 1.165375921       | 1.184502822    |
| c10442_g2_i1  | 1.329826763      | 1.346393736   | 1.177327524       | 1.186481916    |
| c11502_g1_i1  | 1.299695728      | 1.281934764   | 1.132662297       | 1.115871459    |


Down-regulated (65 contigs)

| Contig        | FC-eXpress_EdgeR | FC-RSEM_EdgeR | FC-eXpress_DEseq2 | FC-RSEM_DEseq2 |
| ------------- | ---------------- | ------------- | ----------------- | -------------- |
| c10463_g1_i1  | -4.647911213     | -4.650837614  | -2.577377671      | -2.597380504   |
| c12908_g1_i1  | -4.093852202     | -4.092089388  | -2.774512358      | -2.839744562   |
| c11884_g1_i1  | -3.876257796     | -3.827512007  | -2.663289766      | -2.697446417   |
| c11845_g2_i5  | -3.855873621     | -3.592657982  | -2.356231442      | -2.153797708   |
| c14544_g1_i4  | -3.634831936     | -3.645761854  | -3.257944057      | -3.262375968   |
| c14178_g2_i1  | -3.190677815     | -3.480485993  | -2.378046096      | -2.481154641   |
| c9883_g1_i1   | -2.86774445      | -2.850197273  | -1.993548976      | -2.026394625   |
| c29827_g1_i1  | -2.840205277     | -2.873768907  | -1.837238988      | -1.850474415   |
| c7227_g2_i1   | -2.839010497     | -2.841584019  | -1.754060573      | -1.794388878   |
| c12649_g1_i1  | -2.612330974     | -2.599024547  | -2.205235676      | -2.194709057   |
| c9595_g1_i1   | -2.5988385       | -2.637261541  | -2.040825499      | -2.115297407   |
| c13959_g1_i1  | -2.569090024     | -2.565399729  | -1.847475986      | -1.887635616   |
| c12978_g1_i1  | -2.561182037     | -2.805705543  | -1.672715168      | -1.908891982   |
| c14139_g1_i1  | -2.514938039     | -2.507386733  | -2.218889017      | -2.211905002   |
| c12989_g1_i1  | -2.511936038     | -2.555584107  | -2.066289298      | -2.084453388   |
| c14178_g2_i10 | -2.414607369     | -2.51046318   | -1.92415972       | -1.97189556    |
| c14453_g1_i2  | -2.292898657     | -2.256583841  | -1.482622944      | -1.500351374   |
| c11717_g1_i1  | -2.281021966     | -2.255136588  | -1.555673342      | -1.616259463   |
| c4195_g1_i1   | -2.257216555     | -2.192868272  | -1.457745625      | -1.420531292   |
| c13581_g1_i2  | -2.216994541     | -2.208849853  | -1.731134145      | -1.75200168    |
| c11854_g1_i1  | -2.177535144     | -2.158256416  | -1.792675536      | -1.788818463   |
| c14048_g1_i2  | -2.072363582     | -2.049024534  | -1.66053984       | -1.655403313   |
| c14277_g1_i1  | -2.071168724     | -2.064687455  | -1.652037498      | -1.65241479    |
| c14100_g2_i5  | -2.034399236     | -2.078533356  | -1.552561463      | -1.607284769   |
| c13224_g2_i1  | -2.033895499     | -2.026823907  | -1.380734112      | -1.412359328   |
| c12817_g1_i1  | -1.999028057     | -2.043630649  | -1.474885544      | -1.530709719   |
| c13221_g1_i1  | -1.96275522      | -1.958293392  | -1.490920561      | -1.497839606   |
| c13759_g1_i1  | -1.954015064     | -1.965199921  | -1.456589423      | -1.492316182   |
| c13919_g1_i3  | -1.95167486      | -1.947363737  | -1.430967052      | -1.460581876   |
| c12989_g2_i3  | -1.950604157     | -1.988493615  | -1.525735315      | -1.551197936   |
| c12991_g2_i2  | -1.94756635      | -1.942909528  | -1.67658136       | -1.682382025   |
| c8899_g2_i1   | -1.901182718     | -1.899347752  | -1.541914532      | -1.545045556   |
| c14100_g1_i3  | -1.896383336     | -1.861357612  | -1.542817721      | -1.513930479   |
| c11928_g1_i2  | -1.885715701     | -1.901656701  | -1.592426139      | -1.616613928   |
| c10582_g1_i1  | -1.884639206     | -1.890472023  | -1.467610956      | -1.48577046    |
| c12872_g1_i1  | -1.837547373     | -1.835176495  | -1.454466235      | -1.462360498   |
| c8899_g1_i1   | -1.821926499     | -1.815001618  | -1.450221072      | -1.45288903    |
| c12952_g3_i1  | -1.761319175     | -1.803090303  | -1.313711058      | -1.352571997   |
| c14100_g1_i1  | -1.758822518     | -1.769495714  | -1.415031701      | -1.443781104   |
| c7673_g2_i1   | -1.752051716     | -1.65682227   | -1.451570093      | -1.415429315   |
| c8857_g1_i1   | -1.749405142     | -1.740665125  | -1.348046315      | -1.340876065   |
| c10435_g1_i4  | -1.74702165      | -1.74397999   | -1.507491391      | -1.505136002   |
| c11191_g1_i2  | -1.741618505     | -1.731271005  | -1.299642936      | -1.310399196   |
| c13339_g1_i1  | -1.717062988     | -1.74411251   | -1.502206299      | -1.52706913    |
| c14195_g1_i1  | -1.711054253     | -1.735501487  | -1.367852812      | -1.385936159   |
| c11548_g1_i1  | -1.676593308     | -1.682458032  | -1.384353117      | -1.385705202   |
| c10512_g1_i1  | -1.673215601     | -1.657924331  | -1.374754035      | -1.376043568   |
| c11606_g2_i2  | -1.669468894     | -1.648746681  | -1.461629583      | -1.404147845   |
| c13331_g1_i1  | -1.642770248     | -1.624600994  | -1.388413815      | -1.356253463   |
| c9552_g1_i1   | -1.609464178     | -1.470387777  | -1.307989658      | -1.191918412   |
| c10448_g1_i2  | -1.570746167     | -1.578644037  | -1.315327063      | -1.328717984   |
| c12862_g1_i1  | -1.542398935     | -1.577727744  | -1.338939016      | -1.382263765   |
| c6832_g1_i1   | -1.512990675     | -1.517221956  | -1.270759775      | -1.281183243   |
| c7673_g1_i2   | -1.511754281     | -1.70699115   | -1.217799024      | -1.380225443   |
| c14119_g1_i1  | -1.472018509     | -1.466209681  | -1.291736067      | -1.287227074   |
| c13617_g1_i1  | -1.469310564     | -1.481270781  | -1.259452175      | -1.271342291   |
| c13812_g3_i2  | -1.458372935     | -1.544901084  | -1.172648782      | -1.245096086   |
| c14257_g1_i1  | -1.445667088     | -1.449057748  | -1.266579983      | -1.273550719   |
| c12633_g1_i5  | -1.438264144     | -1.387652238  | -1.233251455      | -1.188871918   |
| c14131_g2_i3  | -1.435565156     | -1.446166204  | -1.20067966       | -1.218850894   |
| c12633_g1_i1  | -1.432225167     | -1.473944887  | -1.209357709      | -1.255041968   |
| c13892_g1_i2  | -1.424750636     | -1.430971449  | -1.17604146       | -1.196699947   |
| c13628_g1_i1  | -1.373942894     | -1.372311893  | -1.156282306      | -1.158573543   |
| c6842_g2_i1   | -1.358488245     | -1.442387768  | -1.172952941      | -1.205797322   |
| c12577_g1_i1  | -1.282382108     | -1.302987601  | -1.0908675        | -1.107610848   |

**9.3. Annotation of the DE contigs (17/09/2019)**
To add the annotation of the DE contigs I used the latest run of FunctionAnnotator (see details here: [+📋 Project overview: C hominivorax resistant to OP: 13.-New-annotation](https://paper.dropbox.com/doc/Project-overview-C-hominivorax-resistant-to-OP-13.-New-annotation-OvduZSfLKJwpEGosi9whH#:uid=243027137736772068931539&amp;h2=13.-New-annotation) )


    head DEcontigsAll.csv
    #Component_ID        FC-eXpress_EdgeR        FDR-eXpress_EdgeR        FC-RSEM_EdgeR        FDR-RSEM_EdgeR        FC-eXpress_DEseq2        FDR-eXpress_DEseq2        FC-RSEM_DEseq2        FDR-RSEM_DEseq2
    c10153_g1_i2        1.73638382821132        0.001267601113711        1.75277102337873        0.000771619728752        1.44300092157556        0.002894538825979        1.46238079663274        0.002070847087375
    c10174_g1_i1        1.53125228358235        0.021444135391842        1.54232660265404        0.013035354698779        1.23808917178955        0.032269552279597        1.24692492079945        0.032077207544911
    c10435_g1_i4        -1.74702164990211        0.000575860470239        -1.74397998987043        0.000273227012071        -1.50749139097551        3.55E-05        -1.50513600182504        4.37E-05
    
    LANG=en_EN sort DEcontigsAll.csv > DEcontigsAll.sorted.txt
    
    join -t '        ' DEcontigsAll.sorted.txt ../../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt > DEcontigsAll-annotated.txt

After annotating the DE contigs I separated the up and down-regulated contigs creating two files (“DEcontigs-UP-annotated.txt” and “DEcontigs-DOWN-annotated.txt”)

DEcontigsAll-annotated.txt:

| #Component_ID | FC-eXpress_EdgeR  | FDR-eXpress_EdgeR | FC-RSEM_EdgeR     | FDR-RSEM_EdgeR    | FC-eXpress_DEseq2 | FDR-eXpress_DEseq2 | FC-RSEM_DEseq2    | FDR-RSEM_DEseq2   | length | best_hit_to_nr                                                                                          | hit_length | E-value            | Bit_score   | GO_Biological_Process                                                                                                                                                                                                                                                                                         | GO_Cellular_Component                                                                                         | GO_Molecular_Function                                                                                                                                                                                         | Enzyme                | species                       | genus      | family        |
| ------------- | ----------------- | ----------------- | ----------------- | ----------------- | ----------------- | ------------------ | ----------------- | ----------------- | ------ | ------------------------------------------------------------------------------------------------------- | ---------- | ------------------ | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- | ----------------------------- | ---------- | ------------- |
| c10153_g1_i2  | 1.73638382821132  | 0.001267601113711 | 1.75277102337873  | 0.000771619728752 | 1.44300092157556  | 0.002894538825979  | 1.46238079663274  | 0.002070847087375 | 877    | gi|906470961|gb|KNC31394.1| hypothetical protein FF38_07224, partial                                    | 116        | 6.42E-46           | 190.316104  | GO:0045087 innate immune response | GO:0042742 defense response to bacterium                                                                                                                                                                                                                                  | GO:0005576 extracellular region                                                                               | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10174_g1_i1  | 1.53125228358235  | 0.021444135391842 | 1.54232660265404  | 0.013035354698779 | 1.23808917178955  | 0.032269552279597  | 1.24692492079945  | 0.032077207544911 | 1621   | gi|906474519|gb|KNC34073.1| putative inorganic phosphate cotransporter                                  | 466        | 8.34E-265          | 905.164283  | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0005316 high affinity inorganic phosphate:sodium symporter activity                                                                                                                                        | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10435_g1_i4  | -1.74702164990211 | 0.000575860470239 | -1.74397998987043 | 0.000273227012071 | -1.50749139097551 | 3.55E-05           | -1.50513600182504 | 4.37E-05          | 3359   | gi|907638386|ref|XP_013100019.1| PREDICTED: kelch-like ECH-associated protein 1 isoform X1              | 815        | 0                  | 1436.86954  | GO:0006200 ATP catabolic process                                                                                                                                                                                                                                                                              | -                                                                                                             | GO:0003676 nucleic acid binding | GO:0008026 ATP-dependent helicase activity | GO:0005524 ATP binding                                                                                                         | -                     | 35570|Stomoxys calcitrans     | Stomoxys   | Muscidae      |
| c10442_g2_i1  | 1.32982676287699  | 0.03044205084054  | 1.3463937357772   | 0.013087363050763 | 1.17732752368182  | 0.009313496314355  | 1.18648191564633  | 0.009427628827031 | 1783   | gi|906459342|gb|KNC22355.1| hypothetical protein FF38_02930                                             | 519        | 0                  | 1146.931118 | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10448_g1_i2  | -1.57074616704312 | 0.012218960150021 | -1.57864403661609 | 0.00886159587428  | -1.31532706304234 | 0.002177948952858  | -1.32871798428907 | 0.001706498403217 | 1118   | gi|906475417|gb|KNC34791.1| hypothetical protein FF38_01919                                             | 230        | 2.38E-116          | 411.632736  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10463_g1_i1  | -4.6479112132852  | 1.49E-07          | -4.65083761360287 | 3.06E-08          | -2.57737767122561 | 7.91E-07           | -2.59738050373477 | 5.03E-07          | 1380   | gi|262357146|gb|ACY56510.1| yolk protein C                                                              | 423        | 1.86E-234          | 804.276619  | GO:0006629 lipid metabolic process                                                                                                                                                                                                                                                                            | GO:0005576 extracellular region                                                                               | GO:0016787 hydrolase activity                                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10504_g1_i1  | 4.58302394986754  | 1.3E-14           | 4.57923170019501  | 4.29E-16          | 3.70317411114045  | 7.93E-23           | 3.75813043913613  | 3.47E-25          | 1088   | gi|906458465|gb|KNC21660.1| hypothetical protein FF38_09761                                             | 299        | 1.97E-149          | 521.609379  | GO:0045087 innate immune response | GO:0019731 antibacterial humoral response | GO:0050829 defense response to Gram-negative bacterium                                                                                                                                                                        | GO:0005576 extracellular region                                                                               | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10512_g1_i1  | -1.67321560106424 | 0.004942521749819 | -1.65792433081931 | 0.004302546030229 | -1.37475403529575 | 0.002340990648546  | -1.37604356835937 | 0.001749237800033 | 1545   | gi|906463224|gb|KNC25335.1| hypothetical protein FF38_13713                                             | 320        | 2.24E-201          | 694.299976  | GO:0006633 fatty acid biosynthetic process | GO:0055114 oxidation-reduction process | GO:0006703 estrogen biosynthetic process | GO:0010811 positive regulation of cell-substrate adhesion | GO:0030198 extracellular matrix organization                                                                     | GO:0016021 integral to membrane | GO:0031012 extracellular matrix | GO:0005789 endoplasmic reticulum membrane | GO:0001968 fibronectin binding | GO:0005518 collagen binding | GO:0004303 estradiol 17-beta-dehydrogenase activity | GO:0008201 heparin binding                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10582_g1_i1  | -1.88463920595834 | 0.009402608472797 | -1.89047202274665 | 0.009316635406884 | -1.46761095558985 | 0.003500377575132  | -1.48577045970015 | 0.00254869911558  | 574    | gi|906473147|gb|KNC33115.1| Larval cuticle protein 5                                                    | 102        | 1.4E-51            | 206.676265  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0042302 structural constituent of cuticle                                                                                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10654_g1_i1  | 4.23227257921257  | 3.93E-14          | 4.16784370350753  | 2.81E-14          | 3.09161114863552  | 3.72E-12           | 3.0580141612974   | 4.27E-12          | 1782   | gi|755876419|ref|XP_011292711.1| PREDICTED: synaptic vesicle glycoprotein 2B-like, partial              | 529        | 1.52E-232          | 797.914334  | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0022857 transmembrane transporter activity                                                                                                                                                                 | -                     | 7370|Musca domestica          | Musca      | Muscidae      |
| c10656_g1_i1  | 6.8420163333172   | 6.82E-25          | 6.84805369603999  | 2.82E-26          | 5.44324036952456  | 5.03E-51           | 5.49686913731415  | 2.64E-54          | 956    | gi|906459961|gb|KNC22839.1| hypothetical protein FF38_08045                                             | 241        | 1.99E-84           | 305.746134  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | GO:0016020 membrane                                                                                           | GO:0004222 metalloendopeptidase activity | GO:0008270 zinc ion binding                                                                                                                                        | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10720_g2_i1  | 2.35678962651081  | 0.005187297251043 | 2.30820504968566  | 0.005417900023366 | 1.56324330052003  | 0.013145803170227  | 1.51369881648098  | 0.021820676031877 | 717    | gi|906474469|gb|KNC34044.1| hypothetical protein FF38_04538                                             | 209        | 7.11E-108          | 383.456902  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c10762_g1_i1  | 2.86463001369922  | 2.59E-07          | 2.86707419852349  | 1.82E-07          | 2.23024316728343  | 5.57E-07           | 2.22538817352193  | 5.03E-07          | 898    | gi|557783382|ref|XP_005191299.1| PREDICTED: serine proteases 1/2-like                                   | 259        | 1.8E-91            | 328.92303   | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | -                     | 7370|Musca domestica          | Musca      | Muscidae      |
| c10801_g1_i1  | 1.92460651252709  | 0.00474432306111  | 1.90648750733674  | 0.003358803835383 | 1.61884493793088  | 0.000158070477037  | 1.62277860713553  | 7.77E-05          | 503    | gi|2565392|gb|AAB81989.1| cuticle 1                                                                     | 118        | 3.41E-65           | 246.667772  | GO:0008363 larval chitin-based cuticle development                                                                                                                                                                                                                                                            | GO:0031012 extracellular matrix                                                                               | GO:0008010 structural constituent of chitin-based larval cuticle                                                                                                                                              | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11147_g1_i1  | 2.74816570703893  | 0.002222384092205 | 2.83262264904158  | 0.001068069533967 | 1.64842968249299  | 0.011022582110553  | 1.66094326690714  | 0.010757007586208 | 574    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c11167_g1_i1  | 2.87297540901924  | 1.45E-06          | 2.85006402666545  | 1.37E-06          | 2.3512760713851   | 7.21E-09           | 2.33572688194647  | 6.62E-09          | 469    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c11191_g1_i2  | -1.7416185045061  | 0.028402387953984 | -1.7312710045699  | 0.021872968149131 | -1.29964293592066 | 0.034758041265964  | -1.31039919574937 | 0.02982369738602  | 1410   | gi|906459609|gb|KNC22552.1| hypothetical protein FF38_09494                                             | 352        | 5.51E-153          | 533.425052  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11218_g1_i1  | 2.68678754048356  | 0.001267601113711 | 2.60527169860071  | 0.002701404917063 | 1.80498175650016  | 0.001623262443506  | 1.74909942579986  | 0.002580189904817 | 413    | gi|195163507|ref|XP_002022591.1| GL13118                                                                | 36         | 2.19E-05           | 58.980361   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7234|Drosophila persimilis    | Drosophila | Drosophilidae |
| c11348_g2_i1  | 1.6434346865582   | 0.0087827620952   | 1.64912515394882  | 0.006923199234754 | 1.46393296076333  | 0.000190532108377  | 1.49205935895617  | 3.21E-05          | 444    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c11502_g1_i1  | 1.29969572751006  | 0.043866237607411 | 1.2819347639888   | 0.025013685043327 | 1.13266229741526  | 0.019421472936318  | 1.11587145930878  | 0.024814840763319 | 1068   | gi|906457459|gb|KNC20912.1| rRNA 2'-O-methyltransferase fibrillarin                                     | 292        | 2.6E-193           | 667.487488  | GO:0007165 signal transduction | GO:0006508 proteolysis                                                                                                                                                                                                                                                       | GO:0005622 intracellular                                                                                      | GO:0030246 carbohydrate binding | GO:0008237 metallopeptidase activity | GO:0005199 structural constituent of cell wall                                                                                       | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11548_g1_i1  | -1.67659330832313 | 0.009636397384923 | -1.6824580318809  | 0.013867242543383 | -1.38435311723247 | 0.001818283315641  | -1.38570520240692 | 0.002070847087375 | 429    | gi|906469296|gb|KNC30110.1| hypothetical protein FF38_09557                                             | 72         | 1.39E-15           | 92.155133   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11606_g2_i2  | -1.66946889436905 | 0.004838050249902 | -1.64874668094951 | 0.007874156648516 | -1.46162958306366 | 3.28E-05           | -1.40414784549802 | 0.00061713329356  | 280    | gi|906467882|gb|KNC29016.1| hypothetical protein FF38_10141, partial                                    | 50         | 7.42E-23           | 113.968682  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11659_g1_i12 | 1.47830619248597  | 0.015376907443488 | 1.44115916769205  | 0.013066278349485 | 1.27368858592831  | 0.005510111164442  | 1.23359203412523  | 0.009591687898781 | 3602   | gi|906460912|gb|KNC23530.1| Sodium channel protein para                                                 | 1085       | 0                  | 2543.452704 | GO:0034765 regulation of ion transmembrane transport | GO:0035725 sodium ion transmembrane transport                                                                                                                                                                                                          | GO:0001518 voltage-gated sodium channel complex                                                               | GO:0005509 calcium ion binding | GO:0005248 voltage-gated sodium channel activity                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11659_g1_i14 | 1.55798356074398  | 0.010428523549856 | 1.55588869026164  | 0.003560989870399 | 1.3278072354155   | 0.004019039728315  | 1.34381843719807  | 0.002268137927426 | 3364   | gi|906460912|gb|KNC23530.1| Sodium channel protein para                                                 | 1077       | 0                  | 2511.186829 | GO:0034765 regulation of ion transmembrane transport | GO:0035725 sodium ion transmembrane transport                                                                                                                                                                                                          | GO:0001518 voltage-gated sodium channel complex                                                               | GO:0005509 calcium ion binding | GO:0005248 voltage-gated sodium channel activity                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11717_g1_i1  | -2.2810219663095  | 0.001571314698673 | -2.25513658785186 | 0.000541474694804 | -1.55567334151259 | 0.010712244583956  | -1.61625946282467 | 0.004368480181668 | 897    | gi|906466047|gb|KNC27574.1| hypothetical protein FF38_05116                                             | 114        | 1.43E-58           | 228.944263  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11820_g1_i1  | 2.24105831010203  | 0.008777256230807 | 2.17671707886054  | 0.009911022075405 | 1.55246988329057  | 0.01061954453776   | 1.53327261418166  | 0.01104347027141  | 2032   | gi|755871685|ref|XP_005182352.2| PREDICTED: long-chain-fatty-acid--CoA ligase bubblegum-like isoform X1 | 638        | 5.24E-273          | 932.431219  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0016874 ligase activity                                                                                                                                                                                    | 6.2.1.3               | 7370|Musca domestica          | Musca      | Muscidae      |
| c11823_g1_i1  | 2.29644934111081  | 5.96E-06          | 2.26765934840923  | 4.01E-06          | 1.87068497701888  | 2.87E-05           | 1.84509458153431  | 4.25E-05          | 1343   | gi|906469063|gb|KNC29946.1| hypothetical protein FF38_08520                                             | 413        | 3.59E-247          | 846.54037   | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0016772 transferase activity, transferring phosphorus-containing groups                                                                                                                                    | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11845_g2_i5  | -3.85587362090522 | 2.95E-07          | -3.59265798194969 | 1.37E-06          | -2.35623144249481 | 8.25E-06           | -2.15379770830528 | 9.78E-05          | 1040   | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c11854_g1_i1  | -2.17753514439891 | 1.14E-05          | -2.15825641594639 | 2.28E-06          | -1.79267553572661 | 1.19E-05           | -1.7888184634285  | 9.49E-06          | 1537   | gi|906469683|gb|KNC30383.1| Phosphoserine phosphatase                                                   | 269        | 1.44E-168          | 585.232231  | GO:0009790 embryo development | GO:0009555 pollen development | GO:0016311 dephosphorylation | GO:0048364 root development | GO:0006564 L-serine biosynthetic process                                                                                                                                         | -                                                                                                             | GO:0004647 phosphoserine phosphatase activity                                                                                                                                                                 | 3.1.3.3               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11879_g1_i1  | 1.37882750143597  | 0.02494642690869  | 1.38364418072702  | 0.01616826431608  | 1.16033343598966  | 0.034939773901903  | 1.16015924434427  | 0.039574221882453 | 1661   | gi|906475742|gb|KNC34999.1| hypothetical protein FF38_02590                                             | 464        | 0                  | 1086.943858 | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11884_g1_i1  | -3.87625779562993 | 4.6E-08           | -3.8275120067833  | 2.86E-08          | -2.66328976576847 | 3.17E-08           | -2.69744641704708 | 7.88E-09          | 1096   | gi|906459100|gb|KNC22176.1| hypothetical protein FF38_01850                                             | 252        | 1.11E-140          | 492.524647  | GO:0055114 oxidation-reduction process | GO:0006066 alcohol metabolic process                                                                                                                                                                                                                                 | -                                                                                                             | GO:0004022 alcohol dehydrogenase (NAD) activity                                                                                                                                                               | 1.1.1.1               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11928_g1_i2  | -1.88571570106103 | 0.003707025498717 | -1.90165670074054 | 0.002868905300248 | -1.59242613905927 | 7.79E-05           | -1.61661392795585 | 4.21E-05          | 1934   | gi|906457846|gb|KNC21189.1| hypothetical protein FF38_00040                                             | 435        | 1.44E-279          | 954.244768  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11959_g1_i1  | 3.01719384168875  | 2.59E-07          | 3.04961458543776  | 7.2E-08           | 2.29829600613873  | 4.07E-07           | 2.32979527078613  | 1.76E-07          | 1243   | gi|906473970|gb|KNC33672.1| hypothetical protein FF38_08753                                             | 390        | 1.01E-212          | 732.019238  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | GO:0016021 integral to membrane                                                                               | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c11961_g1_i1  | 3.8606422892676   | 1.68E-13          | 3.80724805062225  | 6.03E-14          | 3.07637642079331  | 9.52E-15           | 3.03878494260822  | 1.45E-14          | 1930   | gi|906461223|gb|KNC23749.1| hypothetical protein FF38_13814                                             | 541        | 0                  | 1118.300835 | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0022857 transmembrane transporter activity                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12074_g1_i3  | 1.89724039285057  | 0.002774400911239 | 1.97797438346131  | 0.000656726302869 | 1.54328916903995  | 0.002015653007697  | 1.61504189016007  | 0.000730572507228 | 535    | gi|906457654|gb|KNC21046.1| hypothetical protein FF38_02036                                             | 77         | 1.58E-34           | 153.50574   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0003723 RNA binding                                                                                                                                                                                        | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12174_g1_i1  | 3.43692546081618  | 1.31E-06          | 3.42608070397025  | 9.02E-07          | 2.32023380393659  | 4.48E-06           | 2.32821104209456  | 3.31E-06          | 953    | gi|906472329|gb|KNC32479.1| hypothetical protein FF38_00502                                             | 258        | 2.42E-113          | 401.63486   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12466_g1_i1  | 1.79582437546475  | 0.001668024989638 | 1.78697226867551  | 0.001101764062884 | 1.62050394218953  | 3.04E-06           | 1.61254315893875  | 2.81E-06          | 851    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c12476_g1_i1  | 2.22896613954261  | 0.000184852978234 | 2.21796001722559  | 9.76E-05          | 1.74512418323078  | 0.000372125656799  | 1.74363267723085  | 0.000342718015482 | 954    | gi|906474470|gb|KNC34045.1| hypothetical protein FF38_04539                                             | 281        | 3.94E-175          | 607.04578   | GO:0007596 blood coagulation | GO:0006508 proteolysis                                                                                                                                                                                                                                                         | GO:0005576 extracellular region                                                                               | GO:0004252 serine-type endopeptidase activity | GO:0005509 calcium ion binding                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12563_g1_i1  | 2.17951397106408  | 3.12E-05          | 2.18531689280057  | 8.76E-06          | 1.79116344951914  | 4.55E-05           | 1.7955363658282   | 4.37E-05          | 1539   | gi|906464991|gb|KNC26713.1| hypothetical protein FF38_11405                                             | 466        | 3.37E-240          | 823.363474  | GO:0005975 carbohydrate metabolic process | GO:0006032 chitin catabolic process                                                                                                                                                                                                                               | GO:0005576 extracellular region                                                                               | GO:0004568 chitinase activity | GO:0008061 chitin binding                                                                                                                                                     | 3.2.1.14              | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12577_g1_i1  | -1.28238210771028 | 0.038919980762662 | -1.30298760139294 | 0.018408251636005 | -1.09086749967022 | 0.016283326064414  | -1.10761084833512 | 0.016476564967645 | 1605   | gi|906461528|gb|KNC23994.1| Gram-negative bacteria-binding protein 1                                    | 477        | 3.04E-258          | 883.350734  | GO:0005975 carbohydrate metabolic process                                                                                                                                                                                                                                                                     | -                                                                                                             | GO:0004553 hydrolase activity, hydrolyzing O-glycosyl compounds                                                                                                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12633_g1_i1  | -1.43222516702794 | 0.026504872874806 | -1.47394488683384 | 0.011728818922672 | -1.20935770915742 | 0.006803223590521  | -1.25504196849314 | 0.003241513996835 | 689    | gi|906466050|gb|KNC27577.1| hypothetical protein FF38_05113                                             | 150        | 6.12E-62           | 238.487691  | GO:0007165 signal transduction | GO:0006508 proteolysis | GO:0006412 translation | GO:0051301 cell division                                                                                                                                                                                                   | GO:0005840 ribosome                                                                                           | GO:0003735 structural constituent of ribosome | GO:0004197 cysteine-type endopeptidase activity | GO:0000166 nucleotide binding                                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12633_g1_i5  | -1.43826414351641 | 0.036744923029357 | -1.38765223763727 | 0.026752176777488 | -1.23325145490166 | 0.003323455070937  | -1.18887191783962 | 0.005671094839422 | 607    | gi|906466050|gb|KNC27577.1| hypothetical protein FF38_05113                                             | 181        | 1.12E-74           | 275.752504  | GO:0006810 transport | GO:0006412 translation | GO:0051301 cell division                                                                                                                                                                                                                                      | GO:0005840 ribosome                                                                                           | GO:0003735 structural constituent of ribosome | GO:0003824 catalytic activity | GO:0000166 nucleotide binding                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12649_g1_i1  | -2.61233097428549 | 2.56E-08          | -2.59902454659483 | 6.5E-09           | -2.20523567631098 | 1.73E-09           | -2.19470905692974 | 1.73E-09          | 2159   | gi|906465294|gb|KNC26971.1| hypothetical protein FF38_09256                                             | 621        | 0                  | 1408.693706 | GO:0055114 oxidation-reduction process | GO:0006066 alcohol metabolic process                                                                                                                                                                                                                                 | -                                                                                                             | GO:0050660 flavin adenine dinucleotide binding | GO:0008812 choline dehydrogenase activity                                                                                                                    | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12789_g1_i1  | 5.27921243021748  | 6.33E-23          | 5.25291372331722  | 1.57E-24          | 4.58659547549819  | 1.82E-48           | 4.56636427813175  | 5.91E-49          | 1479   | gi|906457223|gb|KNC20706.1| hypothetical protein FF38_07731, partial                                    | 438        | 7.27E-178          | 616.134758  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12817_g1_i1  | -1.99902805680956 | 0.002504728998598 | -2.04363064886408 | 0.000998236221428 | -1.47488554398826 | 0.00893872931846   | -1.53070971929712 | 0.004716961216436 | 1866   | gi|906457643|gb|KNC21035.1| hypothetical protein FF38_02010                                             | 539        | 6.71E-285          | 971.968277  | GO:0008218 bioluminescence | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                           | -                                                                                                             | GO:0016874 ligase activity | GO:0004497 monooxygenase activity                                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12862_g1_i1  | -1.54239893546546 | 0.022605366684263 | -1.5777277444163  | 0.015425999381651 | -1.33893901643247 | 0.00047820808786   | -1.38226376543051 | 0.000138585741444 | 725    | gi|906473940|gb|KNC33642.1| Larval cuticle protein 1                                                    | 138        | 7.39E-49           | 198.950633  | GO:0006810 transport                                                                                                                                                                                                                                                                                          | GO:0005576 extracellular region                                                                               | GO:0042302 structural constituent of cuticle | GO:0005215 transporter activity                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12872_g1_i1  | -1.83754737254902 | 0.004287831299914 | -1.83517649529135 | 0.003231092380711 | -1.45446623545395 | 0.002870967405796  | -1.46236049780855 | 0.00238664048692  | 1047   | gi|906473233|gb|KNC33182.1| hypothetical protein FF38_06781                                             | 278        | 3.94E-175          | 607.04578   | GO:0006807 nitrogen compound metabolic process                                                                                                                                                                                                                                                                | GO:0005737 cytoplasm                                                                                          | GO:0003676 nucleic acid binding | GO:0050152 omega-amidase activity                                                                                                                                           | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12908_g1_i1  | -4.09385220208223 | 4.42E-07          | -4.09208938804777 | 6.7E-07           | -2.77451235837505 | 7.42E-09           | -2.83974456153358 | 1.25E-09          | 1830   | gi|659106067|gb|AID61420.1| cytochrome P450                                                             | 517        | 0                  | 1178.288095 | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding  | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c12952_g3_i1  | -1.76131917505639 | 0.008777256230807 | -1.80309030251167 | 0.004302546030229 | -1.31371105793172 | 0.030475981395372  | -1.35257199704814 | 0.02333346907259  | 1877   | gi|906461796|gb|KNC24209.1| putative fatty acyl-CoA reductase                                           | 547        | 0                  | 1241.456497 | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0080019 fatty-acyl-CoA reductase (alcohol-forming) activity                                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12978_g1_i1  | -2.5611820372452  | 0.000554158250296 | -2.80570554264687 | 1.74E-05          | -1.67271516831545 | 0.006077343650801  | -1.90889198192614 | 0.000536204438525 | 1572   | gi|906474856|gb|KNC34344.1| hypothetical protein FF38_05970                                             | 423        | 8.37E-246          | 841.995881  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0004181 metallocarboxypeptidase activity | GO:0008270 zinc ion binding                                                                                                                                     | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12989_g1_i1  | -2.51193603757056 | 4.45E-07          | -2.55558410668004 | 1.46E-07          | -2.06628929759133 | 1.63E-07           | -2.08445338759086 | 1.97E-07          | 401    | gi|906462496|gb|KNC24741.1| hypothetical protein FF38_05420                                             | 103        | 5.68E-41           | 172.592595  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12989_g2_i3  | -1.9506041565375  | 0.00080271706146  | -1.98849361518819 | 0.000365478972256 | -1.52573531485964 | 0.002084937640004  | -1.55119793619879 | 0.001891875173476 | 295    | gi|906462495|gb|KNC24740.1| hypothetical protein FF38_05421                                             | 86         | 5.1E-37            | 158.050229  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c12991_g2_i2  | -1.94756634979553 | 0.001177675043686 | -1.94290952798235 | 0.00136593506892  | -1.6765813597641  | 5.26E-06           | -1.68238202535246 | 3.07E-06          | 1448   | gi|906457847|gb|KNC21190.1| hypothetical protein FF38_00041                                             | 396        | 5.09E-241          | 826.090168  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0008270 zinc ion binding | GO:0003723 RNA binding                                                                                                                                                          | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13135_g5_i1  | 4.2604632772024   | 2.76E-06          | 4.3731930492462   | 3.05E-06          | 2.50598404624386  | 2.14E-06           | 2.56455332146505  | 7.88E-07          | 315    | gi|63053874|gb|AAY28732.1| heat shock protein 70                                                        | 104        | 8.16E-65           | 240.305487  | GO:0006950 response to stress                                                                                                                                                                                                                                                                                 | -                                                                                                             | GO:0005524 ATP binding                                                                                                                                                                                        | -                     | 265456|Delia antiqua          | Delia      | Anthomyiidae  |
| c13221_g1_i1  | -1.96275521953643 | 0.001023470647962 | -1.95829339234347 | 0.000970149529814 | -1.49092056114221 | 0.005510111164442  | -1.49783960636639 | 0.004696719540519 | 1008   | gi|906459923|gb|KNC22801.1| hypothetical protein FF38_08046                                             | 268        | 1.4E-136           | 478.891179  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0004222 metalloendopeptidase activity | GO:0008270 zinc ion binding                                                                                                                                        | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13224_g2_i1  | -2.03389549853283 | 0.009430441212556 | -2.02682390737175 | 0.007806428345692 | -1.38073411184095 | 0.038926441522003  | -1.41235932763632 | 0.028412558116738 | 954    | gi|906473981|gb|KNC33683.1| Cytochrome P450 6a9                                                         | 306        | 1.67E-179          | 621.588146  | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13293_g1_i5  | 2.01461503117754  | 0.000270852381324 | 1.95498581371514  | 0.000273227012071 | 1.79007527365045  | 3.34E-07           | 1.7316458948363   | 8.05E-07          | 907    | gi|906475171|gb|KNC34585.1| Attacin-A                                                                   | 232        | 5.09E-130          | 457.07763   | GO:0045087 innate immune response | GO:0019731 antibacterial humoral response | GO:0050829 defense response to Gram-negative bacterium                                                                                                                                                                        | GO:0005576 extracellular region                                                                               | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13331_g1_i1  | -1.64277024754668 | 0.015761548077693 | -1.62460099448003 | 0.020043644548121 | -1.38841381476441 | 0.001372550921626  | -1.35625346252413 | 0.003162222828831 | 1024   | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 251        | 8.41E-71           | 265.754627  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13339_g1_i1  | -1.71706298823001 | 0.001458821494163 | -1.74411251033452 | 0.000771619728752 | -1.5022062986345  | 2E-05              | -1.52706912993805 | 1.35E-05          | 989    | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 290        | 2.91E-99           | 354.826619  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13437_g3_i1  | 1.42656222347162  | 0.029992776887381 | 1.44252961450752  | 0.027580364739616 | 1.30230280203533  | 0.000484109793041  | 1.31301760616324  | 0.000428106643212 | 1302   | gi|557780752|ref|XP_005189997.1| PREDICTED: uncharacterized protein LOC101896308                        | 278        | 2.18E-150          | 524.790522  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0008289 lipid binding                                                                                                                                                                                      | -                     | 7370|Musca domestica          | Musca      | Muscidae      |
| c13474_g1_i1  | 1.48140794172463  | 0.021444135391842 | 1.46560573375345  | 0.013129876036914 | 1.27166588868314  | 0.005813570281122  | 1.26027101371691  | 0.006095931895981 | 2018   | gi|907685852|ref|XP_013107184.1| PREDICTED: synaptic vesicle glycoprotein 2A-like                       | 560        | 1.2E-293           | 1001.053009 | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0022857 transmembrane transporter activity                                                                                                                                                                 | -                     | 35570|Stomoxys calcitrans     | Stomoxys   | Muscidae      |
| c13475_g1_i1  | 1.33925716761007  | 0.025829119356291 | 1.42733255304188  | 0.008362545639131 | 1.19971491919372  | 0.004855070824641  | 1.25931790763778  | 0.00400942740608  | 691    | gi|906473189|gb|KNC33148.1| hypothetical protein FF38_03640, partial                                    | 145        | 4.35E-73           | 271.662464  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13581_g1_i2  | -2.21699454102004 | 0.000964522073052 | -2.20884985330105 | 0.000678937319782 | -1.73113414492996 | 0.000195976376867  | -1.75200167970002 | 8.92E-05          | 1862   | gi|906464994|gb|KNC26716.1| hypothetical protein FF38_11370                                             | 294        | 9.8E-181           | 625.678186  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0015020 glucuronosyltransferase activity                                                                                                                                                                   | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13617_g1_i1  | -1.46931056424297 | 0.026483665494914 | -1.48127078114762 | 0.020495316630294 | -1.25945217478108 | 0.001595775670197  | -1.27134229120869 | 0.00145290237151  | 678    | gi|906464505|gb|KNC26324.1| hypothetical protein FF38_00857, partial                                    | 175        | 5.47E-47           | 193.042797  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13626_g1_i1  | 1.40567331091074  | 0.030002891488797 | 1.43479662830972  | 0.016167806159073 | 1.17297607997482  | 0.032908417446907  | 1.19124253315048  | 0.031934114645856 | 1570   | gi|906472519|gb|KNC32642.1| hypothetical protein FF38_13204, partial                                    | 413        | 1.35E-272          | 931.067872  | GO:0007034 vacuolar transport | GO:0015809 arginine transport | GO:0051453 regulation of intracellular pH | GO:0007040 lysosome organization                                                                                                                                                                  | GO:0005765 lysosomal membrane | GO:0031301 integral to organelle membrane                                     | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13628_g1_i1  | -1.37394289384302 | 0.030015384339287 | -1.37231189282833 | 0.028767213947073 | -1.1562823064293  | 0.009750077194816  | -1.1585735426458  | 0.009427628827031 | 926    | gi|906464455|gb|KNC26277.1| hypothetical protein FF38_01137                                             | 290        | 1.11E-121          | 429.356245  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0030246 carbohydrate binding                                                                                                                                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13722_g1_i1  | 1.36104935550643  | 0.025063057532992 | 1.40173734867557  | 0.013867242543383 | 1.14546163368246  | 0.040915961017952  | 1.17371504465198  | 0.036514027407871 | 2222   | gi|906472083|gb|KNC32278.1| hypothetical protein FF38_01816                                             | 708        | 0                  | 1516.398104 | GO:0031167 rRNA methylation                                                                                                                                                                                                                                                                                   | GO:0005634 nucleus                                                                                            | GO:0008649 rRNA methyltransferase activity                                                                                                                                                                    | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13759_g1_i1  | -1.95401506411934 | 0.004218003918341 | -1.96519992108395 | 0.00171757377885  | -1.45658942309525 | 0.008779709966606  | -1.49231618245308 | 0.005117240976119 | 1520   | gi|906463564|gb|KNC25612.1| hypothetical protein FF38_05787                                             | 431        | 1.09E-235          | 808.366659  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0016772 transferase activity, transferring phosphorus-containing groups                                                                                                                                    | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13773_g1_i1  | 1.58652932960427  | 0.006549224608514 | 1.58965091172308  | 0.004636067222834 | 1.44160861393802  | 3.28E-05           | 1.43810198342998  | 4.44E-05          | 1356   | gi|371536099|gb|AEX33294.1| putative salivary trypsin                                                   | 415        | 1.74E-208          | 717.931321  | GO:0006508 proteolysis                                                                                                                                                                                                                                                                                        | GO:0016021 integral to membrane                                                                               | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | -                     | 13632|Lucilia sericata        | Lucilia    | Calliphoridae |
| c13791_g1_i1  | 4.08888139319194  | 1.04E-14          | 4.10706375778439  | 2.59E-16          | 3.30294956539146  | 5.14E-18           | 3.32676492235971  | 9.84E-19          | 1812   | gi|906458780|gb|KNC21921.1| hypothetical protein FF38_12139                                             | 531        | 3.38999999999e-313 | 1066.039207 | GO:0016311 dephosphorylation                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0004035 alkaline phosphatase activity                                                                                                                                                                      | 3.1.3.1               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13792_g2_i2  | 2.33069078626313  | 0.001890439098349 | 2.33334902358282  | 0.001096936605924 | 1.80100162051945  | 0.00022389978678   | 1.83506540303823  | 8.7E-05           | 1337   | gi|906467500|gb|KNC28725.1| hypothetical protein FF38_09446                                             | 159        | 1.34E-88           | 319.379602  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0030246 carbohydrate binding                                                                                                                                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13812_g3_i2  | -1.45837293532399 | 0.013782452857582 | -1.54490108387313 | 0.004078112201383 | -1.1726487824219  | 0.030591453903396  | -1.24509608562165 | 0.017619093172464 | 1857   | gi|906475364|gb|KNC34738.1| Cytochrome P450 4c3                                                         | 531        | 0                  | 1199.192746 | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13820_g1_i1  | 1.86213482887386  | 0.001454265892259 | 1.89564485139335  | 0.001063752954489 | 1.6623958221935   | 2.2E-06            | 1.69930963087497  | 5.06E-07          | 915    | gi|906473923|gb|KNC33625.1| hypothetical protein FF38_08769                                             | 249        | 3.83E-143          | 500.704728  | GO:0010765 positive regulation of sodium ion transport | GO:0006508 proteolysis | GO:0007586 digestion                                                                                                                                                                                                        | GO:0005615 extracellular space | GO:0005886 plasma membrane                                                   | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | 3.4.21.4              | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13820_g1_i7  | 2.01599517523813  | 0.00022917060595  | 2.01890862504708  | 0.000276162987302 | 1.83853976381424  | 3.54E-09           | 1.84129723311076  | 2.19E-09          | 1130   | gi|906473923|gb|KNC33625.1| hypothetical protein FF38_08769                                             | 251        | 9.68E-146          | 509.339258  | GO:0010765 positive regulation of sodium ion transport | GO:0006508 proteolysis | GO:0007586 digestion                                                                                                                                                                                                        | GO:0005615 extracellular space | GO:0005886 plasma membrane                                                   | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | 3.4.21.4              | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13862_g2_i1  | 1.33192562407187  | 0.041097849619341 | 1.35867798994662  | 0.021656446770509 | 1.16537592130789  | 0.010586255955097  | 1.18450282188668  | 0.009427628827031 | 2413   | gi|906457795|gb|KNC21144.1| hypothetical protein FF38_08362                                             | 579        | 0                  | 1287.35584  | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0022857 transmembrane transporter activity                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13892_g1_i2  | -1.42475063587429 | 0.019979847836892 | -1.43097144945038 | 0.008526621060693 | -1.17604146005562 | 0.014960008567893  | -1.19669994738144 | 0.009427628827031 | 1753   | gi|906473183|gb|KNC33142.1| hypothetical protein FF38_03646                                             | 438        | 5.24E-273          | 932.431219  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | GO:0005739 mitochondrion                                                                                      | GO:0047369 succinate-hydroxymethylglutarate CoA-transferase activity                                                                                                                                          | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13919_g1_i3  | -1.95167486038311 | 0.021337431671758 | -1.94736373657219 | 0.015391974872529 | -1.43096705219707 | 0.012478790808125  | -1.46058187633194 | 0.007665023154118 | 2052   | gi|906475254|gb|KNC34646.1| Polypeptide N-acetylgalactosaminyltransferase 8                             | 311        | 3.99E-210          | 723.384708  | GO:0006790 sulfur compound metabolic process                                                                                                                                                                                                                                                                  | GO:0005737 cytoplasm                                                                                          | GO:0004062 aryl sulfotransferase activity | GO:0004027 alcohol sulfotransferase activity                                                                                                                      | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13922_g1_i2  | 2.63042275183663  | 6E-07             | 2.93807575970823  | 4.48E-08          | 2.0929469655755   | 2.43E-06           | 2.26995689580295  | 4.58E-07          | 1967   | gi|906472520|gb|KNC32643.1| hypothetical protein FF38_13205                                             | 553        | 8.61E-278          | 948.336932  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0008484 sulfuric ester hydrolase activity                                                                                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13959_g1_i1  | -2.56909002368909 | 0.000185721432544 | -2.56539972860503 | 9.54E-05          | -1.84747598554116 | 0.00038139696309   | -1.88763561632096 | 0.00015801992326  | 2305   | gi|906472756|gb|KNC32805.1| hypothetical protein FF38_00254                                             | 576        | 0                  | 1187.377073 | GO:0055085 transmembrane transport                                                                                                                                                                                                                                                                            | GO:0016021 integral to membrane                                                                               | GO:0005215 transporter activity                                                                                                                                                                               | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c13993_g1_i2  | 1.42497615668688  | 0.049793718232935 | 1.41843655695081  | 0.048654179025569 | 1.28230710914778  | 0.000702004166093  | 1.27949967622161  | 0.000591312704729 | 2146   | gi|906463872|gb|KNC25845.1| hypothetical protein FF38_08124                                             | 629        | 0                  | 1379.154525 | GO:0006879 cellular iron ion homeostasis | GO:0006826 iron ion transport                                                                                                                                                                                                                                      | GO:0005576 extracellular region                                                                               | GO:0008199 ferric iron binding                                                                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14015_g1_i1  | 2.22921924236142  | 3.5E-05           | 2.22603854788229  | 1.62E-05          | 1.98390315386229  | 7.44E-09           | 1.98907737569431  | 2.79E-09          | 2183   | gi|659105774|gb|AID61342.1| esterase                                                                    | 554        | 8.89E-291          | 991.509581  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0052689 carboxylic ester hydrolase activity                                                                                                                                                                | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c14015_g2_i3  | 4.02486639275931  | 8.36E-07          | 3.97026064256101  | 1.46E-07          | 2.63054581425437  | 1.31E-07           | 2.70395816630924  | 1.65E-08          | 984    | gi|659105774|gb|AID61342.1| esterase                                                                    | 293        | 6.82E-171          | 592.957863  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0052689 carboxylic ester hydrolase activity                                                                                                                                                                | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c14015_g4_i1  | 4.25894197985698  | 0.000183064974892 | 4.33821398921472  | 0.00010018546028  | 2.15090510021991  | 0.000195976376867  | 2.25546405661832  | 5.73E-05          | 947    | gi|906462646|gb|KNC24865.1| hypothetical protein FF38_07854                                             | 203        | 1.02E-117          | 416.177226  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0016787 hydrolase activity                                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14015_g4_i2  | 4.39817452098271  | 5.41E-08          | 4.44832287206977  | 4.57E-08          | 2.74942497560671  | 4.23E-08           | 2.78444689418758  | 1.99E-08          | 886    | gi|659105774|gb|AID61342.1| esterase                                                                    | 269        | 1.7E-157           | 548.421867  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0016787 hydrolase activity                                                                                                                                                                                 | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c14018_g2_i2  | 3.72563367859106  | 6.87E-06          | 3.67114051262213  | 2.8E-06           | 2.48985681030736  | 7.14E-07           | 2.54252475043094  | 1.52E-07          | 879    | gi|119699865|gb|AAD17995.2| 70 kDa heat shock protein ScHSP70                                           | 234        | 3.47E-142          | 497.523586  | GO:0006950 response to stress                                                                                                                                                                                                                                                                                 | -                                                                                                             | GO:0005524 ATP binding                                                                                                                                                                                        | -                     | 59312|Sarcophaga crassipalpis | Sarcophaga | Sarcophagidae |
| c14018_g3_i1  | 3.79972407284079  | 5.57E-06          | 3.78577341124803  | 6.91E-06          | 2.39213161261994  | 4.48E-06           | 2.40064762688616  | 3.31E-06          | 359    | gi|526299367|gb|AGR84223.1| heat shock protein 70-2                                                     | 119        | 7.02E-73           | 267.117974  | GO:0006950 response to stress                                                                                                                                                                                                                                                                                 | -                                                                                                             | GO:0005524 ATP binding                                                                                                                                                                                        | -                     | 113334|Melitaea cinxia        | Melitaea   | Nymphalidae   |
| c14018_g4_i2  | 4.12437108952654  | 6.66E-06          | 4.24121307283065  | 2.12E-06          | 2.55465744384996  | 7.54E-07           | 2.67480293352196  | 1.07E-07          | 531    | gi|907707080|ref|XP_013110576.1| PREDICTED: heat shock protein 68                                       | 177        | 3E-112             | 397.999268  | GO:0006950 response to stress                                                                                                                                                                                                                                                                                 | -                                                                                                             | GO:0005524 ATP binding                                                                                                                                                                                        | -                     | 35570|Stomoxys calcitrans     | Stomoxys   | Muscidae      |
| c14048_g1_i2  | -2.07236358157331 | 6.4E-05           | -2.04902453372323 | 3.72E-05          | -1.66053983996409 | 0.00024208750375   | -1.65540331253694 | 0.000200887815621 | 1586   | gi|906462034|gb|KNC24408.1| hypothetical protein FF38_14164                                             | 328        | 6.88E-225          | 772.465193  | GO:0006783 heme biosynthetic process | GO:0071353 cellular response to interleukin-4 | GO:0051260 protein homooligomerization                                                                                                                                                                                 | -                                                                                                             | GO:0042802 identical protein binding | GO:0032791 lead ion binding | GO:0008270 zinc ion binding | GO:0004655 porphobilinogen synthase activity                                                               | 4.2.1.24              | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14100_g1_i1  | -1.75882251812174 | 0.012218960150021 | -1.76949571396512 | 0.011577009002192 | -1.41503170132394 | 0.00281733067378   | -1.44378110380786 | 0.001462941607984 | 1746   | gi|906462734|gb|KNC24940.1| Protein henna                                                               | 451        | 2.52e-310          | 1056.495779 | GO:0055114 oxidation-reduction process | GO:0006559 L-phenylalanine catabolic process | GO:0071391 cellular response to estrogen stimulus                                                                                                                                                                     | -                                                                                                             | GO:0016597 amino acid binding | GO:0004505 phenylalanine 4-monooxygenase activity | GO:0005506 iron ion binding                                                                                               | 1.14.16.1 | 1.14.16.4 | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14100_g1_i3  | -1.89638333629264 | 0.004999820045418 | -1.86135761235999 | 0.00583414185286  | -1.54281772142091 | 0.000599023643852  | -1.51393047865556 | 0.000913011735624 | 521    | gi|906462734|gb|KNC24940.1| Protein henna                                                               | 130        | 2.15E-74           | 273.480259  | GO:0055114 oxidation-reduction process | GO:0006909 phagocytosis | GO:0006559 L-phenylalanine catabolic process | GO:0007616 long-term memory | GO:0009094 L-phenylalanine biosynthetic process | GO:0071391 cellular response to estrogen stimulus | GO:0006726 eye pigment biosynthetic process             | GO:0005811 lipid particle                                                                                     | GO:0004510 tryptophan 5-monooxygenase activity | GO:0004664 prephenate dehydratase activity | GO:0016597 amino acid binding | GO:0004505 phenylalanine 4-monooxygenase activity | GO:0005506 iron ion binding | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14100_g2_i5  | -2.03439923604265 | 0.009672787245517 | -2.07853335595176 | 0.006305256097815 | -1.5525614633691  | 0.002894538825979  | -1.60728476917581 | 0.001448639879412 | 917    | gi|557769903|ref|XP_005184623.1| PREDICTED: protein henna                                               | 92         | 1.51E-55           | 219.855285  | GO:0055114 oxidation-reduction process | GO:0006559 L-phenylalanine catabolic process | GO:0071391 cellular response to estrogen stimulus                                                                                                                                                                     | -                                                                                                             | GO:0016597 amino acid binding | GO:0004505 phenylalanine 4-monooxygenase activity | GO:0005506 iron ion binding                                                                                               | -                     | 7370|Musca domestica          | Musca      | Muscidae      |
| c14119_g1_i1  | -1.472018508861   | 0.019549543120384 | -1.46620968142371 | 0.014036248655312 | -1.29173606722025 | 0.00024208750375   | -1.28722707445661 | 0.0003076422094   | 1763   | gi|906459360|gb|KNC22368.1| hypothetical protein FF38_08825                                             | 538        | 0                  | 1134.660997 | GO:0009851 auxin biosynthetic process | GO:0009695 jasmonic acid biosynthetic process | GO:0001676 long-chain fatty acid metabolic process                                                                                                                                                                    | GO:0005777 peroxisome                                                                                         | GO:0004321 fatty-acyl-CoA synthase activity | GO:0004467 long-chain fatty acid-CoA ligase activity | GO:0016207 4-coumarate-CoA ligase activity | GO:0008756 o-succinylbenzoate-CoA ligase activity           | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14131_g2_i3  | -1.4355651555142  | 0.045039008432079 | -1.44616620440102 | 0.034568869741593 | -1.20067966043978 | 0.007135445560175  | -1.21885089441973 | 0.005079358717501 | 2256   | gi|755851246|ref|XP_005175488.2| PREDICTED: bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase | 649        | 0                  | 1433.233949 | GO:0016310 phosphorylation | GO:0000103 sulfate assimilation                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0004781 sulfate adenylyltransferase (ATP) activity | GO:0005524 ATP binding | GO:0004020 adenylylsulfate kinase activity                                                                                   | 2.7.1.25 | 2.7.7.4    | 7370|Musca domestica          | Musca      | Muscidae      |
| c14139_g1_i1  | -2.51493803855189 | 7.7E-07           | -2.50738673275192 | 2.09E-07          | -2.21888901682745 | 1.27E-12           | -2.21190500246615 | 1.31E-12          | 2386   | gi|906460284|gb|KNC23077.1| hypothetical protein FF38_10464                                             | 538        | 7.00000000024e-314 | 1068.311452 | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0016874 ligase activity | GO:0004497 monooxygenase activity                                                                                                                                                | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14178_g2_i1  | -3.19067781524265 | 2E-08             | -3.4804859933129  | 4.44E-09          | -2.3780460958031  | 1.79E-07           | -2.48115464067587 | 1.56E-07          | 1880   | gi|906474872|gb|KNC34360.1| hypothetical protein FF38_05930                                             | 268        | 4.82E-158          | 550.239662  | -                                                                                                                                                                                                                                                                                                             | GO:0016459 myosin complex                                                                                     | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14178_g2_i10 | -2.41460736895134 | 9.64E-07          | -2.51046318007865 | 3.47E-07          | -1.92415971953478 | 9.17E-06           | -1.97189555988711 | 9.49E-06          | 2142   | gi|906474872|gb|KNC34360.1| hypothetical protein FF38_05930                                             | 268        | 4.82E-158          | 550.239662  | -                                                                                                                                                                                                                                                                                                             | GO:0016459 myosin complex                                                                                     | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14195_g1_i1  | -1.71105425290883 | 0.001668024989638 | -1.7355014870684  | 0.001063752954489 | -1.36785281238333 | 0.006137383978896  | -1.38593615872327 | 0.005671094839422 | 1681   | gi|906463591|gb|KNC25631.1| putative cytochrome P450 4d14                                               | 503        | 3.82E-292          | 996.054071  | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14257_g1_i1  | -1.44566708787429 | 0.032190953393784 | -1.44905774837578 | 0.029048909017746 | -1.26657998280442 | 0.000842519173365  | -1.27355071917488 | 0.000672369586181 | 937    | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 292        | 1.13E-118          | 419.358368  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14277_g1_i1  | -2.0711687236072  | 6.28E-05          | -2.06468745527091 | 3.98E-05          | -1.65203749773557 | 0.00033313457398   | -1.65241478969104 | 0.00031735171683  | 1723   | gi|906463870|gb|KNC25843.1| hypothetical protein FF38_08180                                             | 520        | 1.39E-304          | 1037.408924 | GO:0016311 dephosphorylation                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0004035 alkaline phosphatase activity                                                                                                                                                                      | 3.1.3.1               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14453_g1_i2  | -2.29289865677168 | 0.008858164148584 | -2.25658384124229 | 0.007925945132529 | -1.48262294371022 | 0.024175903929602  | -1.50035137383714 | 0.019805513208113 | 1828   | gi|659105977|gb|AID61394.1| cytochrome P450                                                             | 498        | 1.64E-274          | 937.430158  | GO:0055114 oxidation-reduction process                                                                                                                                                                                                                                                                        | -                                                                                                             | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding  | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c14544_g1_i4  | -3.63483193566323 | 1.68E-13          | -3.64576185416776 | 2.85E-13          | -3.25794405747177 | 1.11E-31           | -3.26237596802154 | 1.82E-31          | 1906   | gi|906469326|gb|KNC30137.1| hypothetical protein FF38_12389                                             | 546        | 0                  | 1182.378135 | GO:0006212 uracil catabolic process                                                                                                                                                                                                                                                                           | GO:0005737 cytoplasm                                                                                          | GO:0004157 dihydropyrimidinase activity | GO:0051219 phosphoprotein binding                                                                                                                                   | 3.5.2.2               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c14682_g1_i2  | 2.89068189726539  | 6.84E-09          | 2.88689807793142  | 2.45E-09          | 2.53310264426903  | 3.57E-14           | 2.52829437379944  | 2.96E-14          | 4165   | gi|906460053|gb|KNC22899.1| hypothetical protein FF38_00367                                             | 1305       | 0                  | 2130.813068 | GO:0030512 negative regulation of transforming growth factor beta receptor signaling pathway | GO:0010951 negative regulation of endopeptidase activity | GO:0061045 negative regulation of wound healing                                                                                                     | GO:0005615 extracellular space | GO:0005886 plasma membrane | GO:0031225 anchored to membrane                 | GO:0004867 serine-type endopeptidase inhibitor activity | GO:0050431 transforming growth factor beta binding                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c29827_g1_i1  | -2.84020527739403 | 0.000183064974892 | -2.87376890699834 | 0.000137139324197 | -1.83723898807175 | 0.001582628738941  | -1.85047441488812 | 0.00145290237151  | 797    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c3260_g1_i1   | 4.40797885451956  | 0.003146865493979 | 4.355219160043    | 0.002701404917063 | 1.80607144601696  | 0.004548540420174  | 1.8872983941983   | 0.002459085076657 | 982    | gi|906470122|gb|KNC30740.1| hypothetical protein FF38_11653                                             | 276        | 8.67E-110          | 389.819187  | GO:0044070 regulation of anion transport | GO:0015670 carbon dioxide transport | GO:0046903 secretion | GO:0006730 one-carbon metabolic process | GO:2001150 positive regulation of dipeptide transmembrane transport | GO:0002009 morphogenesis of an epithelium | GO:0051453 regulation of intracellular pH | GO:0005829 cytosol | GO:0045177 apical part of cell | GO:0005886 plasma membrane                              | GO:0004089 carbonate dehydratase activity | GO:0008270 zinc ion binding                                                                                                                                       | 4.2.1.1               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c4195_g1_i1   | -2.25721655537508 | 0.005419457774404 | -2.19286827238651 | 0.004811184153724 | -1.45774562514864 | 0.030676220961738  | -1.42053129180736 | 0.042837778092863 | 387    | gi|906466786|gb|KNC28168.1| hypothetical protein FF38_04347                                             | 104        | 8.02E-46           | 187.134961  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c571_g1_i1    | 3.21413709874675  | 0.014975356189602 | 3.11127132864695  | 0.019044907932885 | 1.59167749499018  | 0.020533938201525  | 1.56097332846221  | 0.027520546200771 | 1318   | gi|659105774|gb|AID61342.1| esterase                                                                    | 390        | 2.86E-213          | 733.837033  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0004104 cholinesterase activity                                                                                                                                                                            | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c6173_g1_i1   | 3.35362968121457  | 2.14E-08          | 3.34872360486847  | 6.72E-09          | 2.77539398928248  | 3.18E-13           | 2.78564637050161  | 8.85E-14          | 288    | gi|906468181|gb|KNC29231.1| hypothetical protein FF38_02464                                             | 71         | 1.66E-15           | 90.791786   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c6562_g1_i1   | 2.96796166735311  | 0.000459625014911 | 3.19740059008224  | 5.08E-05          | 1.81530986215646  | 0.002591382706644  | 1.94049315375036  | 0.00086119653709  | 1762   | gi|533135|gb|AAA21258.1| reverse transcriptase                                                          | 497        | 2.84E-159          | 554.329703  | GO:0090305 nucleic acid phosphodiester bond hydrolysis | GO:0006278 RNA-dependent DNA replication                                                                                                                                                                                                             | -                                                                                                             | GO:0046872 metal ion binding | GO:0004519 endonuclease activity | GO:0003964 RNA-directed DNA polymerase activity | GO:0003723 RNA binding                                                                    | -                     | 7216|Drosophila ambigua       | Drosophila | Drosophilidae |
| c6832_g1_i1   | -1.51299067514142 | 0.025852478879473 | -1.5172219562849  | 0.026480237666774 | -1.27075977516493 | 0.003148151182935  | -1.28118324318826 | 0.00254869911558  | 925    | gi|906468195|gb|KNC29245.1| hypothetical protein FF38_02361                                             | 220        | 1.13E-137          | 482.526771  | GO:0006749 glutathione metabolic process                                                                                                                                                                                                                                                                      | GO:0005737 cytoplasm                                                                                          | GO:0004364 glutathione transferase activity                                                                                                                                                                   | 2.5.1.18              | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c6842_g2_i1   | -1.35848824485306 | 0.046952763915074 | -1.44238776815709 | 0.043315549534431 | -1.1729529406158  | 0.002894538825979  | -1.20579732165168 | 0.009591687898781 | 435    | gi|906462855|gb|KNC25052.1| hypothetical protein FF38_11870                                             | 48         | 9.09E-22           | 112.150886  | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | -                                                                                                             | GO:0004364 glutathione transferase activity | GO:0018833 DDT-dehydrochlorinase activity                                                                                                                       | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c7227_g2_i1   | -2.8390104973113  | 0.000210035748827 | -2.84158401879042 | 0.000127333003922 | -1.75406057253234 | 0.004302634758796  | -1.79438887840135 | 0.002762939116429 | 325    | gi|906462046|gb|KNC24420.1| Larval cuticle protein 8, partial                                           | 92         | 6.4E-41            | 170.774799  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0042302 structural constituent of cuticle                                                                                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c7673_g1_i2   | -1.51175428057856 | 0.014869830243938 | -1.70699115022497 | 0.00156060162246  | -1.21779902387132 | 0.020533938201525  | -1.38022544332565 | 0.004611781408965 | 322    | gi|906462838|gb|KNC25035.1| hypothetical protein FF38_11920                                             | 49         | 4.48E-23           | 115.332029  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c7673_g2_i1   | -1.7520517158561  | 0.001494341326906 | -1.65682226950627 | 0.000541474694804 | -1.45157009333782 | 0.000868037991762  | -1.41542931523777 | 0.000374234044544 | 349    | gi|906462838|gb|KNC25035.1| hypothetical protein FF38_11920                                             | 51         | 7.83E-23           | 114.87758   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c7764_g1_i1   | 3.17685155390929  | 1.03E-05          | 3.07923230340451  | 4.96E-05          | 2.27848972076674  | 2.86E-06           | 2.09945370786116  | 5.35E-05          | 325    | gi|906468181|gb|KNC29231.1| hypothetical protein FF38_02464                                             | 81         | 4.38E-17           | 96.245174   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c7834_g1_i3   | 1.4853813964908   | 0.023074413058126 | 1.46257346704369  | 0.019044907932885 | 1.21409739164832  | 0.030676220961738  | 1.19878074678404  | 0.036147167657555 | 1110   | gi|906467433|gb|KNC28665.1| hypothetical protein FF38_04887                                             | 287        | 3.58E-155          | 540.696235  | GO:0030194 positive regulation of blood coagulation | GO:0051788 response to misfolded protein | GO:0051919 positive regulation of fibrinolysis | GO:0016540 protein autoprocessing | GO:0006508 proteolysis | GO:0002542 Factor XII activation | GO:0010756 positive regulation of plasminogen activation    | GO:0005615 extracellular space                                                                                | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                 | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8201_g1_i1   | 1.87723250291759  | 0.000525199698621 | 1.88857593615074  | 0.000249212101359 | 1.71326946726124  | 8.98E-08           | 1.71584600070471  | 1.07E-07          | 506    | gi|2565392|gb|AAB81989.1| cuticle 1                                                                     | 118        | 2.52E-52           | 208.494061  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0042302 structural constituent of cuticle                                                                                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8280_g1_i1   | 6.35442607807669  | 1.3E-14           | 6.31128957568657  | 2.9E-15           | 3.81038286240355  | 1.38E-16           | 3.89010053615415  | 5.55E-18          | 611    | gi|906465099|gb|KNC26810.1| hypothetical protein FF38_06382                                             | 114        | 2.04E-58           | 227.580917  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8591_g1_i1   | 2.0772161147496   | 0.02925979729101  | 2.12302294264834  | 0.017349670191204 | 1.4022544623534   | 0.03738511370928   | 1.44095650733462  | 0.029591309957906 | 366    | gi|906463584|gb|KNC25624.1| hypothetical protein FF38_03381                                             | 91         | 1.6E-32            | 145.780108  | GO:0006030 chitin metabolic process                                                                                                                                                                                                                                                                           | GO:0005576 extracellular region                                                                               | GO:0008061 chitin binding                                                                                                                                                                                     | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8614_g1_i1   | 6.0417462655695   | 9.94E-10          | 5.99527061665752  | 2.78E-10          | 2.87681981765833  | 3.17E-08           | 2.85681462517178  | 3.3E-08           | 546    | gi|906475050|gb|KNC34495.1| hypothetical protein FF38_02738                                             | 128        | 1.38E-63           | 242.577731  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8636_g1_i1   | 2.22617011309009  | 0.004697473095237 | 2.19973148164063  | 0.003560989870399 | 1.53081879802522  | 0.015598230613121  | 1.5330820301091   | 0.01493280272542  | 847    | gi|906468629|gb|KNC29604.1| hypothetical protein FF38_10438                                             | 137        | 5.92E-54           | 214.856346  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8735_g1_i1   | 1.3503375180011   | 0.022198678403738 | 1.31914311921078  | 0.014868239470259 | 1.19660447199998  | 0.006137383978896  | 1.16562599365136  | 0.009427628827031 | 1875   | gi|906472160|gb|KNC32342.1| hypothetical protein FF38_12474                                             | 323        | 1.14E-191          | 662.034101  | GO:0006812 cation transport | GO:0055085 transmembrane transport                                                                                                                                                                                                                                              | GO:0016021 integral to membrane                                                                               | GO:0008324 cation transmembrane transporter activity                                                                                                                                                          | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8766_g1_i1   | 1.87427033674376  | 0.000205906862184 | 1.87662743492253  | 6.62E-05          | 1.64941177826864  | 1.24E-05           | 1.65235767271233  | 1.01E-05          | 749    | gi|109290395|gb|ABG29412.1| small heat shock protein                                                    | 157        | 1.06E-92           | 333.01307   | GO:0006950 response to stress                                                                                                                                                                                                                                                                                 | GO:0005737 cytoplasm | GO:0005634 nucleus                                                                     | GO:0005212 structural constituent of eye lens | GO:0042802 identical protein binding                                                                                                                          | -                     | 7175|Culex pipiens            | Culex      | Culicidae     |
| c8857_g1_i1   | -1.74940514188675 | 0.002784468289927 | -1.74066512542913 | 0.002598272726976 | -1.34804631450641 | 0.014813409150143  | -1.34087606474582 | 0.017619093172464 | 666    | gi|907715916|ref|XP_013112074.1| PREDICTED: acyl-coenzyme A thioesterase 13-like isoform X1             | 145        | 3.33E-70           | 263.027934  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 35570|Stomoxys calcitrans     | Stomoxys   | Muscidae      |
| c8899_g1_i1   | -1.82192649945873 | 0.006870712691504 | -1.8150016179007  | 0.005320134854021 | -1.45022107236722 | 0.002538453598933  | -1.45288903005043 | 0.002268137927426 | 400    | gi|906473940|gb|KNC33642.1| Larval cuticle protein 1                                                    | 119        | 1.32E-48           | 195.769491  | GO:0008363 larval chitin-based cuticle development                                                                                                                                                                                                                                                            | GO:0031012 extracellular matrix                                                                               | GO:0008010 structural constituent of chitin-based larval cuticle                                                                                                                                              | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8899_g2_i1   | -1.90118271835504 | 0.001087884559591 | -1.89934775161815 | 0.000851967661459 | -1.54191453154206 | 0.000526252439424  | -1.54504555583388 | 0.000505066108555 | 288    | gi|906473940|gb|KNC33642.1| Larval cuticle protein 1                                                    | 41         | 8.4E-14            | 85.338399   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | GO:0042302 structural constituent of cuticle                                                                                                                                                                  | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c8918_g1_i1   | 4.36885081512889  | 2.29E-14          | 4.52781283939292  | 3.53E-15          | 3.50962515142693  | 8E-20              | 3.62724104161782  | 1.79E-21          | 377    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c9018_g1_i1   | 3.10614673640719  | 0.000332596359912 | 3.5255929153142   | 2.92E-05          | 1.87523607970049  | 0.001733544401375  | 1.96891049963836  | 0.000916443090001 | 1065   | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c9041_g1_i2   | 3.99199110628237  | 1.97E-09          | 4.07817742542486  | 2.66E-10          | 2.95462175101834  | 2.13E-11           | 3.05251548304725  | 1.31E-12          | 553    | gi|659105647|gb|AID61305.1| odorant binding protein                                                     | 158        | 1.15E-96           | 346.192089  | GO:0007606 sensory perception of chemical stimulus | GO:0035071 salivary gland cell autophagic cell death | GO:0019236 response to pheromone | GO:0042048 olfactory behavior                                                                                                                                  | GO:0005576 extracellular region                                                                               | GO:0005549 odorant binding                                                                                                                                                                                    | -                     | 145453|Calliphora stygia      | Calliphora | Calliphoridae |
| c9103_g1_i1   | 1.40023561652016  | 0.014841583509689 | 1.37544694580559  | 0.009463333005009 | 1.256201152407    | 0.001885598059606  | 1.23026080712179  | 0.002621524099812 | 1742   | gi|906466238|gb|KNC27710.1| H/ACA ribonucleoprotein complex subunit 4                                   | 425        | 4.29E-290          | 989.237336  | GO:0006278 RNA-dependent DNA replication | GO:0060216 definitive hemopoiesis | GO:0006364 rRNA processing | GO:0001522 pseudouridine synthesis                                                                                                                                                                | GO:0019013 viral nucleocapsid | GO:0005730 nucleolus | GO:0005697 telomerase holoenzyme complex               | GO:0003723 RNA binding | GO:0003720 telomerase activity | GO:0009982 pseudouridine synthase activity                                                                                                          | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9212_g1_i1   | 3.52256084661304  | 9.76E-08          | 3.67608154201076  | 3.06E-08          | 2.66020444029458  | 2.13E-09           | 2.74654411817928  | 5.55E-10          | 311    | gi|906461795|gb|KNC24208.1| hypothetical protein FF38_02191                                             | 102        | 1.65E-35           | 153.960189  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9248_g1_i2   | 3.01638065020885  | 0.000659735776064 | 3.00491583223076  | 0.00028988607602  | 1.8964149665371   | 0.001041082531185  | 1.92235578123618  | 0.000704006631805 | 524    | gi|906468148|gb|KNC29206.1| hypothetical protein FF38_01411                                             | 88         | 4.48E-20           | 107.151948  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9327_g1_i1   | 1.42849172658691  | 0.014967272437411 | 1.69355653792798  | 0.000643518914914 | 1.24883555126009  | 0.004868054228332  | 1.46018975140434  | 0.000587536515603 | 452    | gi|906468181|gb|KNC29231.1| hypothetical protein FF38_02464                                             | 84         | 1.14E-16           | 95.790725   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9354_g1_i1   | 1.45811152214025  | 0.03044205084054  | 1.43526648485443  | 0.020250729187056 | 1.23471211588339  | 0.012478790808125  | 1.2182934807609   | 0.014398930717709 | 6334   | gi|906464265|gb|KNC26122.1| hypothetical protein FF38_12214                                             | 2073       | 0                  | 3698.207455 | GO:0006200 ATP catabolic process | GO:0006461 protein complex assembly                                                                                                                                                                                                                                        | GO:0005634 nucleus                                                                                            | GO:0005524 ATP binding | GO:0016887 ATPase activity                                                                                                                                                           | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9502_g1_i1   | 2.92403140808503  | 0.000175462150222 | 2.95940676509344  | 1.66E-06          | 2.07669367265813  | 4.55E-05           | 2.24601005475721  | 8.05E-07          | 349    | gi|906468181|gb|KNC29231.1| hypothetical protein FF38_02464                                             | 85         | 1.03E-18           | 101.698561  | GO:0016560 protein import into peroxisome matrix, docking                                                                                                                                                                                                                                                     | GO:0005777 peroxisome | GO:0016021 integral to membrane                                                       | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9537_g1_i1   | 1.49064223817954  | 0.017903925707647 | 1.48793658638696  | 0.007806428345692 | 1.27860100160208  | 0.006515730242638  | 1.27414167842585  | 0.006681827052287 | 5184   | gi|906462690|gb|KNC24906.1| DNA-directed RNA polymerase I subunit RPA1                                  | 1664       | 0                  | 3508.247799 | GO:0009303 rRNA transcription | GO:0006362 transcription elongation from RNA polymerase I promoter | GO:0006363 termination of RNA polymerase I transcription | GO:0006361 transcription initiation from RNA polymerase I promoter                                                                            | GO:0005736 DNA-directed RNA polymerase I complex | GO:0005654 nucleoplasm                                     | GO:0008270 zinc ion binding | GO:0003677 DNA binding | GO:0005515 protein binding | GO:0001054 RNA polymerase I activity                                                                                      | 2.7.7.6               | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c954_g1_i1    | 1.73485567158455  | 0.003370870230529 | 1.72191237018089  | 0.001974109063469 | 1.54239473248218  | 5.78E-05           | 1.53531034230663  | 4.46E-05          | 876    | gi|557753623|ref|XP_005176556.1| PREDICTED: uncharacterized protein LOC101892259                        | 223        | 3.72E-111          | 394.363677  | GO:0040003 chitin-based cuticle development | GO:0010171 body morphogenesis                                                                                                                                                                                                                                   | GO:0005578 proteinaceous extracellular matrix                                                                 | GO:0005214 structural constituent of chitin-based cuticle                                                                                                                                                     | -                     | 7370|Musca domestica          | Musca      | Muscidae      |
| c9552_g1_i1   | -1.60946417797674 | 0.004530066603431 | -1.47038777742399 | 0.007874156648516 | -1.30798965768129 | 0.006515730242638  | -1.19191841242891 | 0.023587414674748 | 427    | gi|906467977|gb|KNC29079.1| hypothetical protein FF38_05521                                             | 75         | 1.35E-31           | 143.507863  | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9595_g1_i1   | -2.5988385000056  | 1.66E-05          | -2.63726154137244 | 4.99E-06          | -2.04082549904054 | 2.74E-06           | -2.11529740661519 | 2.84E-07          | 441    | gi|906469243|gb|KNC30057.1| hypothetical protein FF38_09551, partial                                    | 88         | 6.44E-12           | 80.339461   | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9754_g1_i1   | 1.81477957744105  | 0.026103883694706 | 1.73958069180336  | 0.028395027823346 | 1.42442808163816  | 0.008403117989443  | 1.38657664098199  | 0.009427628827031 | 304    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c9771_g1_i2   | 3.32516626053949  | 2.3E-08           | 3.55687650350627  | 2.13E-09          | 2.40393325757442  | 7.54E-07           | 2.57702807410149  | 4.23E-08          | 346    | -                                                                                                       | -          | -                  | -           | -                                                                                                                                                                                                                                                                                                             | -                                                                                                             | -                                                                                                                                                                                                             | -                     | -                             | -          | -             |
| c9883_g1_i1   | -2.86774445017259 | 0.000439916104495 | -2.8501972729914  | 0.000391295120582 | -1.9935489756585  | 0.000129913156329  | -2.02639462529601 | 5.73E-05          | 2228   | gi|906470722|gb|KNC31188.1| hypothetical protein FF38_14394                                             | 596        | 0                  | 1279.630208 | GO:0008152 metabolic process                                                                                                                                                                                                                                                                                  | GO:0005739 mitochondrion                                                                                      | GO:0016874 ligase activity                                                                                                                                                                                    | -                     | 7375|Lucilia cuprina          | Lucilia    | Calliphoridae |
| c9928_g1_i1   | 2.09036872347435  | 0.028463006149326 | 2.58210639978114  | 0.003836193684677 | 1.39667042550317  | 0.042436661108639  | 1.65186991186667  | 0.008374329898094 | 2204   | gi|1850812|emb|CAA55706.1| arylphorin receptor                                                          | 688        | 4.08E-169          | 587.050026  | GO:0015032 storage protein import into fat body                                                                                                                                                                                                                                                               | GO:0005811 lipid particle                                                                                     | GO:0008565 protein transporter activity                                                                                                                                                                       | -                     | 7373|Calliphora vicina        | Calliphora | Calliphoridae |



To get DEcontigs-UP-annotated.txt and DEcontigs-DOWN-annotated.txt, I just sorted by the fold change and removed the lines that corresponded to the DE contigs I didn’t want (the DOWN-reg in the DEcontigs-UP-annotated.txt file and vice versa).


**9.4 GO analysis of the DE transcripts (17/09/2019)**
****
To find if specific GO terms are found enriched in the up- and down-regulated genes, I conducted a GO analysis based on the GO annotation of the transcripts that was made by BLAST2GO (from within FunctionAnnotator). Because the file was outputted in a B2G format, I used B2G to do the GO analysis.

For this analysis the reference set was the reduced Trinity transcriptome (after clustering and taking the longest isoform) and tests set were the lists of up-regulated, downregulated and DE transcript (up and down-regulated transcripts combined together)

Fisher’s Exact test was applied to test for enrichment of GO terms in the test set versus the reference set. The list of GO terms found enriched was then reduced to the most specific ones.

The results are the following:


- Ref vs down-regulated (Reduced to specific terms)
| Tags   | GO ID      | GO Name                                | GO Category        | FDR                  | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | -------------------------------------- | ------------------ | -------------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0055114 | oxidation-reduction process            | BIOLOGICAL_PROCESS | 1.40519536817628E-06 | 4.48156711266553E-10 | 14      | 791          | 51             | 32187               |
| [OVER] | GO:0005506 | iron ion binding                       | MOLECULAR_FUNCTION | 1.07075172790063E-05 | 4.26866419989089E-09 | 8       | 184          | 57             | 32794               |
| [OVER] | GO:0071391 | cellular response to estrogen stimulus | BIOLOGICAL_PROCESS | 0.000201629633391    | 1.44687187092762E-07 | 3       | 3            | 62             | 32975               |
| [OVER] | GO:0004505 | phenylalanine 4-monooxygenase activity | MOLECULAR_FUNCTION | 0.001082685422668    | 8.63247825440935E-07 | 3       | 7            | 62             | 32971               |
| [OVER] | GO:0006559 | L-phenylalanine catabolic process      | BIOLOGICAL_PROCESS | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:0042302 | structural constituent of cuticle      | MOLECULAR_FUNCTION | 0.007999824643597    | 1.08433279334354E-05 | 5       | 142          | 60             | 32836               |
| [OVER] | GO:0016597 | amino acid binding                     | MOLECULAR_FUNCTION | 0.014409111666401    | 2.2977374687292E-05  | 3       | 25           | 62             | 32953               |
| [OVER] | GO:0020037 | heme binding                           | MOLECULAR_FUNCTION | 0.021343464463242    | 3.91404626578342E-05 | 5       | 187          | 60             | 32791               |
| [OVER] | GO:0046394 | carboxylic acid biosynthetic process   | BIOLOGICAL_PROCESS | 0.043514499760038    | 9.43769796938872E-05 | 4       | 116          | 61             | 32862               |




- Ref versus up-regulated (Reduced to specific terms) 
| Tags   | GO ID      | GO Name                                     | GO Category        | FDR               | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | ------------------------------------------- | ------------------ | ----------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0004252 | serine-type endopeptidase activity          | MOLECULAR_FUNCTION | 0.000221563791679 | 4.18371091872965E-08 | 8       | 208          | 69             | 32758               |
| [OVER] | GO:0005576 | extracellular region                        | CELLULAR_COMPONENT | 0.000221563791679 | 2.22666144327842E-08 | 14      | 897          | 63             | 32069               |
| [OVER] | GO:0010765 | positive regulation of sodium ion transport | BIOLOGICAL_PROCESS | 0.006111286322295 | 5.35992262360444E-06 | 2       | 0            | 75             | 32966               |
| [OVER] | GO:0055085 | transmembrane transport                     | BIOLOGICAL_PROCESS | 0.008084182577661 | 8.37939511318746E-06 | 10      | 729          | 67             | 32237               |
| [OVER] | GO:0006508 | proteolysis                                 | BIOLOGICAL_PROCESS | 0.012585453969582 | 1.50287899956786E-05 | 11      | 959          | 66             | 32007               |
| [OVER] | GO:0001518 | voltage-gated sodium channel complex        | CELLULAR_COMPONENT | 0.012585453969582 | 1.6055434820069E-05  | 2       | 1            | 75             | 32965               |
| [OVER] | GO:0006950 | response to stress                          | BIOLOGICAL_PROCESS | 0.020753980308028 | 2.97856518533323E-05 | 11      | 1034         | 66             | 31932               |
| [OVER] | GO:0005248 | voltage-gated sodium channel activity       | MOLECULAR_FUNCTION | 0.021164483263117 | 3.20622852813926E-05 | 2       | 2            | 75             | 32964               |




- Ref versus all DE transcripts (Reduced to specific terms) 


I also saved the unreduced list of enriched GO terms:


- Ref versus down (complete result)
| Tags   | GO ID      | GO Name                                                                                                                                                                        | GO Category        | FDR                  | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------ | -------------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0004497 | monooxygenase activity                                                                                                                                                         | MOLECULAR_FUNCTION | 5.0824287840599E-09  | 6.27768719656299E-13 | 10      | 147          | 55             | 32831               |
| [OVER] | GO:0003824 | catalytic activity                                                                                                                                                             | MOLECULAR_FUNCTION | 5.0824287840599E-09  | 8.10465441565922E-13 | 35      | 5056         | 30             | 27922               |
| [OVER] | GO:0016491 | oxidoreductase activity                                                                                                                                                        | MOLECULAR_FUNCTION | 8.88460756202465E-07 | 2.12516525961361E-10 | 14      | 746          | 51             | 32232               |
| [OVER] | GO:0055114 | oxidation-reduction process                                                                                                                                                    | BIOLOGICAL_PROCESS | 1.40519536817628E-06 | 4.48156711266553E-10 | 14      | 791          | 51             | 32187               |
| [OVER] | GO:0005506 | iron ion binding                                                                                                                                                               | MOLECULAR_FUNCTION | 1.07075172790063E-05 | 4.26866419989089E-09 | 8       | 184          | 57             | 32794               |
| [OVER] | GO:0016705 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen                                                                          | MOLECULAR_FUNCTION | 1.49132082905145E-05 | 7.13436850128263E-09 | 8       | 197          | 57             | 32781               |
| [OVER] | GO:0008152 | metabolic process                                                                                                                                                              | BIOLOGICAL_PROCESS | 3.05808866140369E-05 | 1.70679481979157E-08 | 34      | 6797         | 31             | 26181               |
| [OVER] | GO:0043627 | response to estrogen                                                                                                                                                           | BIOLOGICAL_PROCESS | 0.000201629633391    | 1.44687187092762E-07 | 3       | 3            | 62             | 32975               |
| [OVER] | GO:0071391 | cellular response to estrogen stimulus                                                                                                                                         | BIOLOGICAL_PROCESS | 0.000201629633391    | 1.44687187092762E-07 | 3       | 3            | 62             | 32975               |
| [OVER] | GO:0004505 | phenylalanine 4-monooxygenase activity                                                                                                                                         | MOLECULAR_FUNCTION | 0.001082685422668    | 8.63247825440935E-07 | 3       | 7            | 62             | 32971               |
| [OVER] | GO:0006558 | L-phenylalanine metabolic process                                                                                                                                              | BIOLOGICAL_PROCESS | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:0006559 | L-phenylalanine catabolic process                                                                                                                                              | BIOLOGICAL_PROCESS | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:0016714 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced pteridine as one donor, and incorporation of one atom of oxygen | MOLECULAR_FUNCTION | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:1902222 | erythrose 4-phosphate/phosphoenolpyruvate family amino acid catabolic process                                                                                                  | BIOLOGICAL_PROCESS | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:1902221 | erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process                                                                                                  | BIOLOGICAL_PROCESS | 0.002177138323796    | 2.60381716288817E-06 | 3       | 11           | 62             | 32967               |
| [OVER] | GO:0009074 | aromatic amino acid family catabolic process                                                                                                                                   | BIOLOGICAL_PROCESS | 0.004549897529409    | 5.80436616732131E-06 | 3       | 15           | 62             | 32963               |
| [OVER] | GO:0042302 | structural constituent of cuticle                                                                                                                                              | MOLECULAR_FUNCTION | 0.007999824643597    | 1.08433279334354E-05 | 5       | 142          | 60             | 32836               |
| [OVER] | GO:0009072 | aromatic amino acid family metabolic process                                                                                                                                   | BIOLOGICAL_PROCESS | 0.009947273218449    | 1.42761057193502E-05 | 3       | 21           | 62             | 32957               |
| [OVER] | GO:0044281 | small molecule metabolic process                                                                                                                                               | BIOLOGICAL_PROCESS | 0.012587910321031    | 1.90695499999665E-05 | 9       | 760          | 56             | 32218               |
| [OVER] | GO:0016597 | amino acid binding                                                                                                                                                             | MOLECULAR_FUNCTION | 0.014409111666401    | 2.2977374687292E-05  | 3       | 25           | 62             | 32953               |
| [OVER] | GO:0042737 | drug catabolic process                                                                                                                                                         | BIOLOGICAL_PROCESS | 0.019139009498937    | 3.20458618623564E-05 | 4       | 87           | 61             | 32891               |
| [OVER] | GO:0046914 | transition metal ion binding                                                                                                                                                   | MOLECULAR_FUNCTION | 0.02011915850333     | 3.52911407329984E-05 | 12      | 1505         | 53             | 31473               |
| [OVER] | GO:0020037 | heme binding                                                                                                                                                                   | MOLECULAR_FUNCTION | 0.021343464463242    | 3.91404626578342E-05 | 5       | 187          | 60             | 32791               |
| [OVER] | GO:0046906 | tetrapyrrole binding                                                                                                                                                           | MOLECULAR_FUNCTION | 0.022563320937467    | 4.31765031493545E-05 | 5       | 191          | 60             | 32787               |
| [OVER] | GO:1901605 | alpha-amino acid metabolic process                                                                                                                                             | BIOLOGICAL_PROCESS | 0.023276487705212    | 4.63970812175329E-05 | 4       | 96           | 61             | 32882               |
| [OVER] | GO:1901606 | alpha-amino acid catabolic process                                                                                                                                             | BIOLOGICAL_PROCESS | 0.040883328000975    | 8.47525536617235E-05 | 3       | 40           | 62             | 32938               |
| [OVER] | GO:0016053 | organic acid biosynthetic process                                                                                                                                              | BIOLOGICAL_PROCESS | 0.043514499760038    | 0.000100615571124    | 4       | 118          | 61             | 32860               |
| [OVER] | GO:0005198 | structural molecule activity                                                                                                                                                   | MOLECULAR_FUNCTION | 0.043514499760038    | 9.95384339049966E-05 | 7       | 540          | 58             | 32438               |
| [OVER] | GO:0046394 | carboxylic acid biosynthetic process                                                                                                                                           | BIOLOGICAL_PROCESS | 0.043514499760038    | 9.43769796938872E-05 | 4       | 116          | 61             | 32862               |
| [OVER] | GO:0019752 | carboxylic acid metabolic process                                                                                                                                              | BIOLOGICAL_PROCESS | 0.047209389979785    | 0.000112923114287    | 6       | 380          | 59             | 32598               |
| [OVER] | GO:0043436 | oxoacid metabolic process                                                                                                                                                      | BIOLOGICAL_PROCESS | 0.047329032430755    | 0.00012281642485     | 6       | 386          | 59             | 32592               |
| [OVER] | GO:0009063 | cellular amino acid catabolic process                                                                                                                                          | BIOLOGICAL_PROCESS | 0.047329032430755    | 0.000117950269508    | 3       | 45           | 62             | 32933               |
| [OVER] | GO:0006082 | organic acid metabolic process                                                                                                                                                 | BIOLOGICAL_PROCESS | 0.047329032430755    | 0.000124530224064    | 6       | 387          | 59             | 32591               |




- Ref versus up (complete result)
| Tags   | GO ID      | GO Name                                                      | GO Category        | FDR               | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | ------------------------------------------------------------ | ------------------ | ----------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0017171 | serine hydrolase activity                                    | MOLECULAR_FUNCTION | 0.000221563791679 | 8.83287321315109E-08 | 8       | 230          | 69             | 32736               |
| [OVER] | GO:0005576 | extracellular region                                         | CELLULAR_COMPONENT | 0.000221563791679 | 2.22666144327842E-08 | 14      | 897          | 63             | 32069               |
| [OVER] | GO:0004252 | serine-type endopeptidase activity                           | MOLECULAR_FUNCTION | 0.000221563791679 | 4.18371091872965E-08 | 8       | 208          | 69             | 32758               |
| [OVER] | GO:0016825 | hydrolase activity, acting on acid phosphorus-nitrogen bonds | MOLECULAR_FUNCTION | 0.000221563791679 | 8.83287321315109E-08 | 8       | 230          | 69             | 32736               |
| [OVER] | GO:0008236 | serine-type peptidase activity                               | MOLECULAR_FUNCTION | 0.000221563791679 | 8.55185812393667E-08 | 8       | 229          | 69             | 32737               |
| [OVER] | GO:0016787 | hydrolase activity                                           | MOLECULAR_FUNCTION | 0.000541539854701 | 2.59068659560511E-07 | 20      | 2314         | 57             | 30652               |
| [OVER] | GO:0043269 | regulation of ion transport                                  | BIOLOGICAL_PROCESS | 0.000738914574209 | 4.12406475798355E-07 | 5       | 59           | 72             | 32907               |
| [OVER] | GO:0004175 | endopeptidase activity                                       | MOLECULAR_FUNCTION | 0.0022511406162   | 1.43590535238375E-06 | 9       | 458          | 68             | 32508               |
| [OVER] | GO:0070011 | peptidase activity, acting on L-amino acid peptides          | MOLECULAR_FUNCTION | 0.005081810505997 | 3.64665081757096E-06 | 10      | 662          | 67             | 32304               |
| [OVER] | GO:0006814 | sodium ion transport                                         | BIOLOGICAL_PROCESS | 0.005649761275657 | 4.50467331817691E-06 | 4       | 43           | 73             | 32923               |
| [OVER] | GO:0010765 | positive regulation of sodium ion transport                  | BIOLOGICAL_PROCESS | 0.006111286322295 | 5.35992262360444E-06 | 2       | 0            | 75             | 32966               |
| [OVER] | GO:0008233 | peptidase activity                                           | MOLECULAR_FUNCTION | 0.008084182577661 | 8.08835955975312E-06 | 10      | 726          | 67             | 32240               |
| [OVER] | GO:0055085 | transmembrane transport                                      | BIOLOGICAL_PROCESS | 0.008084182577661 | 8.37939511318746E-06 | 10      | 729          | 67             | 32237               |
| [OVER] | GO:0006508 | proteolysis                                                  | BIOLOGICAL_PROCESS | 0.012585453969582 | 1.50287899956786E-05 | 11      | 959          | 66             | 32007               |
| [OVER] | GO:0001518 | voltage-gated sodium channel complex                         | CELLULAR_COMPONENT | 0.012585453969582 | 1.6055434820069E-05  | 2       | 1            | 75             | 32965               |
| [OVER] | GO:0034706 | sodium channel complex                                       | CELLULAR_COMPONENT | 0.012585453969582 | 1.6055434820069E-05  | 2       | 1            | 75             | 32965               |
| [OVER] | GO:0006812 | cation transport                                             | BIOLOGICAL_PROCESS | 0.020753980308028 | 2.83498973181307E-05 | 7       | 369          | 70             | 32597               |
| [OVER] | GO:0006950 | response to stress                                           | BIOLOGICAL_PROCESS | 0.020753980308028 | 2.97856518533323E-05 | 11      | 1034         | 66             | 31932               |
| [OVER] | GO:0005248 | voltage-gated sodium channel activity                        | MOLECULAR_FUNCTION | 0.021164483263117 | 3.20622852813926E-05 | 2       | 2            | 75             | 32964               |
| [OVER] | GO:0006811 | ion transport                                                | BIOLOGICAL_PROCESS | 0.024028236327629 | 3.83164349029328E-05 | 8       | 533          | 69             | 32433               |
| [OVER] | GO:0002028 | regulation of sodium ion transport                           | BIOLOGICAL_PROCESS | 0.030522708951724 | 5.33563042054051E-05 | 2       | 3            | 75             | 32963               |
| [OVER] | GO:0003824 | catalytic activity                                           | MOLECULAR_FUNCTION | 0.030522708951724 | 5.35400731093858E-05 | 26      | 5065         | 51             | 27901               |




- Ref versus DE (complete result)

**GO graphs (based on complete set of GO terms)**


![GO graph representing the Biological Processes terms found enriched in the transcripts up-regulated in the resistant condition](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568812307676_UP-BP-unreduced.png)
![GO graph representing the Cellular Component terms found enriched in the transcripts up-regulated in the resistant condition](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568812308793_UP-CC-unreduced.png)
![GO graph representing the Molecular Processes terms found enriched in the transcripts up-regulated in the resistant condition](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568812639773_UP-MF-unreduced.png)





![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568812817504_DOWN-BP-unreduced.png)
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568812816488_DOWN-MF-unreduced.png)



Barcharts of the GO terms found enriched in the up- and down-reguated transcripts in the resistant condition

![Over-represented GO terms in the transcripts down-regulated  in the resistant condition](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568815529532_OverRep-GO-DOWN-Barchart.png)
![Over-represented GO terms in the transcripts up-regulated  in the resistant condition](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568815529740_OverRep-GO-UP-Barchart.png)

## 10 Annotation 

10.1 Annotation for 142 differentially expressed contains @Pedro M 


[ ] @Gisele A 

10.1.1 @Gisele A Gisele como foi que você chegou na tabela com 81 genes a partir da de 142?

10.1.2 AgBase-GOanna @Gisele A preciso da tabela com os 81 genes que eu joguei nos bancos de dados do AgBase

Link to the databases: http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi

AgBase Community
**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 40
**Query coverage filter (%):** 40
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all

No results

**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/AgBase\ Community/aqvhud1503494149\ gene\ 30%/GOanna2GA_Reformat_ytf6co1503494430.xls


**Program:** blastp
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/AgBase\ Community/aqvhud1503494149\ protein\ 30%/GOanna2GA_Reformat_vefcwd1503496388.xls

Invertebrates
**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 40
**Query coverage filter (%):** 40
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all

No results

**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/Invertebrates/537f9x1503494292\ gene\ 30%/GOanna2GA_Reformat_gvozi91503495336.xls

**Program:** blastp
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/Invertebrates/8ci7wl1503496983\ protein\ 30%/GOanna2GA_Reformat_twi6mt1503497052.xls

Uniprot
**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 40
**Query coverage filter (%):** 40
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/UniProt/twe68z1503492944\ gene\ 40%/GOanna2GA_Reformat_dz9k661503495014.xls

**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/UniProt/vphrhi1503495637\ gene\ 30%/GOanna2GA_Reformat_dlfl1x1503497237.xls

**Program:** blastp
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/UniProt/vgk25a1503496499\ protein\ 30%/GOanna2GA_Reformat_141cx91503496612.xls

SwissProt
**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 40
**Query coverage filter (%):** 40
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all

No results

**Program:** blastx
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/SwissProt/99qv451503494192\ gene\ 30%/GOanna2GA_Reformat_5k7jt11503495149.xls

**Program:** blastp
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 1 (default)
**Pct identity filter (%):** 30
**Query coverage filter (%):** 30
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all


    cd /mnt/HD5/PedroMartins/GOanna/SwissProt/wb8dmt1503496656\ protein\ 30%/GOanna2GA_Reformat_k3623m1503496890.xls

All the resulting tables were merged, and the final results are shown below


    cd /mnt/HD5/PedroMartins/final8_conferido_GO.txt

10.1.3 NCBI-BLAST manual
The merged table did not contain all the sequences that should be analyzed, so we ran a Nucleotide BLAST (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) on the sequences that were still left. But we only took the sequences with an id of 40% or higher.

**Database:** Others - Nucleotide collection (nr/nt)
**Optimize for:** More dissimilar sequences (discontiguous megablast)

We also tried to find the GO terms associated to each aligned sequence (http://www.geneontology.org).

The results are shown below.


    cd /mnt/HD5/PedroMartins/DE_still_missing_annot_blastn_nr_parser_GO-table.txt

10.2 Gisele
Gisele is annotating the contigs using Agbase and blast+.

10.2.1 Agbase - Goanna

McCarthy *et al*. **2006. [AgBase: a functional genomics resource for agriculture](http://www.biomedcentral.com/1471-2164/7/229). *BMC Genomics* 2006, 7:229
[http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi](http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi)
 
For proteins: limited to 30,0000 sequences (I used Transdecoder to get proteins because the transcriptome has more than 40,000 sequences).
For nucleotides: limited to 5,000 sequences.
 
**Program:** blastp
**Database:** AgBase-Uniprot
**e-value:** 10e-4
**Matrix:** BLOSUM62 (default)
**Word size:** 6 (default)
**Number of target seqs**: 5 (default)
**Pct identity filter (%):** 70
**Query coverage filter (%):** 40
**Gap costs:** Existence: 11 extension:1 (default)
**Blast resulta format selection**: both HTML and TSV format (default)
**Type of evidence to return:** check all
 
Link to the original output file: [Chom-resistance_agbase\Chom-resistance_agbase.xls](#)

Goanna2ga ([http://www.agbase.msstate.edu/cgi-bin/tools/GOanna2ga.cgi](http://www.agbase.msstate.edu/cgi-bin/tools/GOanna2ga.cgi))

GOanna2ga is used to convert the output from GOanna into a gene association file format. Use this to add reviewed data obtained from GOanna to the GORetriever annotations. **(website text)**
 
You just have to choose the files and hit submit!
 
Link to the original output files:
 
[Chom-resistance_agbase\Chom-resistance_agbase_GOanna2GA_GOSummary.txt](#)
[Chom-resistance_agbase\Chom-resistance_agbase_GOanna2GA_Reformat.xls](#)
 
With this method we found 7,807 proteins with a linked GO. It seems just a small number compared to the total of sequences (~40,000), but this annotation includes only genes that were functionally tested in previous studies.

**TransDecoder**

With the assembly, it’s possible to identify candidate coding regions within the transcript sequences using TransDecoder (current: v5.0.1)


        > Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q et al. 2011. Full-length transcriptome assembly from RNA- Seq data without a reference genome. Nat Biotechnol **29**: 644-652.
        > 
        > https://github.com/TransDecoder/TransDecoder/wiki → manual and link to download the Pfam-A database

The parameter -m represent the minimum length of protein to be estimated (considering that the minimum contig length is 200 bp)

Upload of the Pfam-A database to server and preparing the database; it’s necessary to have the Hmmer installed (version: 3.1b2)


        > http://hmmer.org
    
    mkdir /mnt/HD5/Chom-resistance/Anotacao/1-TransDecoder


    cd /mnt/HD5/Chom-resistance/Anotacao/1-TransDecoder


            Extracting the long open reading frames (ORFs)


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Trinity.Trinity.fasta


            Identifying ORFs with homology to known proteins via pfam searches


    hmmscan --cpu 10 --domtblout Trinity.pfam /mnt/HD3/Ilp8-evolution/DataBase/Pfam-A.hmm Trinity.Trinity.fasta.transdecoder_dir/longest_orfs.pep


            Predicting the likely coding regions


    TransDecoder.Predict -t /mnt/HD5/Chom-resistance/Trinity.Trinity.fasta --retain_pfam_hits Trinity.pfam

**SignalP**


    mkdir /mnt/HD5/Chom-resistance/Anotacao/2-SignalP


    cd /mnt/HD5/Chom-resistance/Anotacao/2-SignalP

With the peptide predicted by TransDecoder it’s possible to predict the presence and location of signal peptide cleavage sites (version 4.1)


        > Petersen TN, Brunak S, von Heijne G, Nielsen H. 2011. SignalP 4.0: discriminating signal peptides from transmembrane regions. Nature **8**:785-786


            SignalP doesn't accept files with more than 10,000 sequences, so first it’s necessary to split the file with peptide sequences


    perl /mnt/HD3/Ilp8-evolution/Scripts/split_fasta_file.pl /mnt/HD5/Chom-resistance/Anotacao/1-TransDecoder/Trinity.Trinity.fasta.transdecoder.pep 10000 Trinity


            It’s generated 2 files; running signalP using each file:


    signalp -n Trinity_signalp_1 Trinity_1.fasta


    signalp -n Trinity_signalp_2 Trinity_2.fasta



            Only files 1 and 5 had any sequence with signal peptide; putting them together  and removing temporary files


    tail -n +4 HillF_signalp_5 >HillF_signalp_5_tmp


    cat MsexL_signalp_2 MsexL_signalp_4_tmp MsexL_signalp_5_tmp >09-Annotation/02-SignalP/MsexL.signalP | rm MsexL_*

**TMHMM**

With the peptide predicted by TransDecoder it’s possible to predict transmembrane helices domains in the sequences (version 2.0)


        > Sonnhammer EL, von Heijne G, Krogh A. 1998. A hidden Markov model for predicting transmembrane helices in protein sequences. Proc Int Conf Intell Syst Mol Biol **6**:175-182.
        > 
        > Krogh A, Larsson B, von Heijne G, Sonnhammer ELL. 2001. Predicting transmembrane protein topology with a hidden Markov model: Application to complete genomes. Journal of Molecular Biology, **305(3)**:567-580


    mkdir /mnt/HD5/Chom-resistance/Anotacao/3-TMHMM


    cd /mnt/HD5/Chom-resistance/Anotacao/3-TMHMM


    tmhmm /mnt/HD5/Chom-resistance/Anotacao/1-TransDecoder/Trinity.Trinity.fasta.transdecoder.pep >/mnt/HD5/Chom-resistance/Anotacao/3-TMHMM/Trinity.tmhmm


    mv TMHMM_21948 /mnt/HD5/Chom-resistance/Anotacao/3-TMHMM/Trinity_graphics


            Getting only the sequences with transmembrane helices


    grep "Number of predicted TMHs:" /mnt/HD5/Chom-resistance/Anotacao/3-TMHMM/Trinity.tmhmm | awk '!/Number of predicted TMHs:  0/ {print}' >/mnt/HD5/Chom-resistance/Anotacao/3-TMHMM/Trinity_filtred.tmhmm

**GHOSTZ**

Alignment of transcriptome against the ncbi database nr (non-redundant protein sequence database with entries from GenPept, Swissprot, PIR, PDF, PDB, and RefSeq) using GHOSTZ (200 times more efficient than BLAST by using database subsequence clustering)


        > Suzuki S, Kakuta M, Ishida T, Akiyama Y. 2015. Faster sequence homology searches by clustering subsequences. Bioinformatics **31**:1183-1190.
        > 
        > http://www.bi.cs.titech.ac.jp/ghostz/

Url used to download the database: ftp://ftp.ncbi.nlm.nih.gov/blast/db/


    mkdir /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ
            
    cd /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ
            Aligning the transcriptomes


    ghostz aln -i /mnt/HD5/Chom-resistance/Trinity.Trinity.fasta -o /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.dna.ghst -d /mnt/HD3/Ilp8-evolution/DataBase/nr -b 1 -q d -a 10


    ghostz aln -i /mnt/HD5/Chom-resistance/Anotacao/1-TransDecoder/Trinity.Trinity.fasta.transdecoder.pep -o /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.pep.ghst -d /mnt/HD3/Ilp8-evolution/DataBase/nr -b 1 -q p -a 10


            Filtering by id (70%) and e-value (10e-04)


    awk -F "\t" 'NR==1; {if ($3 >= 70 && $11 <= 10e-5) print $0}' /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.dna.ghst >/mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.dna.ghst.filt


    awk -F "\t" 'NR==1; {if ($3 >= 70 && $11 <= 10e-5) print $0}' /mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.pep.ghst >/mnt/HD5/Chom-resistance/Anotacao/4-GHOSTZ/Trinity.pep.ghst.filt

**AgBase - GOanna**

GOanna allows users to quickly add more GO annotations to their data by transferring GO annotations based on sequence homology.


        > McCarthy FM, Wang N, Magee GB, Nanduri B, Lawrence ML, Camon EB, Barrell DG, Hill DP, Dolan ME, Williams WP, Luthe DS, Bridges SM, Burgess SC. 2006. AgBase: a functional genomics resource for agriculture. BMC Genomics, **7:229**
        > 
        > http://www.agbase.msstate.edu/cgi-bin/tools/GOanna.cgi


    mkdir /mnt/HD5/Chom-resistance/Anotacao/5-AgBase


    cd /mnt/HD5/Chom-resistance/Anotacao/5-AgBase

Since input files are limited to 30,000 protein sequences or 5,000 nucleotide sequences, I had to first split the fasta file


    perl /mnt/HD3/Ilp8-evolution/Scripts/split_fasta_file.pl /mnt/HD5/Chom-resistance/Trinity.Trinity.fasta 5000 /mnt/HD5/Chom-resistance/Anotacao/5-AgBase/1-Query/Trinity.dna


            → 9 files

There were less then 30000 protein sequences in our files, so we just copied the TransDecoder file in the Query directory  


    cp /mnt/HD5/Chom-resistance/Anotacao/TransDecoder/Trinity.Trinity.fasta.transdecoder.pep /mnt/HD5/Chom-resistance/Anotacao/5-AgBase/1-Query/Trinity_pep.fasta



All files of each transcriptome were downloaded locally using the command


    scp -r student@143.107.247.2:/mnt/HD5/Chom-resistance/Anotacao/5-AgBase/ Desktop

and uploaded in their website for alignment using the blastx (DNA query against proteins database) and blastp (protein query to protein database) programs. The AgBaseCommunity, Uniprot, SwissProt and Invertebrates database and default parameters were used, except for e-value (10e-4), Identity Filter (40%), Query Coverage Filter (40%),  Nbr. of Target Seqs (1) and type of Evidence to Return (check all). The results were sent by email, downloaded locally 

Since the DNA sequences had to be divided in 9 files, the resulting annotation tables obtained in each database were merged using the commands

    cd /Users/Pedro/Documents/2_Annotation/dna/AgBaseCommunity
    
    cat *_1_*/*.txt *_2_*/*.txt *_3_*/*.txt *_4_*/*.txt *_5_*/*.txt *_6_*/*.txt *_7_*/*.txt *_8_*/*.txt *_9_*/*.txt >Trinity_dna_AgBaseCommunity.txt
    
    tail -n +2 *_9_*/*.tsv >tmp9
    tail -n +2 *_8_*/*.tsv >tmp8
    tail -n +2 *_7_*/*.tsv >tmp7
    tail -n +2 *_6_*/*.tsv >tmp6
    tail -n +2 *_5_*/*.tsv >tmp5
    tail -n +2 *_4_*/*.tsv >tmp4
    tail -n +2 *_3_*/*.tsv >tmp3
    tail -n +2 *_2_*/*.tsv >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_AgBaseCommunity.tsv
    
    rm tmp*
    
    tail -n +5 *_9_*/*.xls >tmp9
    tail -n +5 *_8_*/*.xls >tmp8
    tail -n +5 *_7_*/*.xls >tmp7
    tail -n +5 *_6_*/*.xls >tmp6
    tail -n +5 *_5_*/*.xls >tmp5
    tail -n +5 *_4_*/*.xls >tmp4
    tail -n +5 *_3_*/*.xls >tmp3
    tail -n +5 *_2_*/*.xls >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_AgBaseCommunity.xls
    
    rm tmp*


    cd /Users/Pedro/Documents/2_Annotation/dna/Invertebrates
    
    cat *_1_*/*.txt *_2_*/*.txt *_3_*/*.txt *_4_*/*.txt *_5_*/*.txt *_6_*/*.txt *_7_*/*.txt *_8_*/*.txt *_9_*/*.txt >Trinity_dna_Invertebrates.txt
    
    tail -n +2 *_9_*/*.tsv >tmp9
    tail -n +2 *_8_*/*.tsv >tmp8
    tail -n +2 *_7_*/*.tsv >tmp7
    tail -n +2 *_6_*/*.tsv >tmp6
    tail -n +2 *_5_*/*.tsv >tmp5
    tail -n +2 *_4_*/*.tsv >tmp4
    tail -n +2 *_3_*/*.tsv >tmp3
    tail -n +2 *_2_*/*.tsv >tmp2
    
    cat *_1_*/*.xls tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_Invertebrates.tsv
    
    rm tmp*
    
    tail -n +5 *_9_*/*.xls >tmp9
    tail -n +5 *_8_*/*.xls >tmp8
    tail -n +5 *_7_*/*.xls >tmp7
    tail -n +5 *_6_*/*.xls >tmp6
    tail -n +5 *_5_*/*.xls >tmp5
    tail -n +5 *_4_*/*.xls >tmp4
    tail -n +5 *_3_*/*.xls >tmp3
    tail -n +5 *_2_*/*.xls >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_Invertebrates.xls
    
    rm tmp*


    cd /Users/Pedro/Documents/2_Annotation/dna/SwissProt
    
    cat *_1_*/*.txt *_2_*/*.txt *_3_*/*.txt *_4_*/*.txt *_5_*/*.txt *_6_*/*.txt *_7_*/*.txt *_8_*/*.txt *_9_*/*.txt >Trinity_dna_SwissProt.txt
    
    tail -n +2 *_9_*/*.tsv >tmp9
    tail -n +2 *_8_*/*.tsv >tmp8
    tail -n +2 *_7_*/*.tsv >tmp7
    tail -n +2 *_6_*/*.tsv >tmp6
    tail -n +2 *_5_*/*.tsv >tmp5
    tail -n +2 *_4_*/*.tsv >tmp4
    tail -n +2 *_3_*/*.tsv >tmp3
    tail -n +2 *_2_*/*.tsv >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_SwissProt.tsv
    
    rm tmp*
    
    tail -n +5 *_9_*/*.xls >tmp9
    tail -n +5 *_8_*/*.xls >tmp8
    tail -n +5 *_7_*/*.xls >tmp7
    tail -n +5 *_6_*/*.xls >tmp6
    tail -n +5 *_5_*/*.xls >tmp5
    tail -n +5 *_4_*/*.xls >tmp4
    tail -n +5 *_3_*/*.xls >tmp3
    tail -n +5 *_2_*/*.xls >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_SwissProt.xls
    
    rm tmp*


    cd /Users/Pedro/Documents/2_Annotation/dna/UniProt
    
    cat *_1_*/*.txt *_2_*/*.txt *_3_*/*.txt *_4_*/*.txt *_5_*/*.txt *_6_*/*.txt *_7_*/*.txt *_8_*/*.txt *_9_*/*.txt >Trinity_dna_UniProt.txt
    
    tail -n +2 *_9_*/*.tsv >tmp9
    tail -n +2 *_8_*/*.tsv >tmp8
    tail -n +2 *_7_*/*.tsv >tmp7
    tail -n +2 *_6_*/*.tsv >tmp6
    tail -n +2 *_5_*/*.tsv >tmp5
    tail -n +2 *_4_*/*.tsv >tmp4
    tail -n +2 *_3_*/*.tsv >tmp3
    tail -n +2 *_2_*/*.tsv >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_UniProt.tsv
    
    rm tmp*
    
    tail -n +5 *_9_*/*.xls >tmp9
    tail -n +5 *_8_*/*.xls >tmp8
    tail -n +5 *_7_*/*.xls >tmp7
    tail -n +5 *_6_*/*.xls >tmp6
    tail -n +5 *_5_*/*.xls >tmp5
    tail -n +5 *_4_*/*.xls >tmp4
    tail -n +5 *_3_*/*.xls >tmp3
    tail -n +5 *_2_*/*.xls >tmp2
    
    cat *_1_*/*.tsv tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 >Trinity_dna_UniProt.xls
    
    rm tmp*

Then, we ended up with 12 dna tables (three from each database, one .txt, one .tsv and one .xls)

The .xls tables were reformatted using the Goanna2ga program. GOanna2ga is used to convert the output from GOanna into a gene association file format
GOanna2ga: http://www.agbase.msstate.edu/cgi-bin/tools/GOanna2ga.cgi

Using R to merge the GOanna2ga .xls tables:


    ABC_dna <- read.delim("/Users/Pedro/Documents/2_Annotation/dna/AgBaseCommunity/dna_AgBaseCommunity_reformat.xls",header=TRUE)
    ABC_dna <- ABC_dna[c(2,5,6,10,11)]
    Inv_dna <- read.delim("/Users/Pedro/Documents/2_Annotation/dna/Invertebrates/dna_Invertebrates_reformat.xls",header=TRUE)
    Inv_dna <- Inv_dna[c(2,5,6,10,11)]
    SwP_dna <- read.delim("/Users/Pedro/Documents/2_Annotation/dna/SwissProt/dna_SwissProt_reformat.xls",header=TRUE)
    SwP_dna <- SwP_dna[c(2,5,6,10,11)]
    UnP_dna <- read.delim("/Users/Pedro/Documents/2_Annotation/dna/UniProt/dna_UniProt_reformat.xls",header=TRUE)
    UnP_dna <- UnP_dna[c(2,5,6,10,11)]
    ABC_pep <- read.delim("/Users/Pedro/Documents/2_Annotation/pep/Trinity_pep_AgBaseCommunity/pep_AgBaseCommunity_reformat.xls",header=TRUE)
    ABC_pep <- ABC_pep[c(2,5,6,10,11)]
    Inv_pep <- read.delim("/Users/Pedro/Documents/2_Annotation/pep/Trinity_pep_Invertebrates/pep_Invertebrates_reformat.xls",header=TRUE)
    Inv_pep <- Inv_pep[c(2,5,6,10,11)]
    SwP_pep <- read.delim("/Users/Pedro/Documents/2_Annotation/pep/Trinity_pep_SwissProt/pep_SwissProt_reformat.xls",header=TRUE)
    SwP_pep <- SwP_pep[c(2,5,6,10,11)]
    UnP_pep <- read.delim("/Users/Pedro/Documents/2_Annotation/pep/Trinity_pep_UniProt/pep_UniProt_reformat.xls",header=TRUE)
    UnP_pep <- UnP_pep[c(2,5,6,10,11)]
    
    ABC <- merge(ABC_dna,ABC_pep, all=TRUE)
    Inv <- merge(Inv_dna,Inv_pep, all=TRUE)
    SwP <- merge(SwP_dna,SwP_pep, all=TRUE)
    UnP <- merge(UnP_dna,UnP_pep, all=TRUE)
    ABCeInv <- merge(ABC,Inv, all=TRUE)
    SwPeUnP <- merge(UnP,SwP, all=TRUE)
    Final <- merge(ABCeInv,SwPeUnP, all=TRUE)
    
    write.table(Final,file="/Users/Pedro/Documents/2_Annotation/Final.txt",quote=FALSE,sep="\t",row.names=FALSE)

The results were uploaded to the server using the command 

    scp -r /Users/Pedro/Documents/2_Annotation raquel@143.107.247.2:/mnt/HD5/Chom-resistance/Anotacao/5-AgBase

**GO Enrichment Tests**

The annotated DE contigs were compared to the assembled transcriptome on the GOrilla (http://cbl-gorilla.cs.technion.ac.il/)  websites for the enrichment tests.

Pathway to the annotated DE table

    cd /mnt/HD5/Chom-resistance/DE_Annotation/Final_DE_table.xlsx

Pathway to the assembled transcriptome

    cd /mnt/HD5/Chom-resistance/Anotacao/5-AgBase/2_Annotation/Final.txt


- **GOrilla**

*Parameters*
      - Step 1 (Choose organism): *Drosophila melanogaster*
      - Step 2 (Choose running mode): two unranked lists of genes
      - Step 3 (Paste a ranked list of get/protein names): Target = list of the gene names (DB_object_name column) of the Final_DE_table.xlsx and 
                                                                                               Background = list of the gene names (DB_object_name column) of the Final.txt 
      The repeated gene names on the Final.txt file had to be excluded because the file was too big      to go on GOrilla
      - Step 4 (Choose an ontology): All
      
      Advanced parameters weren’t used
      
  *Results*    
      46,2% of the terms were used in the analysis 
      The system has recognized 4665 genes out of 8314 gene terms entered by the user.
     4665 genes were recognized by gene symbol and 0 genes by other gene IDs .
     672 duplicate genes were removed (keeping the highest ranking instance of each gene) leaving a total of 3993 genes.
     Only 3841 of these genes are associated with a GO term.
          
     
PROCESS    

![](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOPROCESS.png)


[http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOPROCESS.png](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOPROCESS.png)

| GO term    | Description                      | P-value | FDR q-value | Enrichment <br>(N, B, n, b) | Genes                                                                                                                                                                                                                                                       |
| ---------- | -------------------------------- | ------- | ----------- | --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| GO:0000154 | rRNA modification                | 1.96E-5 | 1.21E-1     | 22.86 (3841,14,48,4)        | **nop5** - cg10206 gene product from   transcript cg10206-ra<br>**CG8939** - cg8939 gene product from transcript cg8939-ra<br>**Nop60B** - nucleolar protein at 60b<br>**Fib** - fibrillarin                                                                |
| GO:0009451 | RNA modification                 | 5.97E-4 | 1E0         | 10.00 (3841,32,48,4)        | **nop5** - cg10206 gene product from   transcript cg10206-ra<br>**CG8939** - cg8939 gene product from transcript cg8939-ra<br>**Nop60B** - nucleolar protein at 60b<br>**Fib** - fibrillarin                                                                |
| GO:0006364 | rRNA processing                  | 7.44E-4 | 1E0         | 6.78 (3841,59,48,5)         | **nop5** - cg10206 gene product from transcript cg10206-ra<br>**CG8939** - cg8939 gene product from transcript cg8939-ra<br>**Nop60B** - nucleolar protein at 60b<br>**CG13185** - cg13185 gene product from transcript cg13185-rd<br>**Fib** - fibrillarin |
| GO:0040003 | chitin-based cuticle development | 8.67E-4 | 1E0         | 6.56 (3841,61,48,5)         | **Pcp** - pupal cuticle protein<br>**TwdlD** - tweedled<br>**TwdlF** - tweedlef<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1                                                                                               |
| GO:0016072 | rRNA metabolic process           | 9.35E-4 | 1E0         | 6.45 (3841,62,48,5)         | **nop5** - cg10206 gene product from transcript cg10206-ra<br>**CG8939** - cg8939 gene product from transcript cg8939-ra<br>**Nop60B** - nucleolar protein at 60b<br>**CG13185** - cg13185 gene product from transcript cg13185-rd<br>**Fib** - fibrillarin |


FUNCTION

![](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOFUNCTION.png)


[http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOFUNCTION.png](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOFUNCTION.png)

| GO term    | Description                                           | P-value | FDR q-value | Enrichment <br>(N, B, n, b) | Genes                                                                                                                                                         |
| ---------- | ----------------------------------------------------- | ------- | ----------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| GO:0008010 | structural constituent of chitin-based larval cuticle | 2.94E-7 | 5.8E-4      | 30.78 (3841,13,48,5)        | **Pcp** - pupal cuticle protein<br>**TwdlD** - tweedled<br>**TwdlF** - tweedlef<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1 |
| GO:0005214 | structural constituent of chitin-based cuticle        | 6.74E-7 | 6.64E-4     | 26.67 (3841,15,48,5)        | **Pcp** - pupal cuticle protein<br>**TwdlD** - tweedled<br>**TwdlF** - tweedlef<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1 |
| GO:0042302 | structural constituent of cuticle                     | 6.74E-7 | 4.43E-4     | 26.67 (3841,15,48,5)        | **Pcp** - pupal cuticle protein<br>**TwdlD** - tweedled<br>**TwdlF** - tweedlef<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1 |
| GO:0004497 | monooxygenase activity                                | 9.15E-4 | 4.51E-1     | 15.00 (3841,16,48,3)        | **Cyp6g2** - cg8859 gene product from transcript cg8859-ra<br>**Hn** - henna<br>**Cyp4d21** - cg6730 gene product from transcript cg6730-ra                   |



COMPONENT

![](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOCOMPONENT.png)


[http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOCOMPONENT.png](http://cbl-gorilla.cs.technion.ac.il/GOrilla/ofezxcmb/GOCOMPONENT.png)

| GO term    | Description                               | P-value | FDR q-value | Enrichment <br>(N, B, n, b) | Genes                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ---------- | ----------------------------------------- | ------- | ----------- | --------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| GO:0044421 | extracellular region part                 | 6.98E-6 | 7.3E-3      | 4.53 (3841,212,48,12)       | **Pcp** - pupal cuticle protein<br>**AttA** - attacin-a<br>**Npc2e** - niemann-pick type c-2e<br>**Peritrophin-15a** - cg17814 gene product from transcript cg17814-ra<br>**Tsf1** - transferrin 1<br>**TwdlD** - tweedled<br>**DptB** - diptericin b<br>**Cht4** - chitinase 4<br>**TwdlF** - tweedlef<br>**Obp99b** - odorant-binding protein 99b<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1 |
| GO:0031012 | extracellular matrix                      | 1.17E-4 | 6.11E-2     | 10.00 (3841,40,48,5)        | **Pcp** - pupal cuticle protein<br>**TwdlD** - tweedled<br>**TwdlF** - tweedlef<br>**Lcp2** - larval cuticle protein 2<br>**Lcp1** - larval cuticle protein 1                                                                                                                                                                                                                                                                     |
| GO:0005732 | small nucleolar ribonucleoprotein complex | 2.82E-4 | 9.83E-2     | 21.82 (3841,11,48,3)        | **nop5** - cg10206 gene product from transcript cg10206-ra<br>**Nop60B** - nucleolar protein at 60b<br>**Fib** - fibrillarin                                                                                                                                                                                                                                                                                                      |
| GO:0030684 | preribosome                               | 4.06E-4 | 1.06E-1     | 11.04 (3841,29,48,4)        | **nop5** - cg10206 gene product from transcript cg10206-ra<br>**CG8939** - cg8939 gene product from transcript cg8939-ra<br>**CG13185** - cg13185 gene product from transcript cg13185-rd<br>**Fib** - fibrillarin                                                                                                                                                                                                                |


Table informations:
'P-value' is the enrichment p-value computed according to the mHG or HG model. This p-value is not corrected for multiple testing of 6176 GO terms.

'FDR q-value' is the correction of the above p-value for multiple testing using the Benjamini and Hochberg (1995) method. 
Namely, for the ith term (ranked according to p-value) the FDR q-value is (p-value * number of GO terms) / i. 

Enrichment (N, B, n, b) is defined as follows:
N - is the total number of genes
B - is the total number of genes associated with a specific GO term
n - is the number of genes in the top of the user's input list or in the target set when appropriate
b - is the number of genes in the intersection
Enrichment = (b/n) / (B/N)

Genes: For each GO term you can see the list of associated genes that appear in the optimal top of the list. Each gene name is specified by gene symbol followed by a short description of the gene


08.05.2017

# 11 Variant analysis 
## 11.1 KisSplice + KisSplice2RefTranscriptome + KissDE

**11.1.1 KisSplice**
Version 2.4.0
A local transcriptome assembler for SNPs, indels and AS events.

Gustavo AT Sacomoto, Janice Kielbassa, Rayan Chikhi, Raluca Uricaru, Pavlos Antoniou, Marie-France Sagot, Pierre Peterlongo* and Vincent Lacroix*, KISSPLICE: de-novo calling alternative splicing events from RNA-seq data, BMC Bioinformatics 2012, 13(Suppl 6):S5

Parameters:
-s OUTPUT_SNPS 0, 1 or 2. 
Changes which types of SNPs will be output. If 0 (default), will not output SNPs. If 1, will output Type0a-SNPs. If 2, will output Type0a and Type0b SNPs (warning: this option may increase a lot the running time)
KisSplice may output 2 files for SNPs events, one for Type0a and another for Type0b. Type0a indicates that there is a single SNP within the sequence, while Type0b would indicate the presence of multiple SNPs, i.e. several SNPs separated by less than k nt, or a pattern created by paralogous genes.

-t NBPROCS Number of cores (must be <= number of physical cores)

-o OUT_DIR Path to store the results and the summary log file (default = ./results)

-r READFILES Input fasta/q read files or compressed (.gz) fasta/q files (mutiple, such as "-r file1 -r file2...")


*Running KisSplice (all reads R1 and R2 files concatenated)*

    ulimit -s hard
    
    ./kissplice -s 1 -t 12 -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/KisSplice -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Control-1st-R1-trimmed.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Resistant-1st-R1-trimmed.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R1R2-trimmed.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R1R2-trimmed.fastq 

*Running KisSplice with just R1 (of second replicate)*

    ulimit -s hard
    
    ./kissplice -s 1 -t 12 -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1 -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Control-1st-R1-trimmed.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/1st-replicate/Resistant-1st-R1-trimmed.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Control-2nd-R1-trimmed-paired.fastq -r /Volumes/HD3/analyses/Tatiana/Chom-resistance/2-trimming/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq


KisSplice results *(all reads)*

                      ******** We are done, final coherent results are in files /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/KisSplice/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1R2-trimmed_Resistant-2nd-R1R2-trimmed_k41_coherents_type_*.fa ********** 
                      ******** All non read coherent results are in files /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/KisSplice/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1R2-trimmed_Resistant-2nd-R1R2-trimmed_k41_uncoherent.fa ****** 
    
    
                     TYPES:
                             0a: Single SNPs, Inexact Repeats or sequencing substitution errors, 110755 found
                             0b: Run with -s 2 to also search for Multiple SNPs (warning: this option may increase a lot the running time. You might also want to try the experimental algorithm here)
                             1: Alternative Splicing Events, 7671 found
                             2: Inexact Tandem Repeats, 133 found
                             3: Short Indels (<3nt), 21208 found
                             4: All others, composed by a shorter path of length > 2k not being a SNP, 2569 found
    
    
                      ******** A summary of the execution can be found in the log file: /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/KisSplice/kissplice_log_summary_17-05-13_05-09-2019_853865**********


KisSplice results *(R1)*

           ******** We are done, final coherent results are in files /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1-trimmed-paired_Resistant-2nd-R1-trimmed-paired_k41_coherents_type_*.fa ********** 
                      ******** All non read coherent results are in files /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1-trimmed-paired_Resistant-2nd-R1-trimmed-paired_k41_uncoherent.fa ****** 
    
    
                     TYPES:
                             0a: Single SNPs, Inexact Repeats or sequencing substitution errors, 77087 found
                             0b: Run with -s 2 to also search for Multiple SNPs (warning: this option may increase a lot the running time. You might also want to try the experimental algorithm here)
                             1: Alternative Splicing Events, 4586 found
                             2: Inexact Tandem Repeats, 121 found
                             3: Short Indels (<3nt), 12256 found
                             4: All others, composed by a shorter path of length > 2k not being a SNP, 1599 found
                             
           ******** A summary of the execution can be found in the log file: /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1/kissplice_log_summary_18-28-01_16-05-2017_90729**********

Renaming files 

    for f in *Control*; do mv -v "$f" "${f/Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1-trimmed-paired_Control-2nd-R2-trimmed-paired_Resistant-2nd-R1-trimmed-paired_Resistant-2nd-R2-trimmed-paired/Chom}"; done

Renaming files R1

    for f in *Control*; do mv -v "$f" "${f/Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1-trimmed-paired_Resistant-2nd-R1-trimmed-paired/Chom}"; done



**11.1.2 Transdecoder** 
Version 3.0.1
Transdecoder will predict the ORFs of the transcripts that will be used for the prediction of the functional impact of SNPs. 

*Running Transdecoder*

    ./TransDecoder.LongOrfs -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta 
    
    perl gff3_file_to_bed.pl /Applications/seq_an/TransDecoder-3.0.1/Trinity.Trinity.fasta.transdecoder_dir/Chom-longest_orfs.gff3 >/Applications/seq_an/TransDecoder-3.0.1/Trinity.Trinity.fasta.transdecoder_dir/Chom-longest_orfs.bed

02.09.2019 *Running Transdecoder on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    /Applications/seq_an/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta  
    
    perl /Applications/seq_an/TransDecoder-3.0.1/util/gff3_file_to_bed.pl /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Transdecoder/Trinity.longest.clusters.fasta.transdecoder_dir/longest_orfs.gff3 >/Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Transdecoder/longest_orfs.bed


11.05.2017

**11.1.3 BLAT** 
Standalone BLAT v. 34 fast sequence search command line tool

KisSplice outputs two files for SNPs events, one for Type0a and another for Type0b. Type0a indicates that there is a single SNP within the sequence, while Type0b would indicate the presence of multiple SNPs, i.e. several SNPs separated by less than k nt, or a pattern created by paralogous genes. Only Type0a were used in the analyses.

*Running BLAT (R1_R2)*

    ./blat -minIdentity=80 /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSplice/results_Chom_k41_coherents_type_0a.fa /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Chom-blat-output.psl 

06.09.2019 *Running BLAT (R1_R2) on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    /Applications/seq_an/blatSuite/blat -minIdentity=80 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1R2-trimmed_Resistant-2nd-R1R2-trimmed_k41_coherents_type_0a.fa /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/Chom-blat-outputR1R2.psl 

*Running BLAT (R1)*

    ./blat -minIdentity=80 /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1_only/1-KisSplice/results_Chom_k41_coherents_type_0a.fa /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Chom-blat-outputR1.psl 

02.09.2019 *Running BLAT (R1) on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    /Applications/seq_an/blatSuite/blat -minIdentity=80 /Volumes/HD3/analyses/Tatiana/Chom-resistance/6-read_mapping/databases/Trinity.longest.clusters.fasta /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/KisSpliceResults_Chom_k41_coherents_type_0a.fa /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/Chom-blat-outputR1only.psl


12.05.2017

**11.1.4 KissDE** 
Version 1.5.0
KissDE is a R package which works on pairs of variants, and tests if a variant is enriched in one condition. 

*Running kissDE in R*

    library(kissDE)  
    snp<-kissplice2counts("/Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSplice/results_Chom_k41_coherents_type_0a.fa”, pairedEnd=FALSE)
    conditions<-c("control","resist","control","control","resist","resist") 
    res<-diffExpressedVariants(snp, conditions, pvalue=1) 
    write.table(res$finalTable, file="kissDE_output", sep="\t", quote=FALSE)

06.09.2019 *Running kissDE in R* *on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    library(kissDE)  
    snp<-kissplice2counts("/Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1R2-trimmed_Resistant-2nd-R1R2-trimmed_k41_coherents_type_0a.fa", pairedEnd=FALSE)
    conditions<-c("control","resist","control","resist") 
     
    res<-diffExpressedVariants(snp, conditions, pvalue=1) 
    write.table(res$finalTable, file="kissDE_output", sep="\t", quote=FALSE) 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=0.05) 
    write.table(res$finalTable, file="kissDE_output05", sep="\t", quote=FALSE) 
    

*Running kissDE in R (R1)*

    library(kissDE) 
    snp<-kissplice2counts("/Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSplice/results_Chom_k41_coherents_type_0a.fa”, pairedEnd=FALSE)
    conditions<-c("control","resist","control","resist") 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=1) 
    write.table(res$finalTable, file="kissDE_output", sep="\t", quote=FALSE) 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=0.05) 
    write.table(res$finalTable, file="kissDE_output05", sep="\t", quote=FALSE) 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=0.01) 
    write.table(res$finalTable, file="kissDE_output01", sep="\t", quote=FALSE) 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=0.001) 
    write.table(res$finalTable, file="kissDE_output001", sep="\t", quote=FALSE) 

02.09.2019 *Running kissDE in R* *(R1) on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    library(kissDE)  
    snp<-kissplice2counts("/Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/KisSpliceResults_Chom_k41_coherents_type_0a.fa", pairedEnd=FALSE)
    conditions<-c("control","resist","control","resist") 
     
    res<-diffExpressedVariants(snp, conditions, pvalue=1) 
    write.table(res$finalTable, file="kissDE_output", sep="\t", quote=FALSE) 
    
    res<-diffExpressedVariants(snp, conditions, pvalue=0.05) 
    write.table(res$finalTable, file="kissDE_output05", sep="\t", quote=FALSE) 
    


17.05.2017

**11.1.5 KisSplice2reftranscriptome**
Version 1.2.2
Getting information on the functional impact of a SNP on a protein, without using a reference genome

-b Transdecoder file .bed 
-k KisSplice results
-t Blat output
-s kissDE output
-o path to output file 
-pval p-value threshold
-Q minimum query coverage (default: 90). A query coverage under 70% isn’t recommended.

*Running kisSplice2reftranscriptome (R1+R2)*

    python /Applications/seq_an/kissplice2reftranscriptome/kissplice2reftranscriptome -b /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Transdecoder/Chom-longest_orfs.bed -k /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1R2/results_Chom_k41_coherents_type_0a.fa -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Chom-blat-output.psl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/KisSpliceR1R2/kissDE_output 

06.09.2019 *Running kisSplice2reftranscriptome (R1+R2) on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    python /Applications/seq_an/kissplice2reftranscriptome/kissplice2reftranscriptome -b /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Transdecoder/longest_orfs.bed -k /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/1-KisSplice/results_Control-1st-R1-trimmed_Resistant-1st-R1-trimmed_Control-2nd-R1R2-trimmed_Resistant-2nd-R1R2-trimmed_k41_coherents_type_0a.fa -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/Chom-blat-outputR1R2.psl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1R2/2-KissDE/kissDE_output -o ChomR1R2-k2t.tsv 

*Running kisSplice2reftranscriptome (R1)*

    python /Applications/seq_an/kissplice2reftranscriptome/kissplice2reftranscriptome -b /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1_only/Chom-longest_orfs.bed -k /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1_only/1-KisSplice/results_Chom_k41_coherents_type_0a.fa -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1_only/Chom-blat-outputR1.psl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1_only/2-kissDE/kissDE_output0.001 -o k2t001.tsv

06.09.2019 *Running kisSplice2reftranscriptome (R1 only) on* Trinity.longest.clusters.fasta (same transcripts used for DE analyses)

    python /Applications/seq_an/kissplice2reftranscriptome/kissplice2reftranscriptome -b /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/Transdecoder/longest_orfs.bed -k /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/1-KisSplice/results_Chom_k41_coherents_type_0a.fa -t /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/Chom-blat-outputR1only.psl -s /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/2-KissDE/kissDE_output -o /Volumes/HD3/analyses/Tatiana/Chom-resistance/11-variant_analysis/R1only/3-k2t/ChomR1only-k2t.tsv 

19.05.2017

**11.1.6 KissDE + KisSplice2reftranscriptome**
 **Sorting and joining files

    sed 's/\|Type_0a//' k2t001.tsv >temp.tsv
    
    sort kissDE_output0.001 >kissDE_output001.sorted 
    sort -k 2 temp.tsv >k2t001.sorted
    
    join -1 1 -2 2 -t $'\t'
    

Results: 5060 variants with different frequencies between conditions. 

Condition1 = Control; Condition2 = Resistant

deltaPSI* = frequency_condition2  - frequency_condition1
*for variant1

DeltaPSI is an average between replicates and may be corrected by length (see KissDE manual).

Positive number = increased frequency of variant 1 in the resistant condition

Of the 5060 variants, 4648 were aligned to transdecoder ORFs (412 had no contig aligned).

4012 variants were found in coding sequences. Of these, 981 were non-synonymous. 

06.09.2019

**11.1.7 KissDE + KisSplice2reftranscriptome (***on* Trinity.longest.clusters.fasta**)**
 **Sorting and joining files R1 only

    sed 's/\|Type_0a//' ChomR1only-k2t.tsv >temp.tsv
    
    sort kissDE_output >kissDE_output.sorted 
    sort -k 2 temp.tsv >ChomR1only-k2t.sorted
    
    join -1 1 -2 2 -t $'\t' kissDE_output.sorted ChomR1only-k2t.sorted >ChomR1only-kissDE-k2t.result

Results pval 0.05:

- 14566 variants
- 12209 variants on ORFs
- 3106 non synonymous variants

Results pval 0.01:

- 8928 variants
- 7521 variants on ORFs
- 1892 non synonymous variants

Results pval 0.001:

- 5244 variants
- 4441 variants on ORFs
- 1100 non synonymous variants

**11.1.8 KissDE + KisSplice2reftranscriptome R1 + R2 (***on* Trinity.longest.clusters.fasta**)**
Sorting and joining files R1 and  R2

    sed 's/\|Type_0a//' ChomR1R2-k2t.tsv >temp.tsv
    
    sort kissDE_output >kissDE_output.sorted 
    sort -k 2 temp.tsv >ChomR1R2-k2t.sorted
    
    join -1 1 -2 2 -t $'\t' kissDE_output.sorted ChomR1only-k2t.sorted >ChomR1R2-k2t.sorted >ChomR1R2-kissDE-k2t.result

Results pval 0.05:

- 18752 variants
- 15681 variants on ORFs
- 4025 non synonymous variants
- 2478 with coverage of at least 30 in both conditions

Results pval 0.01:

- 11363 variants
- 9574 variants on ORFs
- 2436 non synonymous variants
- 1777 with coverage of at least 30 in both conditions

Results pval 0.001:

- 6648 variants
- 5620 variants on ORFs
- 1395 non synonymous variants
- 1169 with coverage of at least 30 in both conditions


## 11.2. Analysis of the condition-specific SNPs 

On the new analysis made using the reduced transcriptome and using both R1 and R2 of the second replicate

**11.2.1. Determining thresholds and filters (16/09/2019)**

Using P_value 0.05:

|            | difference in allele usage 15% | difference in allele usage 20% | differnce in allele usage 25% | Difference in allele usage 30% |
| ---------- | ------------------------------ | ------------------------------ | ----------------------------- | ------------------------------ |
| Cov 30     | 1955 (986)                     | 1484 (809)                     | 1087                          | 770                            |
| **Cov 50** | **1423 (700)**                 | **1042 (560)**                 | 747                           | 522                            |

Using P_value 0.01: Number of SNPs (number of transcript)

|            | difference in allele usage 15% | difference in allele usage 20% | differnce in allele usage 25% | Difference in allele usage 30% |
| ---------- | ------------------------------ | ------------------------------ | ----------------------------- | ------------------------------ |
| Cov 30     | 1463 (810)                     | 1150 (681)                     | 867                           | 655                            |
| **Cov 50** | **1112 (607)**                 | **834 (490)**                  | 619                           | 454                            |


=> At first, we thought that we could chose a “per condition”covereage cut-off. However, by looking at the table, I saw that there was a lack of coherence between replicates (one replicate would have very low coverage compared to the other). This lack of conguence is hard to filter for because taking only the difference in coverages would not be correct as the same difference considering small numbers is not comparable to large numbers. Example:

| Cov rep1 | Cov rep2 |                                                                                                            |
| -------- | -------- | ---------------------------------------------------------------------------------------------------------- |
| 0        | 30       | Difference is 30 → Rep1 an Rep2 are different although they should be the same                             |
| 1000     | 1030     | Difference is 30 → Same difference but at this high level of expression 30 of difference is not a big deal |

A simple solution to get rid of the incongruence in lowly expressed locus is to force each replicate to reach a certain leval of coverage. We chose 20. 


=> we chose to analyse 2 sets of condition-specific non-synonymous SNPs with:

|                                        | More lenient                                                                                 | More stringent |
| -------------------------------------- | -------------------------------------------------------------------------------------------- | -------------- |
| P_value                                | 0.01 (considering we only have 2 replicate we are being a bit more stringent on the P-value) | 0.01           |
| Minimun coverage (per replicate)       | 20                                                                                           | 20             |
| difference in allele usage (magnitude) | 15%                                                                                          | 20%            |

 
**SNPs with Pval0.01 , Cov20 and diff15:**

We have 1095 condition-specific SNPs which satisfy all criteria in 614 transcripts


    cut -f 27 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15.csv | grep -v "^#" | wc -l
    1097
    
    cut -f 27 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15.csv | grep -v "^#" | LANG=en_EN sort | uniq | wc -l
    614

**SNPs with Pval0.01 , Cov20 and diff20:**

We have 830 condition-specific SNPs in 492 transcripts


    cut -f 1 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff20-with-annotation.txt | grep -v "^#" | LANG=en_EN sort | wc -l
    830
    > 16:56 sophie@galaxy ~/w/R/1/SNPs-cov20-diff20-P001 $ cut -f 1 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff20-with-annotation.txt | grep -v "^#" | uniq | LANG=en_EN sort | wc -l
    492

**11.2.2. Annotation of the condition-specific SNPs (17/09/2019)**

The annotation of Trinity.longest.clusters.fasta was done using FunctionAnnotator. By sorting and Joining the file of condition-specific SNPs and the annotation table, I built a table of annotated condition-specific SNPs:


    # Getting the contig IF column
    cut -f 27 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15.csv > temp1.txt
    
    # Getting the other columns of interested (I have omited some columns from the KisSplice/DE pipeline which didn't add much info + those which I created to filter)
    cut -f 1,4-11,20,25,30-34,39-41 ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15.csv > temp2.txt
    
    # Pasting to have the contig IDs as first column
    paste temp1.txt temp2.txt > ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-col-subset.csv
    
    # Sorting
    LANG=en_EN sort ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-col-subset.csv > ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-col-subset.sorted
    
    # NOTE: we now have the first part of the table to assemble
    
    # Now we get the important columns of the annotation file + we sort them:
    cut -f 1-10,22-24 ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters.txt | LANG=en_EN sort > ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt
    
    # NOTE: I renamed beforehand the contig ID column of "AnnotationTable_trinity_longest_clusters.txt" to match the header of the first file
    
    # Joining both files to get the final table of condition-specific SNPs with annotation information.
    
    join -t '        ' ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-col-subset.sorted ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt > ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-with-annotation.txt

**11.2.3. GO analysis (17/09/2019)**

****FunctionAnnotator uses BLAST2GO to fetch the GO terms associated with the BlastHits and outputs a GO table that is BLAST2GO compatible. I performed GO analysis using Fisher’s exact test, implemented in BLAST2GO. 

**SNPs with Pval0.01 , Cov20 and diff15:**

over-represented GO terms in the transcript set containing condition-specific SNPs (reduced to most specific):

| Tags   | GO ID      | GO Name                                                                                               | GO Category        | FDR                  | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | ----------------------------------------------------------------------------------------------------- | ------------------ | -------------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0008061 | chitin binding                                                                                        | MOLECULAR_FUNCTION | 2.39576127358761E-18 | 5.15751510021252E-21 | 30      | 136          | 584            | 32293               |
| [OVER] | GO:0005615 | extracellular space                                                                                   | CELLULAR_COMPONENT | 9.10320250679736E-13 | 2.68552457942515E-15 | 31      | 251          | 583            | 32178               |
| [OVER] | GO:0004177 | aminopeptidase activity                                                                               | MOLECULAR_FUNCTION | 1.43149298917552E-09 | 5.25025334891355E-12 | 13      | 37           | 601            | 32392               |
| [OVER] | GO:0004252 | serine-type endopeptidase activity                                                                    | MOLECULAR_FUNCTION | 5.38707816293185E-09 | 2.06171066672563E-11 | 23      | 193          | 591            | 32236               |
| [OVER] | GO:0016490 | structural constituent of peritrophic membrane                                                        | MOLECULAR_FUNCTION | 1.1291572205463E-07  | 5.40180459518244E-10 | 7       | 5            | 607            | 32424               |
| [OVER] | GO:0051287 | NAD binding                                                                                           | MOLECULAR_FUNCTION | 6.63988191500243E-06 | 4.07647031936842E-08 | 9       | 31           | 605            | 32398               |
| [OVER] | GO:0010951 | negative regulation of endopeptidase activity                                                         | BIOLOGICAL_PROCESS | 2.81093680005029E-05 | 1.83779953439742E-07 | 11      | 66           | 603            | 32363               |
| [OVER] | GO:0044262 | cellular carbohydrate metabolic process                                                               | BIOLOGICAL_PROCESS | 0.000220525229608    | 1.84620866758541E-06 | 10      | 68           | 604            | 32361               |
| [OVER] | GO:0045454 | cell redox homeostasis                                                                                | BIOLOGICAL_PROCESS | 0.0003036443977      | 2.66312260779924E-06 | 8       | 40           | 606            | 32389               |
| [OVER] | GO:0032934 | sterol binding                                                                                        | MOLECULAR_FUNCTION | 0.000417837640721    | 3.79791827796329E-06 | 5       | 9            | 609            | 32420               |
| [OVER] | GO:0006108 | malate metabolic process                                                                              | BIOLOGICAL_PROCESS | 0.000427728355085    | 3.95252432005996E-06 | 4       | 3            | 610            | 32426               |
| [OVER] | GO:0033961 | cis-stilbene-oxide hydrolase activity                                                                 | MOLECULAR_FUNCTION | 0.000640675432902    | 6.38529972194932E-06 | 3       | 0            | 611            | 32429               |
| [OVER] | GO:0005581 | collagen trimer                                                                                       | CELLULAR_COMPONENT | 0.000780836205309    | 8.09350236726288E-06 | 7       | 33           | 607            | 32396               |
| [OVER] | GO:0016811 | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides               | MOLECULAR_FUNCTION | 0.000949668042687    | 1.00706306551898E-05 | 8       | 49           | 606            | 32380               |
| [OVER] | GO:0004867 | serine-type endopeptidase inhibitor activity                                                          | MOLECULAR_FUNCTION | 0.001320938434958    | 1.48502885767116E-05 | 8       | 52           | 606            | 32377               |
| [OVER] | GO:0031090 | organelle membrane                                                                                    | CELLULAR_COMPONENT | 0.001501930764422    | 1.7483805741153E-05  | 26      | 526          | 588            | 31903               |
| [OVER] | GO:0030246 | carbohydrate binding                                                                                  | MOLECULAR_FUNCTION | 0.00153567410826     | 1.81214932245628E-05 | 13      | 156          | 601            | 32273               |
| [OVER] | GO:0006096 | glycolytic process                                                                                    | BIOLOGICAL_PROCESS | 0.001707399304209    | 2.11008525077639E-05 | 7       | 39           | 607            | 32390               |
| [OVER] | GO:0005604 | basement membrane                                                                                     | CELLULAR_COMPONENT | 0.001915107527693    | 2.41258961390155E-05 | 6       | 26           | 608            | 32403               |
| [OVER] | GO:0030060 | L-malate dehydrogenase activity                                                                       | MOLECULAR_FUNCTION | 0.001986759596444    | 2.51869539016563E-05 | 3       | 1            | 611            | 32428               |
| [OVER] | GO:0005811 | lipid droplet                                                                                         | CELLULAR_COMPONENT | 0.002082963834019    | 2.65726529614869E-05 | 9       | 75           | 605            | 32354               |
| [OVER] | GO:0006457 | protein folding                                                                                       | BIOLOGICAL_PROCESS | 0.002106995858962    | 2.70472279774299E-05 | 10      | 95           | 604            | 32334               |
| [OVER] | GO:0004568 | chitinase activity                                                                                    | MOLECULAR_FUNCTION | 0.003077138277864    | 4.12182451507807E-05 | 6       | 29           | 608            | 32400               |
| [OVER] | GO:0007586 | digestion                                                                                             | BIOLOGICAL_PROCESS | 0.003277817275463    | 4.41676861388315E-05 | 5       | 17           | 609            | 32412               |
| [OVER] | GO:0006032 | chitin catabolic process                                                                              | BIOLOGICAL_PROCESS | 0.003550330317844    | 4.86889503005261E-05 | 6       | 30           | 608            | 32399               |
| [OVER] | GO:0004471 | malate dehydrogenase (decarboxylating) (NAD+) activity                                                | MOLECULAR_FUNCTION | 0.004302727338206    | 6.2094853150643E-05  | 3       | 2            | 611            | 32427               |
| [OVER] | GO:0030976 | thiamine pyrophosphate binding                                                                        | MOLECULAR_FUNCTION | 0.004302727338206    | 6.2094853150643E-05  | 3       | 2            | 611            | 32427               |
| [OVER] | GO:0042445 | hormone metabolic process                                                                             | BIOLOGICAL_PROCESS | 0.004533415690129    | 6.68698694525442E-05 | 6       | 32           | 608            | 32397               |
| [OVER] | GO:0005506 | iron ion binding                                                                                      | MOLECULAR_FUNCTION | 0.004599240080892    | 6.85163217281208E-05 | 13      | 179          | 601            | 32250               |
| [OVER] | GO:0015035 | protein disulfide oxidoreductase activity                                                             | MOLECULAR_FUNCTION | 0.004611932848741    | 6.91311892491869E-05 | 5       | 19           | 609            | 32410               |
| [OVER] | GO:0006098 | pentose-phosphate shunt                                                                               | BIOLOGICAL_PROCESS | 0.004825393104992    | 7.38698354455874E-05 | 4       | 9            | 610            | 32420               |
| [OVER] | GO:0016853 | isomerase activity                                                                                    | MOLECULAR_FUNCTION | 0.005023367632647    | 7.89031592753489E-05 | 11      | 132          | 603            | 32297               |
| [OVER] | GO:0008010 | structural constituent of chitin-based larval cuticle                                                 | MOLECULAR_FUNCTION | 0.006356212458467    | 0.000103892804496    | 6       | 35           | 608            | 32394               |
| [OVER] | GO:0016747 | transferase activity, transferring acyl groups other than amino-acyl groups                           | MOLECULAR_FUNCTION | 0.007291487980353    | 0.000120924055168    | 12      | 164          | 602            | 32265               |
| [OVER] | GO:0005783 | endoplasmic reticulum                                                                                 | CELLULAR_COMPONENT | 0.007396141266522    | 0.000123249364113    | 17      | 304          | 597            | 32125               |
| [OVER] | GO:0005777 | peroxisome                                                                                            | CELLULAR_COMPONENT | 0.007926795701305    | 0.000135252294696    | 7       | 54           | 607            | 32375               |
| [OVER] | GO:0004035 | alkaline phosphatase activity                                                                         | MOLECULAR_FUNCTION | 0.007986789218895    | 0.000136912747733    | 4       | 11           | 610            | 32418               |
| [OVER] | GO:0005759 | mitochondrial matrix                                                                                  | CELLULAR_COMPONENT | 0.00840405111018     | 0.000144735691261    | 8       | 74           | 606            | 32355               |
| [OVER] | GO:0006094 | gluconeogenesis                                                                                       | BIOLOGICAL_PROCESS | 0.010161962082948    | 0.000179872076416    | 4       | 12           | 610            | 32417               |
| [OVER] | GO:0004497 | monooxygenase activity                                                                                | MOLECULAR_FUNCTION | 0.010185977137796    | 0.000181109304874    | 11      | 146          | 603            | 32283               |
| [OVER] | GO:0016903 | oxidoreductase activity, acting on the aldehyde or oxo group of donors                                | MOLECULAR_FUNCTION | 0.011054739978876    | 0.000200081803158    | 6       | 40           | 608            | 32389               |
| [OVER] | GO:0044445 | cytosolic part                                                                                        | CELLULAR_COMPONENT | 0.011054739978876    | 0.000200081803158    | 6       | 40           | 608            | 32389               |
| [OVER] | GO:0030055 | cell-substrate junction                                                                               | CELLULAR_COMPONENT | 0.011054739978876    | 0.000200081803158    | 6       | 40           | 608            | 32389               |
| [OVER] | GO:0070008 | serine-type exopeptidase activity                                                                     | MOLECULAR_FUNCTION | 0.011478214498933    | 0.000211358811988    | 3       | 4            | 611            | 32425               |
| [OVER] | GO:0004181 | metallocarboxypeptidase activity                                                                      | MOLECULAR_FUNCTION | 0.011478214498933    | 0.000211407076164    | 5       | 25           | 609            | 32404               |
| [OVER] | GO:0006662 | glycerol ether metabolic process                                                                      | BIOLOGICAL_PROCESS | 0.012475712871083    | 0.000231768545604    | 4       | 13           | 610            | 32416               |
| [OVER] | GO:0010038 | response to metal ion                                                                                 | BIOLOGICAL_PROCESS | 0.015080823432298    | 0.000284974896624    | 6       | 43           | 608            | 32386               |
| [OVER] | GO:0007374 | posterior midgut invagination                                                                         | BIOLOGICAL_PROCESS | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0010765 | positive regulation of sodium ion transport                                                           | BIOLOGICAL_PROCESS | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0006533 | aspartate catabolic process                                                                           | BIOLOGICAL_PROCESS | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0007320 | insemination                                                                                          | BIOLOGICAL_PROCESS | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0008445 | D-aspartate oxidase activity                                                                          | MOLECULAR_FUNCTION | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0008756 | o-succinylbenzoate-CoA ligase activity                                                                | MOLECULAR_FUNCTION | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0030478 | actin cap                                                                                             | CELLULAR_COMPONENT | 0.017294572694954    | 0.000344733150513    | 2       | 0            | 612            | 32429               |
| [OVER] | GO:0005773 | vacuole                                                                                               | CELLULAR_COMPONENT | 0.017294572694954    | 0.000338155766711    | 10      | 132          | 604            | 32297               |
| [OVER] | GO:0031967 | organelle envelope                                                                                    | CELLULAR_COMPONENT | 0.017294572694954    | 0.000342038633716    | 16      | 302          | 598            | 32127               |
| [OVER] | GO:0007155 | cell adhesion                                                                                         | BIOLOGICAL_PROCESS | 0.018432837968694    | 0.000374770665127    | 17      | 336          | 597            | 32093               |
| [OVER] | GO:0005201 | extracellular matrix structural constituent                                                           | MOLECULAR_FUNCTION | 0.018698533238568    | 0.000383154444452    | 7       | 65           | 607            | 32364               |
| [OVER] | GO:0008283 | cell population proliferation                                                                         | BIOLOGICAL_PROCESS | 0.020199638561242    | 0.000423577176017    | 18      | 372          | 596            | 32057               |
| [OVER] | GO:0006633 | fatty acid biosynthetic process                                                                       | BIOLOGICAL_PROCESS | 0.020876437776966    | 0.000439433868053    | 6       | 47           | 608            | 32382               |
| [OVER] | GO:0016328 | lateral plasma membrane                                                                               | CELLULAR_COMPONENT | 0.021363491596777    | 0.000451389353623    | 4       | 16           | 610            | 32413               |
| [OVER] | GO:0003779 | actin binding                                                                                         | MOLECULAR_FUNCTION | 0.021633549959093    | 0.000459648220304    | 13      | 220          | 601            | 32209               |
| [OVER] | GO:0009396 | folic acid-containing compound biosynthetic process                                                   | BIOLOGICAL_PROCESS | 0.022916795236834    | 0.000493345137454    | 3       | 6            | 611            | 32423               |
| [OVER] | GO:0009081 | branched-chain amino acid metabolic process                                                           | BIOLOGICAL_PROCESS | 0.022916795236834    | 0.000493345137454    | 3       | 6            | 611            | 32423               |
| [OVER] | GO:0016705 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | MOLECULAR_FUNCTION | 0.022916795236834    | 0.000491453836876    | 12      | 193          | 602            | 32236               |
| [OVER] | GO:0061783 | peptidoglycan muralytic activity                                                                      | MOLECULAR_FUNCTION | 0.029746894583061    | 0.000661727283422    | 4       | 18           | 610            | 32411               |
| [OVER] | GO:0001666 | response to hypoxia                                                                                   | BIOLOGICAL_PROCESS | 0.031625542539641    | 0.00071647605143     | 6       | 52           | 608            | 32377               |
| [OVER] | GO:0051082 | unfolded protein binding                                                                              | MOLECULAR_FUNCTION | 0.031625542539641    | 0.00071647605143     | 6       | 52           | 608            | 32377               |
| [OVER] | GO:0051259 | protein complex oligomerization                                                                       | BIOLOGICAL_PROCESS | 0.036100835771384    | 0.000843369867726    | 7       | 75           | 607            | 32354               |
| [OVER] | GO:0048514 | blood vessel morphogenesis                                                                            | BIOLOGICAL_PROCESS | 0.038490270852264    | 0.000906532740881    | 7       | 76           | 607            | 32353               |
| [OVER] | GO:0099512 | supramolecular fiber                                                                                  | CELLULAR_COMPONENT | 0.038490270852264    | 0.000908397398523    | 15      | 300          | 599            | 32129               |
| [OVER] | GO:0006873 | cellular ion homeostasis                                                                              | BIOLOGICAL_PROCESS | 0.039473745484371    | 0.000937902739144    | 6       | 55           | 608            | 32374               |
| [OVER] | GO:0070938 | contractile ring                                                                                      | CELLULAR_COMPONENT | 0.039536029966975    | 0.000942534919481    | 3       | 8            | 611            | 32421               |
| [OVER] | GO:0007375 | anterior midgut invagination                                                                          | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0031581 | hemidesmosome assembly                                                                                | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0038095 | Fc-epsilon receptor signaling pathway                                                                 | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0019343 | cysteine biosynthetic process via cystathionine                                                       | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0035074 | pupation                                                                                              | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0006535 | cysteine biosynthetic process from serine                                                             | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0071346 | cellular response to interferon-gamma                                                                 | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0009251 | glucan catabolic process                                                                              | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0035357 | peroxisome proliferator activated receptor signaling pathway                                          | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0031670 | cellular response to nutrient                                                                         | BIOLOGICAL_PROCESS | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0004467 | long-chain fatty acid-CoA ligase activity                                                             | MOLECULAR_FUNCTION | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0004122 | cystathionine beta-synthase activity                                                                  | MOLECULAR_FUNCTION | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0003878 | ATP citrate synthase activity                                                                         | MOLECULAR_FUNCTION | 0.040412494210268    | 0.001021428852229    | 2       | 1            | 612            | 32428               |
| [OVER] | GO:0044433 | cytoplasmic vesicle part                                                                              | CELLULAR_COMPONENT | 0.043706384626816    | 0.001118621389348    | 7       | 79           | 607            | 32350               |
| [OVER] | GO:0044242 | cellular lipid catabolic process                                                                      | BIOLOGICAL_PROCESS | 0.045467853215533    | 0.001170954918563    | 5       | 38           | 609            | 32391               |
| [OVER] | GO:0003073 | regulation of systemic arterial blood pressure                                                        | BIOLOGICAL_PROCESS | 0.047977486883365    | 0.001239412035577    | 3       | 9            | 611            | 32420               |

=> As this analysis gave quite a lot of results I was concerned that using a 15% difference in allele usage threshold was too low (that we were getting uninterresting SNPs)

**SNPs with Pval0.01 , Cov20 and diff20:**

Reduced set:

| Tags   | GO ID      | GO Name                                                                                               | GO Category        | FDR                  | P-Value              | Nr Test | Nr Reference | Non Annot Test | Non Annot Reference |
| ------ | ---------- | ----------------------------------------------------------------------------------------------------- | ------------------ | -------------------- | -------------------- | ------- | ------------ | -------------- | ------------------- |
| [OVER] | GO:0006030 | chitin metabolic process                                                                              | BIOLOGICAL_PROCESS | 1.13820986357857E-17 | 1.72428539371655E-20 | 28      | 151          | 464            | 32400               |
| [OVER] | GO:0008061 | chitin binding                                                                                        | MOLECULAR_FUNCTION | 2.53322918194451E-15 | 5.04949207053204E-18 | 25      | 141          | 467            | 32410               |
| [OVER] | GO:0005615 | extracellular space                                                                                   | CELLULAR_COMPONENT | 1.47359499528362E-12 | 3.28979906457833E-15 | 28      | 254          | 464            | 32297               |
| [OVER] | GO:0004177 | aminopeptidase activity                                                                               | MOLECULAR_FUNCTION | 4.12780925284847E-08 | 1.57977072346298E-10 | 11      | 39           | 481            | 32512               |
| [OVER] | GO:0004252 | serine-type endopeptidase activity                                                                    | MOLECULAR_FUNCTION | 1.82356073931711E-07 | 7.560609029221E-10   | 19      | 197          | 473            | 32354               |
| [OVER] | GO:0042302 | structural constituent of cuticle                                                                     | MOLECULAR_FUNCTION | 1.32100597767254E-06 | 6.21426827321639E-09 | 15      | 132          | 477            | 32419               |
| [OVER] | GO:0016490 | structural constituent of peritrophic membrane                                                        | MOLECULAR_FUNCTION | 1.89261612439187E-06 | 9.05413550179492E-09 | 6       | 6            | 486            | 32545               |
| [OVER] | GO:0010951 | negative regulation of endopeptidase activity                                                         | BIOLOGICAL_PROCESS | 3.91912132015324E-06 | 1.96862257351024E-08 | 11      | 66           | 481            | 32485               |
| [OVER] | GO:0032934 | sterol binding                                                                                        | MOLECULAR_FUNCTION | 0.000183150275275    | 1.2850601358774E-06  | 5       | 9            | 487            | 32542               |
| [OVER] | GO:0030246 | carbohydrate binding                                                                                  | MOLECULAR_FUNCTION | 0.000237514126308    | 1.68543750927969E-06 | 13      | 156          | 479            | 32395               |
| [OVER] | GO:0051287 | NAD binding                                                                                           | MOLECULAR_FUNCTION | 0.000264098183462    | 1.89513925303375E-06 | 7       | 33           | 485            | 32518               |
| [OVER] | GO:0004553 | hydrolase activity, hydrolyzing O-glycosyl compounds                                                  | MOLECULAR_FUNCTION | 0.000290760654919    | 2.10964914667664E-06 | 12      | 134          | 480            | 32417               |
| [OVER] | GO:0008235 | metalloexopeptidase activity                                                                          | MOLECULAR_FUNCTION | 0.000364373992658    | 2.67281193784983E-06 | 7       | 35           | 485            | 32516               |
| [OVER] | GO:0004867 | serine-type endopeptidase inhibitor activity                                                          | MOLECULAR_FUNCTION | 0.000391710169554    | 2.96259175874413E-06 | 8       | 52           | 484            | 32499               |
| [OVER] | GO:0033961 | cis-stilbene-oxide hydrolase activity                                                                 | MOLECULAR_FUNCTION | 0.0004242665152      | 3.28128304691423E-06 | 3       | 0            | 489            | 32551               |
| [OVER] | GO:0004180 | carboxypeptidase activity                                                                             | MOLECULAR_FUNCTION | 0.000814490903444    | 6.75387130905724E-06 | 7       | 41           | 485            | 32510               |
| [OVER] | GO:0007586 | digestion                                                                                             | BIOLOGICAL_PROCESS | 0.00171539138721     | 1.53184368814839E-05 | 5       | 17           | 487            | 32534               |
| [OVER] | GO:0042445 | hormone metabolic process                                                                             | BIOLOGICAL_PROCESS | 0.002124386056348    | 1.94789026056455E-05 | 6       | 32           | 486            | 32519               |
| [OVER] | GO:0044262 | cellular carbohydrate metabolic process                                                               | BIOLOGICAL_PROCESS | 0.002323360059406    | 2.14885797234208E-05 | 8       | 70           | 484            | 32481               |
| [OVER] | GO:0005581 | collagen trimer                                                                                       | CELLULAR_COMPONENT | 0.002807044995774    | 2.64097679398296E-05 | 6       | 34           | 486            | 32517               |
| [OVER] | GO:0005783 | endoplasmic reticulum                                                                                 | CELLULAR_COMPONENT | 0.003088027610687    | 3.00382210575479E-05 | 16      | 305          | 476            | 32246               |
| [OVER] | GO:0004471 | malate dehydrogenase (decarboxylating) (NAD+) activity                                                | MOLECULAR_FUNCTION | 0.003179914378633    | 3.20886775608155E-05 | 3       | 2            | 489            | 32549               |
| [OVER] | GO:0031012 | extracellular matrix                                                                                  | CELLULAR_COMPONENT | 0.003179914378633    | 3.21997389639884E-05 | 13      | 209          | 479            | 32342               |
| [OVER] | GO:0005506 | iron ion binding                                                                                      | MOLECULAR_FUNCTION | 0.003324288656178    | 3.41917745691992E-05 | 12      | 180          | 480            | 32371               |
| [OVER] | GO:0005811 | lipid droplet                                                                                         | CELLULAR_COMPONENT | 0.003538963144584    | 3.69641342641158E-05 | 8       | 76           | 484            | 32475               |
| [OVER] | GO:0004035 | alkaline phosphatase activity                                                                         | MOLECULAR_FUNCTION | 0.005365666600042    | 5.81829578700181E-05 | 4       | 11           | 488            | 32540               |
| [OVER] | GO:0016705 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | MOLECULAR_FUNCTION | 0.005909892448019    | 6.45555147008983E-05 | 12      | 193          | 480            | 32358               |
| [OVER] | GO:0006108 | malate metabolic process                                                                              | BIOLOGICAL_PROCESS | 0.009435308141185    | 0.000109835352305    | 3       | 4            | 489            | 32547               |
| [OVER] | GO:0004497 | monooxygenase activity                                                                                | MOLECULAR_FUNCTION | 0.010827399004446    | 0.000130356980519    | 10      | 147          | 482            | 32404               |
| [OVER] | GO:0006457 | protein folding                                                                                       | BIOLOGICAL_PROCESS | 0.014579595538693    | 0.000180181574589    | 8       | 97           | 484            | 32454               |
| [OVER] | GO:0016811 | hydrolase activity, acting on carbon-nitrogen (but not peptide) bonds, in linear amides               | MOLECULAR_FUNCTION | 0.016010748482324    | 0.000201698155016    | 6       | 51           | 486            | 32500               |
| [OVER] | GO:0007374 | posterior midgut invagination                                                                         | BIOLOGICAL_PROCESS | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0010765 | positive regulation of sodium ion transport                                                           | BIOLOGICAL_PROCESS | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0006533 | aspartate catabolic process                                                                           | BIOLOGICAL_PROCESS | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0007320 | insemination                                                                                          | BIOLOGICAL_PROCESS | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0008445 | D-aspartate oxidase activity                                                                          | MOLECULAR_FUNCTION | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0008756 | o-succinylbenzoate-CoA ligase activity                                                                | MOLECULAR_FUNCTION | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0030478 | actin cap                                                                                             | CELLULAR_COMPONENT | 0.016323702557707    | 0.000221258924797    | 2       | 0            | 490            | 32551               |
| [OVER] | GO:0005759 | mitochondrial matrix                                                                                  | CELLULAR_COMPONENT | 0.016561432819864    | 0.000225801707239    | 7       | 75           | 485            | 32476               |
| [OVER] | GO:0009396 | folic acid-containing compound biosynthetic process                                                   | BIOLOGICAL_PROCESS | 0.018798719042997    | 0.000257804152081    | 3       | 6            | 489            | 32545               |
| [OVER] | GO:0005777 | peroxisome                                                                                            | CELLULAR_COMPONENT | 0.020800650029868    | 0.000293550873488    | 6       | 55           | 486            | 32496               |
| [OVER] | GO:0007155 | cell adhesion                                                                                         | BIOLOGICAL_PROCESS | 0.021233797992279    | 0.000303049740123    | 15      | 338          | 477            | 32213               |
| [OVER] | GO:0005829 | cytosol                                                                                               | CELLULAR_COMPONENT | 0.023561403482919    | 0.000338148032764    | 13      | 268          | 479            | 32283               |
| [OVER] | GO:0019318 | hexose metabolic process                                                                              | BIOLOGICAL_PROCESS | 0.024276059435578    | 0.000350340197563    | 6       | 57           | 486            | 32494               |
| [OVER] | GO:0042737 | drug catabolic process                                                                                | BIOLOGICAL_PROCESS | 0.028449787881907    | 0.000428720292591    | 7       | 84           | 485            | 32467               |
| [OVER] | GO:0070938 | contractile ring                                                                                      | CELLULAR_COMPONENT | 0.032522118198718    | 0.000495273845954    | 3       | 8            | 489            | 32543               |
| [OVER] | GO:0007375 | anterior midgut invagination                                                                          | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0038095 | Fc-epsilon receptor signaling pathway                                                                 | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0019343 | cysteine biosynthetic process via cystathionine                                                       | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0035074 | pupation                                                                                              | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0006535 | cysteine biosynthetic process from serine                                                             | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0019395 | fatty acid oxidation                                                                                  | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000649366597563    | 4       | 23           | 488            | 32528               |
| [OVER] | GO:0071346 | cellular response to interferon-gamma                                                                 | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0009251 | glucan catabolic process                                                                              | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0035357 | peroxisome proliferator activated receptor signaling pathway                                          | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0031670 | cellular response to nutrient                                                                         | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0009062 | fatty acid catabolic process                                                                          | BIOLOGICAL_PROCESS | 0.038698500478984    | 0.000649366597563    | 4       | 23           | 488            | 32528               |
| [OVER] | GO:0004467 | long-chain fatty acid-CoA ligase activity                                                             | MOLECULAR_FUNCTION | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0004122 | cystathionine beta-synthase activity                                                                  | MOLECULAR_FUNCTION | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0020037 | heme binding                                                                                          | MOLECULAR_FUNCTION | 0.038698500478984    | 0.000650510065875    | 10      | 182          | 482            | 32369               |
| [OVER] | GO:0003878 | ATP citrate synthase activity                                                                         | MOLECULAR_FUNCTION | 0.038698500478984    | 0.000657214208422    | 2       | 1            | 490            | 32550               |
| [OVER] | GO:0005201 | extracellular matrix structural constituent                                                           | MOLECULAR_FUNCTION | 0.041804794296862    | 0.000719967753797    | 6       | 66           | 486            | 32485               |
| [OVER] | GO:0010038 | response to metal ion                                                                                 | BIOLOGICAL_PROCESS | 0.045658181827508    | 0.000797252576959    | 5       | 44           | 487            | 32507               |
| [OVER] | GO:0050680 | negative regulation of epithelial cell proliferation                                                  | BIOLOGICAL_PROCESS | 0.047223115099143    | 0.000839639185705    | 3       | 10           | 489            | 32541               |
| [OVER] | GO:0070633 | transepithelial transport                                                                             | BIOLOGICAL_PROCESS | 0.047223115099143    | 0.000839639185705    | 3       | 10           | 489            | 32541               |
| [OVER] | GO:0001889 | liver development                                                                                     | BIOLOGICAL_PROCESS | 0.047637349641279    | 0.00085839906067     | 4       | 25           | 488            | 32526               |
| [OVER] | GO:0044425 | membrane part                                                                                         | CELLULAR_COMPONENT | 0.047637349641279    | 0.000853779834716    | 54      | 2232         | 438            | 30319               |
| [OVER] | GO:0006875 | cellular metal ion homeostasis                                                                        | BIOLOGICAL_PROCESS | 0.048351471667672    | 0.000875122314508    | 5       | 45           | 487            | 32506               |
| [OVER] | GO:0006091 | generation of precursor metabolites and energy                                                        | BIOLOGICAL_PROCESS | 0.049640335352882    | 0.000902407627209    | 9       | 157          | 483            | 32394               |



## 11.2.4. Control **of the differential allele usage using Chom genomic reads**


The pipeline used to identify the condition-specific SNPs returned interresting results: some SNPs were found “homozygous” (in expression) in the control condition and “heterozygous” in the resistant conditions. 
We hypothesize that both alleles are present in the genome of the control population but that only one expressed itself (hence we see the locus as homozygous in the control). In the resistant, both alleles are being expressed. To test this hypothesis we used genomic reads from the original population.
The other explanation is that new alleles appeared in the resistant group. However, this does not seems plausible as the experiment was performed on the same population (ther
Another explanation is that when the population was sampled for RNA extraction some alleles were lost (bottleneck effect). However, more than 40 individuals were sampled for each condition. 

The genomic reads of the original population (composed of susceptible and resistant individuals) can be found trimmed on the Mendel server here: 
 /mnt/HD4/Genoma/1-Trimmomatic/
 

    # Get the sequenced of the transcripts containing condition-specific SNPs. To do this we need a file with the complete header (without the prompt) of the contigs of interest and the fasta file which had the sequences
    
    ## Get the contig names of the cs Transcripts with the prompt (>) + sorting
    sed 's/^/>/' ID-SNPs-Pval01-Cov20-diff20-uniq-sorted.txt | LANG=en_EN sort > ID-SNPs-Pval01-Cov20-diff20-uniq-sorted-with-prompt.txt
    
    ## Get the headers of the contig sequences
    grep "^>" ../../Transcriptome-and-annot/Trinity.longest.clusters.fasta | LANG=en_EN sort > headers-Trinity.longest.clusters.sorted.txt
    
    # Join the list of partial header wit hthe one with the complete headers (we get a file with the complete headers of the sequences of interest)
    LANG=en_EN join -t ' ' ID-SNPs-Pval01-Cov20-diff20-uniq-sorted-with-prompt.txt headers-Trinity.longest.clusters.sorted.txt > headers-SNPs-Pval01-Cov20-diff20.txt
    
    sed 's/>//' headers-SNPs-Pval01-Cov20-diff20.txt > headers-SNPs-Pval01-Cov20-diff20-2.txt
    
    ~/Programs/seqtk/seqtk subseq ../../Transcriptome-and-annot/Trinity.longest.clusters.fasta headers-Trinity.longest.clusters.sorted.txt > seqs-SNPs-Pval01-Cov20-diff20.fasta
    
    rsync -avz -P seqs-SNPs-Pval01-Cov20-diff20.fasta sophie@143.107.247.2:/mnt/HD6/Sophie/resistance-chom/
    
    /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2-build -f seqs-SNPs-Pval01-Cov20-diff20.fasta seqs-SNPs-Pval01-Cov20-diff20.BWT2
    
    /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x seqs-SNPs-Pval01-Cov20-diff20.BWT2 -1 /mnt/HD4/Genoma/1-Trimmomatic/Chom-R1-paired.fastq -2 /mnt/HD4/Genoma/1-Trimmomatic/Chom-R2-paired.fastq -S bowtie2-alignments-genomicReads-to-csSNP-transcripts.sam
    154067376 reads; of these:
      154067376 (100.00%) were paired; of these:
        153752513 (99.80%) aligned concordantly 0 times
        301240 (0.20%) aligned concordantly exactly 1 time
        13623 (0.01%) aligned concordantly >1 times
        ----
        153752513 pairs aligned concordantly 0 times; of these:
          1987 (0.00%) aligned discordantly 1 time
        ----
        153750526 pairs aligned 0 times concordantly or discordantly; of these:
          307501052 mates make up the pairs; of these:
            307460309 (99.99%) aligned 0 times
            38352 (0.01%) aligned exactly 1 time
            2391 (0.00%) aligned >1 times
    0.22% overall alignment rate

 




## **12. 1-to-1 Orthologs search**

Version: 2014 
Paper: Yang Y, Smith SA. 2014. Orthology inference in non-model organisms using transcriptomes and low coverage genomes: Improving accuracy and matrix occupancy for phylogenomics. Molecular Biology and Evolution 31(11): 3081–3092. 
Homepage: [https://bitbucket.org/yangya/phylogenomic_dataset_construction](https://bitbucket.org/yangya/phylogenomic_dataset_construction) 

**12.****1** **Getting the sequences**

***Cochliomyia hominivorax*** **(Chom.tra.fasta):** transcriptome assembled on this project; getting the longest isoform:


    perl /mnt/HD2/seq_an/trinityrnaseq_r20140717/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /mnt/HD5/Chom-resistance/Trinity.Trinity.fasta >/mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Chom.tra.fasta

***Glossin******g morsitans*** **(Gmor.tra.fasta):** Initiative IGG, Attardo GM, Abila PP *et al.* (2014) Genome Sequence of the Tsetse Fly (Glossina morsitans): Vector of African Trypanosomiasis. *Science (New York, NY)*, **344**, 380–386.

            Downaloded from VectorBase: https://www.vectorbase.org/organisms/glossina-morsitans
            Version: 1.8
            21 Nov 2017

***Musca domestica*** ********(Mdom.tra.fasta):** Scott JG, Warren WC, Beukeboom LW *et al.* (2014) Genome of the house fly, Musca domestica L., a global vector of diseases with adaptations to a septic environment.

            Downloaded from VectorBase:  https://www.vectorbase.org/organisms/musca-domestica
            Version: 1.3
            21 Nov 2017

***Stomoxys calcitrans*** ********(Scal.tra.fasta):** ******Not published.

            Downloaded from VectorBase:  https://www.vectorbase.org/organisms/stomoxys-calcitrans
            Version: 1.3
            21 Nov 2017

***Lucil******la cuprina*** **(Lcup.cds.fasta):** ****Anstead C a., Korhonen PK, Young ND *et al.* (2015b) Lucilia cuprina genome unlocks parasitic fly biology to underpin future interventions. *Nature Communications*, **6**, 7344.

            Downaloded from The sheep Blowfly Project FTP: ftp://ftp.hgsc.bcm.edu/I5K-pilot/Sheep_blowfly/maker_annotation/version_0.5.3/\
            21 Nov 2017

***Chrysomya megacephala*** ********(Cmeg.tra.fasta):** ******Xiaoyun W, Mei X, Chaoliang L *e*t *al* (2015)The developmental transcriptome of the synanthropic fly *Chrysomya megacephala* and insights into olfactory proteins. *BMC Genomics*, 16(1): 20.

            Downloaded from NCBI (via SRA Toolkit): https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR1663113
            Run: SRR1663113 
            Experiment: SRX768308
            Biosample: SAMN03200162 (SRS752515)
            SRA  study: SRP050024
            Condition: Third instar feeding larvae
            17 Nov 2015 - Assembled by Marina

getting the longest isoform:


    perl /mnt/HD2/seq_an/trinityrnaseq_r20140717/util/misc/get_longest_isoform_seq_per_trinity_gene.pl /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Cmeg.fasta >/mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Cmeg.tra.fasta

**12.2** ******Change contigs/transcripts name to species@**


    sed -i 's|>|>Chom@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Chom.tra.fasta


    sed -i 's|>|>Cmeg@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Cmeg.tra.fasta


    sed -i 's|>|>Gmor@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Gmor.tra.fasta


    sed -i 's|>|>Lcup@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Lcup.cds.fasta


    sed -i 's|>|>Mdom@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Mdom.tra.fasta


    sed -i 's|>|>Scal@|' /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Scal.tra.fasta

**12.3** ******Getting CDS using TransDecoder**

The step is not necessary for *L. cuprina*, since we already have the cds sequences.


    mkdir /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder | cd /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Chom.tra.fasta


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Cmeg.tra.fasta


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Gmor.tra.fasta


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Mdom.tra.fasta


    TransDecoder.LongOrfs -t /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Scal.tra.fasta

Summary



**12.4 Correcting contigs/transcripts names from TransDecoder**


    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/fix_names_from_transdecoder.py /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Chom.tra.fasta.transdecoder_dir/ /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Chom.tra.fasta.transdecoder_dir/


    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/fix_names_from_transdecoder.py /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Cmeg.tra.fasta.transdecoder_dir/ /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Cmeg.tra.fasta.transdecoder_dir/ 

 

    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/fix_names_from_transdecoder.py /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Gmor.tra.fasta.transdecoder_dir/ /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Gmor.tra.fasta.transdecoder_dir/

 

    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/fix_names_from_transdecoder.py /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Mdom.tra.fasta.transdecoder_dir/ /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Mdom.tra.fasta.transdecoder_dir/

 

    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/fix_names_from_transdecoder.py /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Scal.tra.fasta.transdecoder_dir/ /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Scal.tra.fasta.transdecoder_dir/

 
**12.5 Clustering contigs/transcripts using cd-hit-est**


    mkdir /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Chom.tra.fasta.transdecoder_dir/longest_orfs.fa.cds -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Chom.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Cmeg.tra.fasta.transdecoder_dir/longest_orfs.fa.cds -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Cmeg.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Gmor.tra.fasta.transdecoder_dir/longest_orfs.fa.cds -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Gmor.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Mdom.tra.fasta.transdecoder_dir/longest_orfs.fa.cds -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Mdom.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/2-TransDecoder/Scal.tra.fasta.transdecoder_dir/longest_orfs.fa.cds -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Scal.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 

    cd-hit-est -i /mnt/HD5/Chom-resistance/Orthologs/1-Sequences/Lcup.cds.fasta -o /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/Lcup.fa.cds.cd-hit-est -c 0.99 -n 10 -r 0 -T 15

 
 Chom.fa.cds.cd-hit-est.clstr:15519
Cmeg.fa.cds.cd-hit-est.clstr:18466
Gmor.fa.cds.cd-hit-est.clstr:16407
Lcup.fa.cds.cd-hit-est.clstr:18396
Mdom.fa.cds.cd-hit-est.clstr:25422
Scal.fa.cds.cd-hit-est.clstr:26076
 
**12.6** ******All-by-all-blast/cut sequences by alignment**
 

    mkdir /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all

 

    cat /mnt/HD5/Chom-resistance/Orthologs/3-cd-hit-est/*.cd-hit-est >/mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa

 

    makeblastdb -in /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa -parse_seqids -dbtype nucl -out /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa

 

    blastn -task blastn -db /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa -query /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa -evalue 0.0001 -num_threads 20 -max_target_seqs 1000 -out /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.raw.blast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'


    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/all-by-all_blast_to_mcl.py /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.raw.blast 0.4

**12.7** ******Markov clustering (mcl)**


    mkdir /mnt/HD5/Chom-resistance/Orthologs/5-mcl


    mcl /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.raw.blast.hit-frac0.4.minusLogEvalue --abc -te 10 -tf 'gq(5)' -I 1.4 -o /mnt/HD5/Chom-resistance/Orthologs/5-mcl/hit-frac0.4_I1.4_e5


    python /mnt/HD2/seq_an/yangya-phylogenomic_dataset_construction-61c9ee3932b2/write_fasta_files_from_mcl.py /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa /mnt/HD5/Chom-resistance/Orthologs/5-mcl/hit-frac0.4_I1.4_e5 6 /mnt/HD5/Chom-resistance/Orthologs/5-mcl/


    mkdir /mnt/HD5/Chom-resistance/Orthologs/5-mcl/MCL_trees


    mv /mnt/HD5/Chom-resistance/Orthologs/5-mcl/clus* /mnt/HD5/Chom-resistance/Orthologs/5-mcl/MCL_trees

**12.8** ******Tree constructing and sequences extraction**


    cd /mnt/HD5/Chom-resistance/Orthologs/5-mcl/MCL_trees/


    python /mnt/HD2/seq_an/yangya_novo2/fasta_to_tree.py 10 dna


    mkdir /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs


    python /mnt/HD2/seq_an/yangya_novo2/filter_one-to-one_orthologs.py /mnt/HD5/Chom-resistance/Orthologs/5-mcl/MCL_trees 6 /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs


            5725 files read, 3180 written to /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs/


    python /mnt/HD2/seq_an/yangya_novo2/write_fasta_files_from_trees.py /mnt/HD5/Chom-resistance/Orthologs/4-All-by-all/all.fa /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs .tre /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs/Sequences


    cd /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs/Sequences


    perl /mnt/HD4/Scripts/get_sequences_from_multiple_files.pl -s Chom,Cmeg,Gmor,Lcup,Mdom,Scal -e .fa -o /mnt/HD5/Chom-resistance/Orthologs/6-Orthologs/Sequences

It was found 3180 1-to-1 orthologos among the species *Co. hominivorax, Ch. megacephala, G. morsitans, L. cuprina, M. domestica* and  *S. calcitrans*.


# 13. New annotation

The reduced transcriptome (longest isoform + clustering) was annotated using FunctionAnnotator on the 16/09/2019 using the following input parameters:

ANALYSIS PARAMETERS:
1. Select sequence type:

-  Eukaryotic

2. Select analysis module:

-  Best hit (in NCBI-nr databsae),
-       & Taxonomic distribution for best hits, 
-       & GO function annotation (Blast2GO)
-  Enzyme prediction (PRIAM database)
-  Domain region identification (Domain finder)
-  Transmembrane domain prediction (TMHMM)
-  Subcellular localization prediction (Psort)
-  Secretory protein prediction (SignalP)

Here are some basic information outputed by FunctionAnnotator, which could be interesting:

| Job ID            | 15683825642878                 |
| ----------------- | ------------------------------ |
| Fasta file        | Trinity.longest.clusters.fasta |
| File size         | 29,897,208 bytes               |
| Number of Entries | 33,038 entries                 |
| Uploaded on       | Fri, 13 Sep 19 21:49:24 +0800  |

Filtered:

| aaseq with length <= 66 | 12,694 |

![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730559896_image.png)

![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730573854_image.png)




![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730936874_AnnotationTypes.png)


**Enzyme identification**
There were **1089** entries identified to have at least one enzyme hit.

**Domain identification**
There were **9010** entries identified to have at least one domain with domain coverage > 50%. 

**Transmembrane domain identification**

| # of predicted TM (TransMembrane) domain | Count |
| ---------------------------------------- | ----- |
| Entries with only 1 TM Domain            | 2749  |
| Entries with multiple TM Domains         | 1998  |
| all                                      | 4747  |

There were **4747** entries identified to have at least one transmembrane domain.

**Signal peptide identification**
There were **1115** entries predicted to have signal peptides.

**hits against NCBI non-redundant protein database**
There are 16541 entries found hit to NCBI-nr database.

**Mapping to Gene Ontoloty**
There are **11055** entries found mapping to Gene Ontology terms by blast2go.
9639 entries found mapping to Biological Process terms on this level
7442 entries found mapping to Cellular Component terms on this level
9688 entries found mapping to Molecular Function terms on this level


![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730697288_image.png)

![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730712156_image.png)



![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730728550_image.png)


Taxomonic distribution of the annotations:

- Species
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730763294_image.png)

- Genus
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730815183_image.png)

- Family
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730833627_image.png)



- Order
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568730856073_image.png)






Gisele’s notes
I used a different approach to annoted the whole transcriptome (located in Darwin: /Volumes/HD3/analyses/Tatiana/Chom-resistance/4-assembly/Trinity.Trinity.fasta) and the DEs. I did not erased the last annotation protocol just in case we need to go back for some reason.

Both whole transciptome and DE genes were annotated with a web-based tool called FunctionAnnotator (http://fa.cgu.edu.tw/index.php). The results were in the following tables.


https://www.dropbox.com/s/7n1hbxqvu21elum/Transcriptome_annotation.txt?dl=0



https://www.dropbox.com/s/b7gengud69628jx/DE_annotation.txt?dl=0




|  |



# 14. Selection test

We used only the CDS to search for orthologs among 7 fly species. Thus, the first thing to do was to eliminate stop codons. You just need to enter your directory and run the following script:


[x] perl /mnt/HD4/Scripts/eliminate_stop_codons_seqs.pl -e .cds

Now, we can simple run translatorX that translate CDS to protein, do the aminoacid alignment with muscle and for last back-translate the alignment to fasta. Also, all alignments are “corrected” with gblocks leaving all gaps in the alignments. I made a script to run tranlatorX in all disered seqeunces, you just need to enter your directory and hit:


[x] perl /mnt/HD4/Scripts/translatorX.pl -c .fa -g -b5=a

Then, we need to put all sequences name just with three letters ( for more information, search the features of the phylip format). After this correction, we transform the fasta alignment to phylip.


[x] for i in *.align-cln; do sed -i ‘/>Chom.*/>Chom/’ $i
[x] for i in *.align-cln; do sed -i ‘/>Scal.*/>Scal/’ $i
[x] for i in *.align-cln; do sed -i ‘/>Lcup.*/>Lcup/’ $i
[x] for i in *.align-cln; do sed -i ‘/>Gmor.*/>GMor/’ $i
[x] for i in *.align-cln; do sed -i ‘/>Mdom.*/>Mdom/’ $i
[x] for i in *.align-cln; do sed -i ‘/>Preg.*/>Preg/’ $
[x] for i in *.align-cln; do sed -i ‘/>Ccap*/>Ccap/’ $i


[x] perl /mnt/HD4/Scripts/fasta_to_phylip.pl -e .align-cln

With all alignments in phylip format, we can run codeml for each ortholog with a perl script. This script is not running very well (I did not had the courage to improve it) and it is not working to change the codeml.ctl file. SO, we need to change manually for the parameters we want. This time, I am trying to use the NSsite features, that means that omega can vary among sites. Also, I followed the suggestions in the PAML manual. For the whole tree, I used 7 models. We need to set all of them to model =0 and the NSsites according to what we want:


| Model               | NSsites |
| ------------------- | ------- |
| Model 0 (one ratio) | 0       |
| M1a (neutral)       | 1       |
| M2a (selection)     | 2       |
| M2a_ref             | 22      |
| M3 (discrete)       | 3       |
| M7 (beta)           | 7       |
| M8 (beta&ω)         | 8       |

Make a new directory for each model you want to run and put all sequences, tree file and codeml.ctl  in each directory. Then, make the chanes you need in the codeml.ctl and run the codeml script:


[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M0.log 2> tree_M0.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M1a.log 2> tree_M1a.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M2a.log 2> tree_M2a.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M2a_ref.log 2> tree_M2a_ref.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M3.log 2> tree_M3.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M7.log 2> tree_M7.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > tree_M8.log 2> tree_M8.err

For model= 2 taging *C. hominivorax* we can use only three models (I am almost pretty sure that’s it!):

| Model               | NSsites |
| ------------------- | ------- |
| Model 0 (one ratio) | 0       |
| M2a (selection)     | 2       |
| M3 (discrete)       | 3       |

Do the same protocol from before and run the codeml script. Don’t forget to tag the branch you want with a # followed by a number in the end of the species name. Also, tspecies name in the tree file must be exactly the same ass the alignment files.


[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > branch_M0.log 2> branch_M0.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > branch_M2a.log 2> branch_M2a.err
[x] perl /mnt/HD4/Scripts/codeml.pl -s .phy -t tree.tre > branch_M3.log 2> branch_M3.err


# 15. Comparing the DE genes to the condition-specific polymorphic genes

142 genes were found DE (analysis based on the reduced transcriptome)
614 transcripts had SNPs which were condition-specific.

Eye-balling the data it seems like the DE genes and the polymorphic gees are associated with the same processes. Could it be due to the fact that the same genes are found DE and polymorphic? 

**15.1. How many genes are DE and have condition-specific SNPs?**


    Quick Venn diagram (made online using Venny 2.1. http://bioinfogp.cnb.csic.es/tools/venny/):
![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568735788926_image.png)


Only 18 DE transcripts also have condition-specific SNPs. Who are they?


    join -t '     ' ID-SNPs-R1R2-P01-cov20-diff15-uniq.sorted.txt ID-DE-contigs.sorted | wc -l
    18
    join -t '     ' ID-SNPs-R1R2-P01-cov20-diff15-uniq.sorted.txt .ID-DE-contigs.sorted > ID-common-contigs-DE-and-SNPs.txt
    
    LANG=en_EN sort ID-common-contigs-DE-and-SNPs.txt > ID-common-contigs-DE-and-SNPs.sorted.txt
    
    # I manually added the header to the file of sorted common IDs
    
    join -t '        ' ID-common-contigs-DE-and-SNPs.sorted.txt ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt > common-contigs-DE-and-SNPs-annotated.txt

common-contigs-DE-and-SNPs-annotated.txt:

| #Component_ID | length | best_hit_to_nr                                                                                          | hit_length | E-value   | Bit_score   | GO_Biological_Process                                                                                                                      | GO_Cellular_Component                                       | GO_Molecular_Function                                                                                                                                                                                        | Enzyme             | species                  | genus      | family        |
| ------------- | ------ | ------------------------------------------------------------------------------------------------------- | ---------- | --------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------ | ------------------------ | ---------- | ------------- |
| c11348_g2_i1  | 444    | -                                                                                                       | -          | -         | -           | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | -                        | -          | -             |
| c11548_g1_i1  | 429    | gi|906469296|gb|KNC30110.1| hypothetical protein FF38_09557                                             | 72         | 1.39E-15  | 92.155133   | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c12633_g1_i1  | 689    | gi|906466050|gb|KNC27577.1| hypothetical protein FF38_05113                                             | 150        | 6.12E-62  | 238.487691  | GO:0007165 signal transduction | GO:0006508 proteolysis | GO:0006412 translation | GO:0051301 cell division                                | GO:0005840 ribosome                                         | GO:0003735 structural constituent of ribosome | GO:0004197 cysteine-type endopeptidase activity | GO:0000166 nucleotide binding                                                                              | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c12633_g1_i5  | 607    | gi|906466050|gb|KNC27577.1| hypothetical protein FF38_05113                                             | 181        | 1.12E-74  | 275.752504  | GO:0006810 transport | GO:0006412 translation | GO:0051301 cell division                                                                   | GO:0005840 ribosome                                         | GO:0003735 structural constituent of ribosome | GO:0003824 catalytic activity | GO:0000166 nucleotide binding                                                                                                | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c12908_g1_i1  | 1830   | gi|659106067|gb|AID61420.1| cytochrome P450                                                             | 517        | 0         | 1178.288095 | GO:0055114 oxidation-reduction process                                                                                                     | -                                                           | GO:0020037 heme binding | GO:0016705 oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen | GO:0004497 monooxygenase activity | GO:0005506 iron ion binding | -                  | 145453|Calliphora stygia | Calliphora | Calliphoridae |
| c13331_g1_i1  | 1024   | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 251        | 8.41E-71  | 265.754627  | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13339_g1_i1  | 989    | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 290        | 2.91E-99  | 354.826619  | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13437_g3_i1  | 1302   | gi|557780752|ref|XP_005189997.1| PREDICTED: uncharacterized protein LOC101896308                        | 278        | 2.18E-150 | 524.790522  | -                                                                                                                                          | -                                                           | GO:0008289 lipid binding                                                                                                                                                                                     | -                  | 7370|Musca domestica     | Musca      | Muscidae      |
| c13617_g1_i1  | 678    | gi|906464505|gb|KNC26324.1| hypothetical protein FF38_00857, partial                                    | 175        | 5.47E-47  | 193.042797  | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13628_g1_i1  | 926    | gi|906464455|gb|KNC26277.1| hypothetical protein FF38_01137                                             | 290        | 1.11E-121 | 429.356245  | -                                                                                                                                          | -                                                           | GO:0030246 carbohydrate binding                                                                                                                                                                              | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13820_g1_i1  | 915    | gi|906473923|gb|KNC33625.1| hypothetical protein FF38_08769                                             | 249        | 3.83E-143 | 500.704728  | GO:0010765 positive regulation of sodium ion transport | GO:0006508 proteolysis | GO:0007586 digestion                                     | GO:0005615 extracellular space | GO:0005886 plasma membrane | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                | 3.4.21.4           | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13820_g1_i7  | 1130   | gi|906473923|gb|KNC33625.1| hypothetical protein FF38_08769                                             | 251        | 9.68E-146 | 509.339258  | GO:0010765 positive regulation of sodium ion transport | GO:0006508 proteolysis | GO:0007586 digestion                                     | GO:0005615 extracellular space | GO:0005886 plasma membrane | GO:0004252 serine-type endopeptidase activity                                                                                                                                                                | 3.4.21.4           | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c13993_g1_i2  | 2146   | gi|906463872|gb|KNC25845.1| hypothetical protein FF38_08124                                             | 629        | 0         | 1379.154525 | GO:0006879 cellular iron ion homeostasis | GO:0006826 iron ion transport                                                                   | GO:0005576 extracellular region                             | GO:0008199 ferric iron binding                                                                                                                                                                               | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c14119_g1_i1  | 1763   | gi|906459360|gb|KNC22368.1| hypothetical protein FF38_08825                                             | 538        | 0         | 1134.660997 | GO:0009851 auxin biosynthetic process | GO:0009695 jasmonic acid biosynthetic process | GO:0001676 long-chain fatty acid metabolic process | GO:0005777 peroxisome                                       | GO:0004321 fatty-acyl-CoA synthase activity | GO:0004467 long-chain fatty acid-CoA ligase activity | GO:0016207 4-coumarate-CoA ligase activity | GO:0008756 o-succinylbenzoate-CoA ligase activity          | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c14131_g2_i3  | 2256   | gi|755851246|ref|XP_005175488.2| PREDICTED: bifunctional 3'-phosphoadenosine 5'-phosphosulfate synthase | 649        | 0         | 1433.233949 | GO:0016310 phosphorylation | GO:0000103 sulfate assimilation                                                                               | -                                                           | GO:0004781 sulfate adenylyltransferase (ATP) activity | GO:0005524 ATP binding | GO:0004020 adenylylsulfate kinase activity                                                                                  | 2.7.1.25 | 2.7.7.4 | 7370|Musca domestica     | Musca      | Muscidae      |
| c14257_g1_i1  | 937    | gi|906468463|gb|KNC29468.1| hypothetical protein FF38_13303                                             | 292        | 1.13E-118 | 419.358368  | -                                                                                                                                          | -                                                           | -                                                                                                                                                                                                            | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c14544_g1_i4  | 1906   | gi|906469326|gb|KNC30137.1| hypothetical protein FF38_12389                                             | 546        | 0         | 1182.378135 | GO:0006212 uracil catabolic process                                                                                                        | GO:0005737 cytoplasm                                        | GO:0004157 dihydropyrimidinase activity | GO:0051219 phosphoprotein binding                                                                                                                                  | 3.5.2.2            | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |
| c8201_g1_i1   | 506    | gi|2565392|gb|AAB81989.1| cuticle 1                                                                     | 118        | 2.52E-52  | 208.494061  | -                                                                                                                                          | -                                                           | GO:0042302 structural constituent of cuticle                                                                                                                                                                 | -                  | 7375|Lucilia cuprina     | Lucilia    | Calliphoridae |



**15.2. How many unique GO terms are in common between the DE set and the SNP set**


    # Getting the uniq GO terms associated with the DE transcripts
    grep -o "GO:[0-9]*" ../9-testing/Annotating-DE-contigs/DEcontigsAll-annotated.txt | LANG=en_EN sort | uniq > Uniq-GO-terms-of-DEGs.sorted.txt
    
    # Getting the unique GO terms associated with the transcripts containing condition-specific SNPs
    grep -o "GO:[0-9]*" ../11-Variant-analysis/ChomR1R2-Pval01-nonSyn-condSpe-Cov20-diff15-with-annotation.txt | LANG=en_EN sort | uniq > Uniq-GO-terms-of-polymorphic-genes.sorted.txt
    
    # Counting
    wc -l Uniq-GO-terms-of-*
      208 Uniq-GO-terms-of-DEGs.sorted.txt
     1191 Uniq-GO-terms-of-polymorphic-genes.sorted.txt

First observation: there are many more different GO terms in the “transcripts with SNPs” set then in the DE gene set. This is directly due to the higher number of transcripts in the SNP set (614 vs 142)

A quick Venn diagram shows that 63.5% of unique GO terms found in the DE set is also in the SNP set. Not sure this is very interesting as these GO termscould be very common.

![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568737578443_image.png)


Which are these common GO terms? 


    # getting the 132 common GO terms
    
    join -t '        ' Uniq-GO-terms-of-DEGs.sorted.txt Uniq-GO-terms-of-polymorphic-genes.sorted.txt > Common-GO-terms-DE-SNP.txt
    
    


    # Getting the description of the GO code names from the annotated transcriptome
    
    # Sanity Check: How many GO terms are present in the annotation?
    grep -o "GO:[0-9]*" ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt | wc -l
    72244
    
    # Getting the unique GO terms plus their description + adding tab + sorting
    
    grep -o "GO:[0-9]*.*     " ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt | tr '\|' '\n' | tr '      ' '\n' | grep -o "GO:[0-9]* .*" | sed 's/ /    /' | LANG=en_EN sort | sed 's/ $//' | uniq > ../Transcriptome-and-annot/fareport_15683825642878/GO-terms-with-descriptions-uniq-sorted.txt
    
    # Sanity Check before uniq
    grep -o "GO:[0-9]*.*     " ../Transcriptome-and-annot/fareport_15683825642878/AnnotationTable_trinity_longest_clusters-simplified.sorted.txt | tr '\|' '\n' | tr '      ' '\n' | grep -o "GO:[0-9]* .*" | sed 's/ /    /' | LANG=en_EN sort | sed 's/ $//' | wc -l
    72244


    # Common GO terms with description
    join -t '        ' Common-GO-terms-DE-SNP-sorted.txt ../Transcriptome-and-annot/fareport_15683825642878/GO-terms-with-descriptions-uniq-sorted.txt > Common-GO-terms-DE-SNP-with-description.txt
    
     wc -l Common-GO-terms-DE-SNP-with-description.txt
    132 Common-GO-terms-DE-SNP-with-description.txt


| GO:0000103 | sulfate assimilation                                                                                  |
| ---------- | ----------------------------------------------------------------------------------------------------- |
| GO:0000166 | nucleotide binding                                                                                    |
| GO:0001522 | pseudouridine synthesis                                                                               |
| GO:0001676 | long-chain fatty acid metabolic process                                                               |
| GO:0003676 | nucleic acid binding                                                                                  |
| GO:0003677 | DNA binding                                                                                           |
| GO:0003723 | RNA binding                                                                                           |
| GO:0003735 | structural constituent of ribosome                                                                    |
| GO:0003824 | catalytic activity                                                                                    |
| GO:0003964 | RNA-directed DNA polymerase activity                                                                  |
| GO:0004020 | adenylylsulfate kinase activity                                                                       |
| GO:0004022 | alcohol dehydrogenase (NAD) activity                                                                  |
| GO:0004035 | alkaline phosphatase activity                                                                         |
| GO:0004157 | dihydropyrimidinase activity                                                                          |
| GO:0004181 | metallocarboxypeptidase activity                                                                      |
| GO:0004197 | cysteine-type endopeptidase activity                                                                  |
| GO:0004222 | metalloendopeptidase activity                                                                         |
| GO:0004252 | serine-type endopeptidase activity                                                                    |
| GO:0004321 | fatty-acyl-CoA synthase activity                                                                      |
| GO:0004364 | glutathione transferase activity                                                                      |
| GO:0004467 | long-chain fatty acid-CoA ligase activity                                                             |
| GO:0004497 | monooxygenase activity                                                                                |
| GO:0004519 | endonuclease activity                                                                                 |
| GO:0004553 | hydrolase activity, hydrolyzing O-glycosyl compounds                                                  |
| GO:0004568 | chitinase activity                                                                                    |
| GO:0004781 | sulfate adenylyltransferase (ATP) activity                                                            |
| GO:0004867 | serine-type endopeptidase inhibitor activity                                                          |
| GO:0005212 | structural constituent of eye lens                                                                    |
| GO:0005215 | transporter activity                                                                                  |
| GO:0005506 | iron ion binding                                                                                      |
| GO:0005509 | calcium ion binding                                                                                   |
| GO:0005515 | protein binding                                                                                       |
| GO:0005518 | collagen binding                                                                                      |
| GO:0005524 | ATP binding                                                                                           |
| GO:0005549 | odorant binding                                                                                       |
| GO:0005576 | extracellular region                                                                                  |
| GO:0005578 | proteinaceous extracellular matrix                                                                    |
| GO:0005615 | extracellular space                                                                                   |
| GO:0005622 | intracellular                                                                                         |
| GO:0005634 | nucleus                                                                                               |
| GO:0005730 | nucleolus                                                                                             |
| GO:0005737 | cytoplasm                                                                                             |
| GO:0005739 | mitochondrion                                                                                         |
| GO:0005765 | lysosomal membrane                                                                                    |
| GO:0005777 | peroxisome                                                                                            |
| GO:0005789 | endoplasmic reticulum membrane                                                                        |
| GO:0005811 | lipid particle                                                                                        |
| GO:0005829 | cytosol                                                                                               |
| GO:0005840 | ribosome                                                                                              |
| GO:0005886 | plasma membrane                                                                                       |
| GO:0005975 | carbohydrate metabolic process                                                                        |
| GO:0006030 | chitin metabolic process                                                                              |
| GO:0006032 | chitin catabolic process                                                                              |
| GO:0006066 | alcohol metabolic process                                                                             |
| GO:0006200 | ATP catabolic process                                                                                 |
| GO:0006212 | uracil catabolic process                                                                              |
| GO:0006278 | RNA-dependent DNA replication                                                                         |
| GO:0006412 | translation                                                                                           |
| GO:0006508 | proteolysis                                                                                           |
| GO:0006559 | L-phenylalanine catabolic process                                                                     |
| GO:0006629 | lipid metabolic process                                                                               |
| GO:0006633 | fatty acid biosynthetic process                                                                       |
| GO:0006730 | one-carbon metabolic process                                                                          |
| GO:0006749 | glutathione metabolic process                                                                         |
| GO:0006807 | nitrogen compound metabolic process                                                                   |
| GO:0006810 | transport                                                                                             |
| GO:0006812 | cation transport                                                                                      |
| GO:0006826 | iron ion transport                                                                                    |
| GO:0006879 | cellular iron ion homeostasis                                                                         |
| GO:0006950 | response to stress                                                                                    |
| GO:0007165 | signal transduction                                                                                   |
| GO:0007586 | digestion                                                                                             |
| GO:0007596 | blood coagulation                                                                                     |
| GO:0007616 | long-term memory                                                                                      |
| GO:0008010 | structural constituent of chitin-based larval cuticle                                                 |
| GO:0008026 | ATP-dependent helicase activity                                                                       |
| GO:0008061 | chitin binding                                                                                        |
| GO:0008152 | metabolic process                                                                                     |
| GO:0008199 | ferric iron binding                                                                                   |
| GO:0008218 | bioluminescence                                                                                       |
| GO:0008237 | metallopeptidase activity                                                                             |
| GO:0008270 | zinc ion binding                                                                                      |
| GO:0008289 | lipid binding                                                                                         |
| GO:0008324 | cation transmembrane transporter activity                                                             |
| GO:0008484 | sulfuric ester hydrolase activity                                                                     |
| GO:0008756 | o-succinylbenzoate-CoA ligase activity                                                                |
| GO:0009695 | jasmonic acid biosynthetic process                                                                    |
| GO:0009790 | embryo development                                                                                    |
| GO:0009851 | auxin biosynthetic process                                                                            |
| GO:0009982 | pseudouridine synthase activity                                                                       |
| GO:0010171 | body morphogenesis                                                                                    |
| GO:0010765 | positive regulation of sodium ion transport                                                           |
| GO:0010951 | negative regulation of endopeptidase activity                                                         |
| GO:0016020 | membrane                                                                                              |
| GO:0016021 | integral to membrane                                                                                  |
| GO:0016207 | 4-coumarate-CoA ligase activity                                                                       |
| GO:0016310 | phosphorylation                                                                                       |
| GO:0016311 | dephosphorylation                                                                                     |
| GO:0016459 | myosin complex                                                                                        |
| GO:0016705 | oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen |
| GO:0016772 | transferase activity, transferring phosphorus-containing groups                                       |
| GO:0016787 | hydrolase activity                                                                                    |
| GO:0016874 | ligase activity                                                                                       |
| GO:0016887 | ATPase activity                                                                                       |
| GO:0018833 | DDT-dehydrochlorinase activity                                                                        |
| GO:0019013 | viral nucleocapsid                                                                                    |
| GO:0020037 | heme binding                                                                                          |
| GO:0022857 | transmembrane transporter activity                                                                    |
| GO:0030198 | extracellular matrix organization                                                                     |
| GO:0030246 | carbohydrate binding                                                                                  |
| GO:0031012 | extracellular matrix                                                                                  |
| GO:0031167 | rRNA methylation                                                                                      |
| GO:0031225 | anchored to membrane                                                                                  |
| GO:0035071 | salivary gland cell autophagic cell death                                                             |
| GO:0035725 | sodium ion transmembrane transport                                                                    |
| GO:0040003 | chitin-based cuticle development                                                                      |
| GO:0042048 | olfactory behavior                                                                                    |
| GO:0042302 | structural constituent of cuticle                                                                     |
| GO:0042742 | defense response to bacterium                                                                         |
| GO:0042802 | identical protein binding                                                                             |
| GO:0045087 | innate immune response                                                                                |
| GO:0046872 | metal ion binding                                                                                     |
| GO:0050660 | flavin adenine dinucleotide binding                                                                   |
| GO:0050829 | defense response to Gram-negative bacterium                                                           |
| GO:0051219 | phosphoprotein binding                                                                                |
| GO:0051260 | protein homooligomerization                                                                           |
| GO:0051301 | cell division                                                                                         |
| GO:0052689 | carboxylic ester hydrolase activity                                                                   |
| GO:0055085 | transmembrane transport                                                                               |
| GO:0055114 | oxidation-reduction process                                                                           |
| GO:0071391 | cellular response to estrogen stimulus                                                                |
| GO:0090305 | nucleic acid phosphodiester bond hydrolysis                                                           |


15.3. which GO terms are found over represented in both the DE transcript set (Up and down) and the condition-specific polymorphic transcript set?

On reduced GO term set

    cut -f 2 GO_overRep_DOWN-reduced.txt | grep -v "ID" | LANG=en_EN sort | uniq > GO-IDs-DOWN-uniq-sorted-reduced.txt
    
    cut -f 2 GO_overRep_UP-reduced.txt | grep -v "ID" | LANG=en_EN sort | uniq > GO-IDs-UP-uniq-sorted-reduced.txt
    
    cat GO-IDs-UP-uniq-sorted-reduced.txt GO-IDs-DOWN-uniq-sorted-reduced.txt | LANG=en_EN sort | uniq > GO-IDs-ALL-uniq-sorted-reduced.txt
    
    cut -f 2 GOs-enriched-in-polymorphic-transcripts-Pval01-Cov20-diff20-reduced.txt | grep -v "ID" | LANG=en_EN sort | uniq > GO-IDs-SNPset-uniq-sorted-reduced.txt


![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568928753158_image.png)


On the unreduced GO terms


    # Getting unique list of GO terms enriched in either the UP- or down (or both) transcripts.
    cut -f 2 GO_overRep_DOWN-unreduced.txt | grep -v "ID" | LANG=en_EN sort | uniq > GO-IDs-DOWN-uniq-sorted.txt
    
    cut -f 2 GO_overRep_UP-unreduced.txt | grep -v "ID" | LANG=en_EN sort | uniq > GO-IDs-UP-uniq-sorted.txt
    
    cat GO-IDs-UP-uniq-sorted.txt GO-IDs-DOWN-uniq-sorted.txt | LANG=en_EN sort | uniq > GO-IDs-ALL-uniq-sorted.txt



![](https://paper-attachments.dropbox.com/s_E4B4A65B9CCE1DECDD860280A6E312767E250AC530F2A15681F7BF9C29FF5B07_1568928338813_image.png)



# 16. Analysis on the Esterase E3

The Esterase E3 has been extensively studied for its role in insecticide resistance. In Co. hominivorax, the resistant individuals possess 2 point mutations, which confer resistance: G137D and W251S.



## 16.1. Frequencies of different haplotypes at positions 137 and 251

In the study, we expect to see a mixture of both alleles at both positions (G and D at position 137 and W and S at position 25Selvagem 1) in the control condition. However, in the resistant condition, we expect that all individual will carry the 137D and 251S mutations that confer resistance.

Contrary to our expectations, the 137 position does not appear in our condition-specific SNP analysis (we should see it homozygous for 137D in the resistant but not in the control). To look into this discrepancy, I have aligned genomic reads from the control population against the contig of the esterase E3 as well as the RNA-seq reads of the different conditions. With these informations we will check:
- that we have both alleles at positions 137 and 251 in the DNA of the control population
- that both alleles are expressed in the control population
- that only the mutations conferring resistance are expressed in the resistant population.

This analysis is being conducted in Mendel, in the folder ‘/mnt/HD6/Sophie/resistance-chom/esteraseE3’


    # Getting the contig corresponding to the Esterase E3
    
    grep -A 1 "c13624_g1" seqs-SNPs-Pval01-Cov20-diff20.fasta > esteraseE3/c13624_g1_i1_esteraseE3.fasta
    
    Take the reverse complement of the contig to have the transcript in 5'-3'. 
    
    The sequence of the Esterase E3 used for analysis is:
    
    >c13624_g1_i1_reversed len=2159
    TAATTTTAAACAAAACAGTTGTTAGTTGCTGTTGCTTTAAGCTTCTTGGGGATTTGAACTTTAATTATTCTTAAAATAAAACTGTGCGTTTGCGTGTTTTGCAGTGAACAAATTAAAAAAAGAAAAAGAAATAAAAACATAAACTTTATAATTTTCCAAGAAATCGTTTTAAATTGTATCTAGATACAAAGAAAGTTTATATATTTTATTGGCTTGTTTCTTTTTTTACTAGCAATTTCTAAAATTTTACCTCTTGTGCTTTAAAGCGAACAAAACGCCAAAATGAATTTCAACGTCAGTTTAATGGAGAAATTAAAATGGAAAATTAAATGTTTTGAAAATAAATTTCTAAATTATCGTTTGAGTACAAATGAAACGGCTGTGGCTGAAACAGAATATGGCAAAGTGAAAGGTATTAAACGTTTAACAGTTTATGATGATTCATATTATAGTTTTGAGGGTATACCATACGCCCAACCCCCCTTGGGTGAATTGAGATTTAAGGCACCCCAGCGACCAACACCTTGGGATGGTGTACGTGATTGTTGCAATAACAAAGATAAATCAGTGCAAGTTGATTTTATAACGGGTAAAACATGTGGTTCAGAGGATTGTTTATACTTGAGCGTCTATACGAATAATCTGACTCCAGAAACTAAACGTCCAGTTTTGGTATACATTCATGGTGGTGGCTTCGTTATCGGTGAAAATCATCGTGAATATTATGGACCTGATTATTTTATTAAGAAAGATGTTGTATTAATTACCATACAATATCGTTTGGGAGTTCTTGGTTTCTTGAGCTTAAATTCTGAAGAGCTTAATGTACCCGGTAATGCTGGCCTTAAGGATCAAGTTATGGCTTTACGTTGGATTAAAAATAATTGTGCCAATTTCGGTGGTAATCCTGATAATATCACTGTCTTTGGTGAGAGTGCTGGTGGAGCCTCTGCCCATTACATGATGTTAACTGAACAAACACGTGGTCTCTTCCATCGTGGTATTTTAATGTCTGGCAATGCTGTATGTCCTTCGGCCATTAGTCAAAATCAACATCGTGCTTATGCTATAGCTAAATTGACTGGTTATAAGGGTGAAAATAATGATAAGGATGTTTTGGAATTCTTAATGAAAGCTAAAGCTCATGATTTAATCAAATTGGAAGACAAAGTTTTGACACCCGAAGAACATGTGAATAAAGTGATGTTTGCTTTTGGTCCCACTGTGGAACCTTATCAGACTGCTGACTGTGTTTTACCTAAACATCCTAGAGAAATGGTGAAGACTGCTTGGGGTAATTCTATACCCACAATGATGGGTAATACTTCGTACGAGGGATTACTGTTTACACCGATTGTTAAACAAATGCCGGCACTTTTGAAAGAGTTGGAAACTTGTGCTAATTTTGTTCCTACTGAATTAGCCGATTCCGAACGCAGTTCTGCTGAAACCTTGGAATTGGGTGCCAAAATTAAAAAAGCCCATGTTACAGGAGAAACACCAACTAATGATAACTTTTTGGATCTTTGCTCACACTTTTACTTCTGGTTCCCCATGCATCGTCTACTCCAAATACGTTTCAAACATACATCTGGAACTCCTGTCTATTTATATCGTTTTGATTTTGATTCTGAAGAAATAATTAATCCTTATCGTATTATGCGTCATGGCCGTGGCGTTAAGGGTGTAAGTCATGCTGATGAATTAACATACCTTTTCTGGAATGCATTGGCCAAACGTTTGCCTAAAGAATCCCGAGAATACAAGACAATTGAACGTATGGTGGGTATATGGACTCAGTTTGCTACTACCGGTAATCCCTATAGCAATGAAATTGATGGCCTGGAAAATATTATTTGGGACCCTATCAAGAAATCTGATGAAGTCTATAAATGTTTAAATATAAGTGATGAATTGAAGATTATCGATGTGCCTGAAATGGAAAAAATTAAACAATGGGAATCACTGTATGAAAAACGCAAAGATTTATTTTAGAAATCTGATTTGTTTGTAACCAAAGATTGAGTGTGGAAGGAATTATGTTATTATTGTAATTTTATAATCTTTTTTTATGTATGTATATTTATATAATTATTAAAAAAAACAAATAAAAAATAGTGCATGTAATATTTTTTGCTTTTATTCTACTATAAAGTCGT

Nessa sequencia o codon 137 esta na posição 691 e o codon 251 esta na posição 1033.

**RNA alignments**


    [sophie@mendel esteraseE3]$ /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x c13624_g1_i1_esteraseE3_reversed.BWT2 -U ../trimmed-RNA-reads/1st-replicate/Control-1st-R1-trimmed.fastq -S c13624_g1_i1_esteraseE3-RNA-control-Rep1.sam
    11055926 reads; of these:
      11055926 (100.00%) were unpaired; of these:
        11053017 (99.97%) aligned 0 times
        2909 (0.03%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.03% overall alignment rate
    [sophie@mendel esteraseE3]$ /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x c13624_g1_i1_esteraseE3_reversed.BWT2 -U ../trimmed-RNA-reads/1st-replicate/Resistant-1st-R1-trimmed.fastq -S c13624_g1_i1_esteraseE3-RNA-resist-Rep1.sam
    11279633 reads; of these:
      11279633 (100.00%) were unpaired; of these:
        11277758 (99.98%) aligned 0 times
        1875 (0.02%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    0.02% overall alignment rate
    
    /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x c13624_g1_i1_esteraseE3_reversed.BWT2 -1 ../trim
    med-RNA-reads/2nd-replicate/Control-2nd-R1-trimmed-paired.fastq -2 ../trimmed-RNA-reads/2nd-replicate/Control-2nd-R2-trimmed-pa
    ired.fastq -S c13624_g1_i1_esteraseE3-RNA-control-Rep2.sam -p 10
    26277392 reads; of these:
      26277392 (100.00%) were paired; of these:
        26269843 (99.97%) aligned concordantly 0 times
        7549 (0.03%) aligned concordantly exactly 1 time
        0 (0.00%) aligned concordantly >1 times
        ----
        26269843 pairs aligned concordantly 0 times; of these:
          119 (0.00%) aligned discordantly 1 time
        ----
        26269724 pairs aligned 0 times concordantly or discordantly; of these:
          52539448 mates make up the pairs; of these:
            52539368 (100.00%) aligned 0 times
            80 (0.00%) aligned exactly 1 time
            0 (0.00%) aligned >1 times
    0.03% overall alignment rate
    /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x c13624_g1_i1_esteraseE3_reversed.BWT2 -1 ../trimmed-RNA-reads/2nd-replicate/Resistant-2nd-R1-trimmed-paired.fastq -2 ../trimmed-RNA-reads/2nd-replicate/Resistant-2nd-R2-trimmed-paired.fastq -S c13624_g1_i1_esteraseE3-RNA-resist-Rep2.sam -p 10
    24091747 reads; of these:
      24091747 (100.00%) were paired; of these:
        24087017 (99.98%) aligned concordantly 0 times
        4730 (0.02%) aligned concordantly exactly 1 time
        0 (0.00%) aligned concordantly >1 times
        ----
        24087017 pairs aligned concordantly 0 times; of these:
          85 (0.00%) aligned discordantly 1 time
        ----
        24086932 pairs aligned 0 times concordantly or discordantly; of these:
          48173864 mates make up the pairs; of these:
            48173796 (100.00%) aligned 0 times
            68 (0.00%) aligned exactly 1 time
            0 (0.00%) aligned >1 times
    0.02% overall alignment rate

**DNA alignment**


    /mnt/HD2/seq_an/bowtie2-2.2.6/bowtie2 --no-unal -x c13624_g1_i1_esteraseE3_reversed.BWT2 -1 /mnt/HD4/Genoma/1-Trimmomatic/Chom-R1-paired.fastq -2 /mnt/HD4/Genoma/1-Trimmomatic/Chom-R2-paired.fastq -S c13624_g1_i1_esteraseE3-DNA-control.sam
    154067376 reads; of these:
      154067376 (100.00%) were paired; of these:
        154066825 (100.00%) aligned concordantly 0 times
        551 (0.00%) aligned concordantly exactly 1 time
        0 (0.00%) aligned concordantly >1 times
        ----
        154066825 pairs aligned concordantly 0 times; of these:
          3 (0.00%) aligned discordantly 1 time
        ----
        154066822 pairs aligned 0 times concordantly or discordantly; of these:
          308133644 mates make up the pairs; of these:
            308133557 (100.00%) aligned 0 times
            87 (0.00%) aligned exactly 1 time
            0 (0.00%) aligned >1 times
    0.00% overall alignment rate
    

Now we need to download the sam alignments and open them with tablet.

Expressão alelo específica
Reference/consensus file:
c13624_g1_i1_esteraseE3_reversed.fasta

*C. hominivorax* (Carvalho et al., 2006)
Wild type: codon 137 (GGG, Gly); codon 251 (TGG, Trp)
Mutant: codon 137 (GAC, Asp); codon 251 (TCG, Ser)

| file                                               | codon 137 (position 691-693)                                                                                                                                                  | codon 251 (position 1033-1035)                                                                                                                                                                                    |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| (Reference) c13624_g1_i1_esteraseE3_reversed.fasta |                                                                                                                                                                               |                                                                                                                                                                                                                   |
| c13624_g1_i1_esteraseE3-DNA-control.sam            | GGG (53/98.15%), GGT (1/1.85%) Gly                                                                                                                                            | TCG (148/87.06%) Ser<br>TGG (22/12.94%) Trp                                                                                                                                                                       |
| c13624_g1_i1_esteraseE3-RNA-control-Rep1.sam       | GGG (23/43.40%) Gly<br>GAC (30/56.60%) Asp                                                                                                                                    | TCG (74/48.37%) Ser<br>TGG (78/50.98%) Trp<br>TAG (1/0.65%) Stop                                                                                                                                                  |
| c13624_g1_i1_esteraseE3-RNA-control-Rep2.sam       | GGG (403/53.24%), GGT (4/0.53%), GGC (1/0.13%) Gly<br>GAC (341/45.05%) Asp<br>GAG (2/0.26%) Glu<br>CTC (2)/0.26% Leu<br>GTG (2/0/26%), GTC (1/0.13%) Val<br>AGG (1/0.13%) Arg | TCG (503/51.86%) Ser<br>TGG (444/45.77%) Trp<br>TGT (12/1.24%) Cys<br>ACG (2/0.21%) Thr<br>CCG (2/0.21%) Pro<br>TGA (2/0.21%) Stop<br>GGG (1/0.10%), GGT (1/0.10%) Gly<br>AGG (1/0.10%) Arg<br>*GG (1)<br>T*G (1) |
| c13624_g1_i1_esteraseE3-RNA-resist-Rep1.sam        | GGG ****(34/80.95%), GGC (7/16.67%) Gly<br>GAC (1/2.38%) Asp                                                                                                                  | TCG (115/95.04%) Ser<br>TTG (5/4.13%) Leu<br>TGC (1/0.83%) Cys                                                                                                                                                    |
| c13624_g1_i1_esteraseE3-RNA-resist-Rep2.sam        | GGG (482/94.32%), GGC (19/3.72%) Gly<br>GAC (5/0.98%) Asp<br>GTG (3/0.59%) Val <br>GCG (1/0.20%) Ala<br>AGG (1/0.20%) Arg                                                     | TCG (558/89.86%), TCT (1/0.16%), TCA (1/0.16%), TCC (1/0.16%) Ser<br>TGG (34/5.48%), Trp<br>TTG (17/2.74%) Leu<br>CCG (5/0.81%) Pro<br>GCG (2/0.32%) Gly<br>CGG (1/0.16%) Arg<br>TAG (1/0.16%) Stop               |




## Get info on Esterase E3 sequences of other Calliphorids

Starting with the species used by Gi




[x] Achar a tabela dos 60 genes anotado @Gisele A 
[x] Passar o pipeline para anotação @Gisele A 
[x] Anotaçao e teste de enriquecimento de GO @Pedro M @Raquel D @Tatiana T 
[x] Encontrar espécies para o dn/ds @Raquel D @Pedro M @Tatiana T 
[x] Passar o pipeline dos ortólogos @Gisele A 
[x] Encontrar os ortólogos @Pedro M 
[x] dn/ds para todos os ortologos @Gisele A 
[x] Anotação dos SNPs com frequências diferente entre Controle e Resistente @Sophie T
[x] Busca da Esterase na anotação e nos SNPs @Sophie T


[x] GO analysis @Sophie T


[ ] Draft paper @Sophie T
[ ] 

