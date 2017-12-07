#!/usr/bin/env perl
my $command="$0 ";
my $i=0;
while (defined($ARGV[$i])){
  $command .= "$ARGV[$i] ";
  $i++;
}
use Getopt::Std;
use Cwd;
use File::Basename;
use lib dirname (__FILE__);
use ChildManager;
use strict;

#
# columns in fasta output
#
my $COLUMNS=60;

# Default parameters
*LOG=*STDERR;
my $verbose=0;
my $Verbose=0;
my $evalue_cut='1e-30';
my $percId_min=98;
my $tmpDir = "assembleContigs_$$";
my $overlap_min=100;
my $extend_min=50;
my $fastaLen_min=200;
my $wordSize_min=11;
my $cm;
my $cleanup=1;
my @rmFiles=();
my $early_stop='stop_after_this_iteration';
#
# For debugging its good to see overlapping fragments like
# Q----------XXXXX
# S          XXXXX------   
# so the $Xoffset=5 will show and additional 5 nucleotides from the subject sequence S ($Xoffset should be 0 in production mode)
my $Xoffset=0;
my $extend_new;
my $fh_in1;
my $fh_in2;
my $fh_out;
my $cores=1;
my $N_subsets=100;
my $thisDir=cwd();
#
# Process command line
#
getopts('hi:j:o:vVl:e:p:b:t:x:L:w::c:sn:S:C:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i filename] [-j filename] [-c cores] [-o filename] [-e number] [-p number] [-x number] [-w number] [-L number] [-t string]  [-S string][-l filename] [-v]\n");
  print ("Description:\n");
  print ("$0 - assemble entries from 2 fasta files\n");
  print ("### Require blastn and makeblastdb i.e. ncbi-blast/2.6.0+ ###\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input fasta file 1 (Query sequences)\n");
  print ("  -j  : input fasta file 2 (Subject sequences)\n");
  print ("  -c  : number of cores [$cores]\n");
  print ("  -o  : output fasta file [STDOUT]\n");
  print ("  -e  : max evalue [$evalue_cut]\n");
  print ("  -p  : min percent id from blastn alignments [$percId_min]\n");
  print ("  -b  : min number of bp in overlap [$overlap_min]\n");
  print ("  -x  : min number of bp to extend a sequence [$extend_min]\n");
  print ("  -w  : blastn minimum word size [$wordSize_min]\n");
  print ("  -L  : minimum length of query and subject fasta sequences [$fastaLen_min]\n");
  print ("  -t  : tmp dir [$tmpDir]\n");
  print ("  -n  : number of subsets [$N_subsets]\n");
  print ("  -S  : file apearing in temp dir resulting in early stop after current interation [$early_stop]\n");
  print ("  -C  : Cleanup [on]\n");
  print ("  -l  : logfile [STDERR]\n");
  print ("  -v  : Verbose [off]\n");
  print ("  -V  : Extended verbose [off]\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
#
if (not defined($Getopt::Std::opt_i)){
    die "option -i not defined\n";
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open($fh_in1,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open($fh_in1,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
}

if (not defined($Getopt::Std::opt_j)){
    die "option -j not defined\n";
} 
else{
  # Read from file
  if (($Getopt::Std::opt_j=~/\.gz$/) || ($Getopt::Std::opt_j=~/\.Z$/)){
    open($fh_in2,"gunzip -c $Getopt::Std::opt_j |") || die ("can't open file $Getopt::Std::opt_j: $!");
  }
  else{
    open($fh_in2,"<$Getopt::Std::opt_j") || die ("can't open file $Getopt::Std::opt_j: $!");
  }
}
#
# If not file name is given, use standard output
#
if (not defined($Getopt::Std::opt_o)){
  # Output goes to std output
  *OUT = *STDOUT;
} else {
  # Open file to write to
  open(OUT, ">$Getopt::Std::opt_o") || die ("can't open file $Getopt::Std::opt_o: $!");
}
if (defined($Getopt::Std::opt_c)){
    $cores=$Getopt::Std::opt_c;
}
if (defined($Getopt::Std::opt_n)){
    $N_subsets=$Getopt::Std::opt_n;
}
if (defined($Getopt::Std::opt_l)){
    open(LOG,">$Getopt::Std::opt_l");
}
if (defined($Getopt::Std::opt_v)){
    $verbose=1;
}
if (defined($Getopt::Std::opt_V)){
    $Verbose=1;
}
if (defined($Getopt::Std::opt_e)){
    $evalue_cut=$Getopt::Std::opt_e;
}
if (defined($Getopt::Std::opt_p)){
    $percId_min=$Getopt::Std::opt_p;
}
if (defined($Getopt::Std::opt_b)){
    $overlap_min=$Getopt::Std::opt_b;
}
if (defined($Getopt::Std::opt_x)){
    $extend_min=$Getopt::Std::opt_x;
}
if (defined($Getopt::Std::opt_t)){
    $tmpDir=$Getopt::Std::opt_t;
    if (($tmpDir eq $thisDir) || ($tmpDir eq '.')){
	print LOG "temp directory can not be the same as working directory\n";
	exit;
    }
}
if (defined($Getopt::Std::opt_L)){
    $fastaLen_min=$Getopt::Std::opt_L;
}
if (defined($Getopt::Std::opt_w)){
    $wordSize_min=$Getopt::Std::opt_w;
}
if (defined($Getopt::Std::opt_S)){
    $early_stop=$Getopt::Std::opt_S;
}
if (defined($Getopt::Std::opt_C)){
    $cleanup=0;
}
###############################################################################
# Main
#
###############################################################################
#
# Dependincies
#

my $datestring = localtime();
if ($verbose || $Verbose){
    print LOG "## Local date and time $datestring - Start program\n";
    print LOG "# $command\n";
    print LOG "# working dir: $thisDir\n\n";
    my $str="# -i : input fasta file 1:";
    printf LOG ("%-40s\t$Getopt::Std::opt_i\n", $str);
    $str="# -j : input fasta file 2";
    printf LOG ("%-40s\t$Getopt::Std::opt_j\n", $str);
    $str="# -o : output fasta file";
    printf LOG ("%-40s\t$Getopt::Std::opt_o\n", $str);
    $str="# -e : max evalue";
    printf LOG ("%-40s\t$evalue_cut\n", $str);
    $str="# -p : min percent id from blastn alignments";
    printf LOG ("%-40s\t$percId_min\n", $str);
    $str="# -b : min number of bp in overlap";
    printf LOG ("%-40s\t$overlap_min\n", $str);
    $str="# -x : min number of bp to extend a sequence";
    printf LOG ("%-40s\t$extend_min\n", $str);
    $str="# -L : min length of input sequences";
    printf LOG ("%-40s\t$fastaLen_min\n", $str);
    $str="# -t : tmp dir";
    printf LOG ("%-40s\t$tmpDir\n", $str);
    $str="# -c : number of cores";
    printf LOG ("%-40s\t$cores\n", $str);
    $str="# -n : number of subsets";
    printf LOG ("%-40s\t$N_subsets\n", $str);
    $str="# -l : logfile";
    if (defined($Getopt::Std::opt_l)){
	printf LOG ("%-40s\t$Getopt::Std::opt_l\n", $str);
    }
    else{
	printf LOG ("%-40s\tSTDERR\n", $str);
    }
    if (defined($Getopt::Std::opt_v)){
	$str="# -v : Verbose";
	printf LOG ("%-40s\ton\n", $str);
    }
    if (defined($Getopt::Std::opt_V)){
	$str="# -V : Extended verbose";
	printf LOG ("%-40s\ton\n", $str);
    }
    $str="# -S : file in tmpDir for early stopping";
    printf LOG ("%-40s\t$early_stop\n", $str);

    $str="# -C : cleanup";
    if (! defined($Getopt::Std::opt_C)){
	printf LOG ("%-40s\ton\n", $str);
    }
    else{
	printf LOG ("%-40s\toff\n", $str);
    }
    print LOG "\n\n";
}


if (! -d $tmpDir){
    print LOG "# Making tmp dir: $tmpDir\n" if ($Verbose);
    system("mkdir -p $tmpDir");
}

my $len;
my $k=0;
#
# Fasta file with query sequences
#

my $Qfasta;
my $Sfasta;
my $entries_S;
my $entries_Q;
my $entries_ok_S;
my $entries_ok_Q;
my $entries_ok;
my $iterDir;
my $cmd;

my @w=();
my $Qname;
my $Sname;
my $Qlen;
my $Qstart;
my $Qend;
my $Slen;
my $Sstart;
my $Send;
my $overlap;
my $extend;
my $QS;
my $type;

my %rec=();
my %printout=();
my %idxQ=();
my %idxS=();

my $Qseq='';
my $Sseq='';
my $Qheader='';
my $Sheader='';
my $string='';
my $wordSize;
my $rh_idxQ;
my $rh_idxS;


$cm = new ChildManager($cores);




my $wordSize=&calcWordSize();


my $contigFsa="$iterDir/contigs.fsa";
my $deletedSequences=0;
my $extended_sequences=0;
my $iteration=1;
my $first_iteration_flag=1;
my $makePrelimFsa_flag=1;


my $statusFile="$tmpDir/status.tab";
my $fh_status;
if (-e $statusFile){
    ($makePrelimFsa_flag, $iteration) = &checkStatus($statusFile);
}
else{
    $iterDir="$tmpDir" . "/" . "$iteration";
    if (! -d "$iterDir"){
	system("mkdir -p $iterDir");
    }
    $Qfasta="$iterDir/Q.fsa";
    $Sfasta="$iterDir/S.fsa";
}


my $fh_status;
if ($makePrelimFsa_flag){
    if ($iteration > 1){
	open($fh_status, '>>' , $statusFile) || die "Can not open $statusFile for writing: $!\n";    
    }
    else{
	open($fh_status, '>' , $statusFile) || die "Can not open $statusFile for writing: $!\n";
    }

    #
    # make 1 combined contig fasta file by iterative blastn runs to remove/extend contig sequences
    #
    &makePrelimFsa();
    close($fh_status);
}


if (defined($Getopt::Std::opt_o)){
    system("cp $contigFsa $Getopt::Std::opt_o");
}
else{
    system("cat $contigFsa");
}

if (! -e "$tmpDir/$early_stop"){
    if ($cleanup){
	system("rm -r $tmpDir");
    }
}

if ($verbose || $Verbose){
    $datestring = localtime();
    print LOG "## Local date and time $datestring - End program\n";
}
exit(0);

##########################################
#
# End main prog
#
##########################################
sub makePrelimFsa{

    while (($extended_sequences || $first_iteration_flag) && (! -e "$tmpDir/$early_stop")){
    
	$extended_sequences=0;
	$first_iteration_flag=0;
	
	#
	# Make a bp lenght selection fasta entries - only iteration 1
	#
	if ($iteration == 1){
	    print LOG "# Selecting fasta entries longer than $fastaLen_min bp\n" if ($Verbose);
	    #
	    # Fasta file with Query sequences
	    #
	    if (! -e $Qfasta){
		print LOG "# Making $Qfasta\n" if ($Verbose);
		open($fh_out,">$Qfasta");
		($entries_Q, $entries_ok_Q)=&parseFSA($fh_in1, $fh_out);
		print LOG "# Sequences in $Getopt::Std::opt_i:\t$entries_Q\n" if ($verbose || $Verbose);
		print LOG "# Subset where len >= $fastaLen_min: $entries_ok_Q\n" if ($verbose || $Verbose);
		close($fh_in1);
		close($fh_out);
		if ($verbose){
		    my $N;
		    my $Min;
		    my $Max;
		    my $N50;
		    my $deleted=$entries_Q-$entries_ok_Q;
		    ($N, $Min, $Max, $N50)=&n50($Qfasta);
		    print LOG "input=$Getopt::Std::opt_i\tSequences= $N\tShortest= $Min\tLongest= $Max\tDeleted= $deleted\tExtended= - $N50\n" if ($verbose || $Verbose);
		    
		}
	    }
	    else{
		print LOG "# File already made: $Qfasta\n" if ($Verbose);
	    } 
	    
	    
	    #
	    # Fasta file with subject sequences
	    #
	    if (! -e $Sfasta){
		print LOG "# Making $Sfasta\n" if ($Verbose);
		open($fh_out,">$Sfasta");
		($entries_S, $entries_ok_S)=&parseFSA($fh_in2, $fh_out);
		print LOG "# Sequences in $Getopt::Std::opt_j:\t$entries_S\n" if ($verbose || $Verbose);
		print LOG "# Subset where len >= $fastaLen_min: $entries_ok_S\n" if ($verbose || $Verbose);
		close($fh_in2);
		close($fh_out);
		
		
		if ($Getopt::Std::opt_i eq $Getopt::Std::opt_j){
		    print LOG "# Query and Subject fasta file names are identical\n" if ($Verbose);
		}
		else{
		    if ($verbose){
			my $N;
			my $Min;
			my $Max;
			my $N50;
			my $deleted=$entries_S-$entries_ok_S;
			($N, $Min, $Max, $N50)=&n50($Sfasta);
			print LOG "input=$Getopt::Std::opt_j\tSequences= $N\tShortest= $Min\tLongest= $Max\tDeleted= $deleted\tExtended= - $N50\n" if ($verbose || $Verbose);
			
		    }

		    print LOG "# Query and Subject fasta file names are different\n" if ($Verbose);
		    print LOG "# Making one join subject fasta file\n";
		    system("cat $Qfasta >> $Sfasta");
		    if ($verbose || $Verbose){
			my $N;
			my $Min;
			my $Max;
			my $N50;
			my $deleted=$entries_S-$entries_ok_S + $entries_Q-$entries_ok_Q;
			print LOG "# Input fasta file 1+2 with length > $fastaLen_min bp\n" if ($verbose || $Verbose);
			($N, $Min, $Max, $N50)=&n50($Sfasta);
			print LOG "Combined fasta\tSequences= $N\tShortest= $Min\tLongest= $Max\tDeleted= $deleted\tExtended= - $N50\n" if ($verbose || $Verbose);
		    }
		}
	    }
	    else{
		print LOG "# File already made: $Sfasta\n" if ($Verbose);
	    }
	}
	if ($cleanup){
	    push(@rmFiles, $Sfasta);
	    push(@rmFiles, $Qfasta);
	}
	#
	# Split Query sequences into $N_subsets chunks
	#
	$datestring = localtime();
	print LOG "## Local date and time $datestring - Starting split $Sfasta into $N_subsets sets\n" if ($verbose);
	&splitFSA($Sfasta, $N_subsets);
	$datestring = localtime();
	print LOG "## Local date and time $datestring - Finished split $Sfasta into $N_subsets sets\n" if ($verbose);
	
	my $in= $Sfasta;
	my $db= "$Sfasta" . '.nin';
	if (! -e $db){
	    $datestring = localtime();
	    print LOG "## Local date and time $datestring - Starting formating blast databases\n" if ($verbose);
	    $cmd= "makeblastdb -dbtype 'nucl' -in $in -out $in";
	    print LOG "# $cmd\n" if ($verbose);
	    system("$cmd");
	    $datestring = localtime();
	    print LOG "## Local date and time $datestring - Finished formating blast databases\n" if ($verbose);
	    
	    #
	    # cleanup blastdb files
	    #
	    my $blast_db_file="$in" . '.nsq';
	    if (! -e $blast_db_file){
		$blast_db_file="$in" . '.*.nsq';
	    }
	    push(@rmFiles, $blast_db_file);
	    
	    $blast_db_file="$in" . '.nhr';
	    if (! -e $blast_db_file){
		$blast_db_file="$in" . '.*.nhr';
	    }
	    push(@rmFiles, $blast_db_file);
	    
	    $blast_db_file="$in" . '.nin';
	    if (! -e $blast_db_file){
		$blast_db_file="$in" . '.*.nin';
	    }
	    push(@rmFiles, $blast_db_file);

	    $blast_db_file="$in" . '.nal';
	    if (-e $blast_db_file){
		push(@rmFiles, $blast_db_file);
	    }
	}
			
	my @commands=&prepareBlastCommands($Sfasta, $Sfasta, $iterDir);
	#
	# run blast
	#
	$datestring = localtime();
	print LOG "## Local date and time $datestring - Starting blastn iteration=$iteration\n" if ($verbose);
	
	foreach my $command (@commands){
	    print LOG "# $command\n" if ($Verbose);
	    $cm->start("$command");
	}
	#
	# Wait for all jobs to finish
	#
	$cm->wait_all_children;
	$datestring = localtime();
	print LOG "## Local date and time $datestring - Finished blastn iteration=$iteration\n" if ($verbose);
	
	
	open(BLAST,"cat $iterDir/blast.*.out |");
	# my $fh_deleted;
	# open($fh_deleted,">$iterDir/deleted_entries.txt");
	
	
	while (<BLAST>){
	    chomp;
	    @w=split(/\s+/);
	    $Qname=$w[0];
	    $Sname=$w[1];
	    $Qlen=$w[2];
	    $Qstart=$w[3];
	    $Qend=$w[4];
	    $Slen=$w[5];
	    $Sstart=$w[6];
	    $Send=$w[7];
	    
	    if ($Qname eq $Sname){
		next;
	    }
	    
	    if (exists($printout{$Qname})){
		if ($printout{$Qname} == 0){
		    next;
		}
	    }
	    my $complete_subset=0;
	    $complete_subset=&Q_subset_of_S();
	    
	    if ($complete_subset){
		if (! exists($printout{$Qname})){		
		    $printout{$Qname}=0;
#		print $fh_deleted "delQ\t$_\n";
		    next;
		}
	    }
	    
	    if (exists($printout{$Sname})){
		if ($printout{$Sname} == 0){
		    next;
		}
	    }
	    
	    $complete_subset=&S_subset_of_Q();
	    if ($complete_subset){
		if (! exists($printout{$Sname})){		
		    $printout{$Sname}=0;
#		print $fh_deleted "delS\t$_\n";
		    next;
		}
	    }
	    
	    &get_type();
	    #
	    # subject on plus strand i.e. sstart<send
	    #
	    my $id='5-end';
	    if (($type == 1) || ($type==2)){
		$id='3-end';
	    }
	    
	    # extend the longest of the sequences Query or Subject and only keep longest of (type1 or type2) and (type3 or type4)
	    if ($type){
		$QS="$Qname" . '_' . "$id";
		if (exists ($rec{extend}{$QS})){
		    if ($extend > $rec{extend}{$QS}){
			$rec{extend}{$QS}=$extend;
			$rec{line}{$QS}="$type\t$_" . "\t$overlap";
		    }
		}
		else{
		    $rec{extend}{$QS}=$extend;
		    $rec{line}{$QS}="$type\t$_" . "\t$overlap";
		}
	    }
	}
	close(BLAST);
	open(SET,"| sort -k2 -k1nr  > $iterDir/seq_to_be_extended.txt");
	#print SET "# Type\tQname\tSname\tQlen\tQstart\tQend\tSlen\tSstart\tSend\tOverlap\tExtend\n";
	#print SET "# Type 1: Extend query sequence in 3'-end with plus strand of unaligned subject sequence\n";
	#print SET "# Type 2: Extend query sequence in 3'-end with revComp strand of unaligned subject sequence\n";
	#print SET "# Type 3: Extend query sequence in 5'-end with plus strand of unaligned subject sequence\n";
	#print SET "# Type 4: Extend query sequence in 5'-end with revComp strand of unaligned subject sequence\n";
	
	my %tst=();
	foreach my $key (keys (%{$rec{line}})){
	    my @w=split(/\s+/);
	    my $Qname=$w[1];
	    my $Sname=$w[2];
	    if (exists($printout{$Qname})){
		if ($printout{$Qname} == 0){
		    next;
		}
	    }
	    if (exists($printout{$Sname})){
		if ($printout{$Sname} == 0){
		    next;
		}
	    }
	    
	    #
	    # not allowed to extend a with b and later on extend b with a as that will extend until INFINITY
	    #
	    if (! exists($tst{$Qname})){
		$tst{$Qname}=$Sname;
		if (exists($tst{$Sname})){
		    if ($tst{$Sname} eq "$Qname"){
			next;
		    }
		}
	    }
	    print SET "$rec{line}{$key}\t$rec{extend}{$key}\n";
	}
	close(SET);
	%tst=();
	
	#
	# make index for the $iterDir/S.fsa file
	#
	print LOG "# Index fasta $Sfasta and store in hash\n" if ($Verbose);
	$rh_idxS=&makeFastaIndex($Sfasta);
	
	#
	# make index for the $iterDir/Q.fsa file
	#
	print LOG "# Index fasta $Qfasta and store in hash\n" if ($Verbose);
	$rh_idxQ=&makeFastaIndex($Qfasta);
	
	%rec=();
	open(SET,"<$iterDir/seq_to_be_extended.txt");
	my $i;
	print LOG "# Storing alignment values in hash\n" if ($Verbose);
	
	while (<SET>){
	    chomp;
	    @w=split(/\s+/);
	    $Qname=$w[1];
	    if (exists($rec{$Qname}[0])){
		$rec{$Qname}[0]++;
		my $i = $rec{$Qname}[0];
		$rec{$Qname}[$i]=$_;
	    }
	    else{
		$i=1;
		$rec{$Qname}[0]=1;
		$rec{$Qname}[$i]=$_;	
	    }
	}
	close(SET);
	
	open(Q,"<$Qfasta");
	open(S,"<$Sfasta");
	
	my $extended_fsa="$iterDir/extended.fsa";
	my $not_extended_fsa="$iterDir/not_extended.fsa";
	open(my $fh_extended,">$extended_fsa");
	open(my $fh_not_extended,">$not_extended_fsa");
	push(@rmFiles, $extended_fsa);
	push(@rmFiles, $not_extended_fsa);
	
	print LOG "# Assemble and print output contig fasta file\n\n" if ($Verbose);
	
	foreach my $key (keys %{$rh_idxQ}){
	    if (exists($printout{$key})){
		if ($printout{$key} == 0){
		    next;
		}
	    }
	    seek(Q,$rh_idxQ->{$key},0);
	    $_=<Q>;
	    chomp;
	    if (m/^>(\S+)/){
		$Qheader=$1;
	    }
	    else{
		die "error in header '$_' in file: $Qfasta trying to extract '$key' at position $rh_idxQ->{$key}\n";
	    }
	    
	    #
	    # Are there alignments of Q such that it can be extended, then $rec{$key}[0] is 1 or 2 i.e. longestextension of 5' and 3' ends
	    # thus get the query sequence such that is can be extended with Subject hit based on alignment
	    if (exists($rec{$key}[0])){
		$Qseq='';
		# store sequence in $seq;
		while (<Q>){
		    chomp;
		    if (m/^>/){
			last;
		    }	    
		    $Qseq .= $_
		}
		#
		# now get subject sequence to extend 5' and/or 3' end of query sequence
		# the value 'type' define which end to extend and if its plus or minus strand of subject sequence
		#
		$extend_new=0;
		for (my $i=1;$i<=$rec{$key}[0];$i++){
		    @w=split(/\s+/,$rec{$key}[$i]);
		    $type= $w[0];
		    $Qname= $w[1];
		    $Sname= $w[2];
		    $Slen= $w[6];
		    $Sstart=$w[7];
		    $Send=$w[8];
		    $extend=$w[-1];
		    if ($type == 1){
			#
			# Type 1: Extend query sequence in 3'-end with plus strand of unaligned subject sequence
			#
			&extend_type1($Sname);
		    }
		    elsif ($type == 2){
			#
			# Type 2: Extend query sequence in 3'-end with revComp strand of unaligned subject sequence
			#
			&extend_type2($Sname);
		    }
		    elsif ($type==3){
			#
			# Type 3: Extend query sequence in 5'-end with plus strand of unaligned subject sequence
			#
			&extend_type3($Sname);
		    }
		    elsif ($type==4){
			#
			# Type 4: Extend query sequence in 5'-end with revComp strand of unaligned subject sequence
			#
			&extend_type4($Sname);
		    }
		    else{
			die "unrecognized type has been read (type can be 1-4): '$rec{$key}[$i]'\n";
		    }
		    
		    #
		    # Now mark subject sequence to be deleted - NB! a query can not be extended in both ends by the same subject sequence, as then the query seq
		    # would have been deleted as is would be a complete subset of the subject sequence
		    #
		    $printout{$Sname}=0;
		}
		#
		# print out the extended Query sequence, and mark the old query sequence to be deleted
		#
		if ($type){
		    $printout{$Qname}=0;
		    $Qheader .= "_" . "E$extend_new";
		    $printout{$Qheader}=1;
		    &string2fasta($fh_extended, $Qheader, $Qseq);
		    $extended_sequences++;
		}
	    }
	    else{
		#
		# The query sequence that has not been extended so just write it out
		#
		if (! exists($printout{$Qheader})){
		    $printout{$Qheader}=1;
		    print $fh_not_extended ">$Qheader\n";
		    while (<Q>){
			if (m/^>/){
			    last;
			}
			print $fh_not_extended "$_";
		    }
		}
	    }
	}
	#
	# printout all database sequences that were not already used to extend a query sequence
	# only do this in iteration 1
	
	if ($iteration == 1){
	    foreach my $key (keys %$rh_idxS){
		seek(S,$rh_idxS->{$key},0);
		if (exists($printout{$key})){
		    if ($printout{$key}==0){
			next;
		    }
		}
		else{
		    $_=<S>;
		    chomp;
		    print $fh_not_extended "$_\n";
		    while (<S>){
			if (m/^>/){
			    last;
			}	    
			print $fh_not_extended "$_"
		    }
		}
	    }
	}
	close(Q);
	close(S);
	close(OUT);
	close($fh_extended);
	close($fh_not_extended);
#    close($fh_deleted);
	$contigFsa="$iterDir/contigs.fsa";
	
	system("cat $extended_fsa $not_extended_fsa > $contigFsa");
	
	my $N;
	my $Min;
	my $Max;
	my $N50;
	if ($verbose || $Verbose){
	    foreach my $key (keys %printout){
		if ($printout{$key} == 0){
		    $deletedSequences++;
		}
	    }
	    my ($N, $Min, $Max, $N50)=&n50($contigFsa);
	    print LOG "$iteration\tSequences= $N\tShortest= $Min\tLongest= $Max\tDeleted= $deletedSequences\tExtended= $extended_sequences $N50\n";
	}
	print $fh_status "Extended\t$extended_sequences\t$iteration\n";
	
	if ($extended_sequences == 0){
	    last;
	}
	else{
	    if ($iteration > 1){
		push(@rmFiles, $Sfasta);
		push(@rmFiles, $Qfasta);
	    }
	    
	    $iteration++;
	    %idxQ=();
	    %idxS=();
	    %rec=();
	    %printout=();
	    $iterDir="$tmpDir" . '/' . "$iteration";
	    if (! -d "$iterDir"){
		system("mkdir -p $iterDir");
	    }
	    $Qfasta="$iterDir" . '/' . 'Q.fsa';
	    $Sfasta="$iterDir" . '/' . 'S.fsa';	
	    
	    system("cp $contigFsa $Sfasta");
	    system("cp $contigFsa $Qfasta");

	    $deletedSequences=0;
	    if (($cleanup) && (! -e "$tmpDir/$early_stop")){
		foreach my $file (@rmFiles){
		    if (-e $file){
			system("rm $file")
		    }
		    else{
			print LOG "Can not remove file: $file\n" if ($Verbose);
		    }
		}
	    }
	    @rmFiles=();
	}
    }
    return;
}

sub get_type{

    #
    # Query is always on plus strand
    #
    $type=0;
    $overlap=0;
    $extend=0;

    if ($Qname eq $Sname){
	return($type);
    }
    if ($Qend == $Qlen){
	#
	# Add to 3'-end of query sequence
	#
	if ($Sstart < $Send){
	    # Subject on plus strand
	    
	    if (($Send != $Slen) && ($Sstart == 1)){
		$overlap = $Send - $Sstart +1;
		$extend = $Slen - $Send +1;
		if (($overlap >= $overlap_min) && ($extend >= $extend_min)){
		    # Extend query sequence in 3'-end with plus strand of unaligned subject sequence
		    $type =1;
		}
	    }
	}
	else{
	    # Subject on minus strand
	    if ( ($Send != 1) && ($Sstart==$Slen)){
		$overlap = $Sstart - $Send +1;
		$extend = $Slen - $overlap;
		if (($overlap >= $overlap_min) && ($extend >= $extend_min)){
		    # Extend query sequence in 3'-end with revComp strand of unaligned subject sequence
		    $type=2;
		}
	    }
	}
    }
    elsif ($Qstart == 1){
	#
	# Add to 5'-end of query sequence
	#
	if ($Sstart < $Send){
	    # Subject on plus strand
	    
	    if (($Send == $Slen) && ($Sstart != 1)){
		$overlap = $Send - $Sstart + 1;
		$extend = $Slen - $overlap;;
		if (($overlap >= $overlap_min) && ($extend >= $extend_min)){
		    # Extend query sequence in 5'-end with plus strand of unaligned subject sequence
		    $type =3;
		}
	    }
	}
	else{
	    # Subject on minus strand
	    if (($Send == 1) && ($Sstart != $Slen)){
		$overlap = $Sstart - $Send +1;
		$extend = $Slen - $overlap;
		if (($overlap >= $overlap_min) && ($extend >= $extend_min)){
		    # Extend query sequence in 5'-end with revComp strand of unaligned subject sequence
		    $type=4;
		}
	    }
	}
	
    }
    return($type);
}

sub parseFSA{
    my ($fh_in, $fh_out)=@_;

    my $line = '';
    my $header;
    my $len=0;
    my $seq;
    my $entries=0;
    my $entries_ok=0;
    # Finding the header line ">....."
    $line = <$fh_in> while (defined $line and $line !~ m/^>/);


    while (defined $line) {
	# line starts with '>'
	chomp $line;
	$header=$line;
	$seq = '';
	while (defined ($line = <$fh_in>) and $line !~ m/^>/) {
	    chomp $line;
	    $seq .= $line;
	}
	$entries++;
	$len=length($seq);
	if ($len >= $fastaLen_min){
	    $entries_ok++;
	    print $fh_out "$header\n";
	    for (my $i = 0; $i < $len; $i += 60) {
		print $fh_out substr($seq, $i, 60), "\n";
	    }
	}
    }
    return($entries, $entries_ok);
}


sub string2fasta{
    my ($fh, $header, $seq) = @_;

    my $len = length($seq);
    my $lines= int($len/$COLUMNS);

    print $fh ">$header\n";

    for (my $i=0;$i<=$lines;$i++){
	print $fh substr($seq,$i*$COLUMNS,$COLUMNS);
	if ($i != $lines){
	    print $fh  "\n";
	}
	else{
	    if ($len%$COLUMNS !=0){
		print $fh "\n";
	    }
	}
    }		
    return(0);
}

sub extend_type1{
    my ($name) =@_;
    #
    # Extend Query with 3'-end of Subject sequence on plus strand
    # Q -----...........b
    # S XXXXXXXa........c........d
    # position: b= Qlen in Query sequence
    # position: a= 1    in Subj sequence
    # position. d= Slen in Subj sequence
    # positions c+1 -> d can be extended to 3'-end of query
    seek(S,$rh_idxS->{$name},0);
    $_=<S>;
    chomp;
    if (! m/^>(\S+)/){
	die "error in header '$_' in file: $Sfasta trying to extract '$name' at position $rh_idxS->{$name}\n";
    }
    #
    # length by which the query sequence is extended
    #
    $extend_new += $extend;
    $Sseq='';
    while (<S>){	    
	chomp;
	if (m/^>/){
	    last;
	}
	$Sseq .= $_
    }
#    $string = lc(substr($Sseq,$Send-$Xoffset,$Slen-1));
    $string = substr($Sseq,$Send-$Xoffset,$Slen-1);
    $Qseq .= $string;
    return(0);
}

sub extend_type2{    
    my ($name)=@_;

    seek(S,$rh_idxS->{$name},0);
    $_=<S>;
    chomp;
    if (! m/^>(\S+)/){
	die "error in fasta\n";
    }
    
    #
    # length by which the query sequence is extended
    #
    $extend_new += $extend;
    $Sseq='';
    while (<S>){	    
	chomp;
	if (m/^>/){
	    last;
	}
	$Sseq .= $_
    }
    my $rev='';
#    $string = lc(substr($Sseq,0,$Send-1+$Xoffset));
    $string = substr($Sseq,0,$Send-1+$Xoffset);
    $rev = reverse($string);
    $rev =~ tr/ACGTatgc/TGCAtacg/;
    $Qseq .= $rev;
    return(0);
}

sub extend_type3{
    my ($name)=@_;
    seek(S,$rh_idxS->{$name},0);
    $_=<S>;
    chomp;
    if (! m/^>(\S+)/){
	die "Error in fasta\n";
    }
    
    #
    # length by which the query sequence is extended
    #
    $extend_new += $extend;
    $Sseq='';
    while (<S>){	    
	chomp;
	if (m/^>/){
	    last;
	}
	$Sseq .= $_
    }
#    $string = lc(substr($Sseq,0,$Sstart-1+$Xoffset));
    $string = substr($Sseq,0,$Sstart-1+$Xoffset);
    $Qseq = $string . $Qseq;
    return(0);    
}

sub extend_type4{
    my ($name) = @_;
    seek(S,$rh_idxS->{$name},0);
    $_=<S>;
    chomp;
    if (! m/^>(\S+)/){
	die "error in fasta\n";
    }
    
    #
    # length by which the query sequence is extended
    #
    $extend_new += $extend;
    $Sseq='';
    while (<S>){	    
	chomp;
	if (m/^>/){
	    last;
	}
	$Sseq .= $_
    }
    my $rev='';
#    $string = lc(substr($Sseq,$Sstart-$Xoffset));
    $string = substr($Sseq,$Sstart-$Xoffset);
    $rev = reverse($string);
    $rev =~ tr/ACGTatgc/TGCAtacg/;
    $Qseq = $rev . $Qseq;
    return(0);
}

sub n50{
    my ($fasta)=@_;
    
    if (($fasta=~/\.gz$/) || ($fasta=~/\.Z$/)){
	open(FSA,"gunzip -c $fasta |");
    }
    else{
	open(FSA,"<$fasta");
    }
    
    my $sum=0;
    my %bp=();
    my $entries=0;
    my $max = -99999;
    my $min = 999999;
    my $seq='';;
    my $len=0;
    my $first=1;
    while (defined($_=<FSA>)){
	chomp;
	if (m/^>/){
	    $entries++;
	    if ($first){
		$first=0;
		next;
	    }
	    $len=length($seq);
	    if ($len > $max){
		$max=$len;
	    }
	    if ($len < $min){
		$min=$len;
	    }
	    $sum += $len;
	    $bp{$len}++;
	    $seq='';
	    next;
	}
	$seq .= $_;
    }
    # add for the last entry
    $len=length($seq);
    $sum += $len;
    $bp{$len}++;
    if ($len > $max){
	$max=$len;
    }
    if ($len < $min){
	$min=$len;
    }
    $seq='';
    
    my $val=0;
    my @nx_values=();
    my %nx=();
    my $n10=int($sum*1/10 +0.5);
    my $n20=int($sum*2/10 +0.5);
    my $n30=int($sum*3/10 +0.5);
    my $n40=int($sum*4/10 +0.5);
    my $n50=int($sum*5/10 +0.5);

    foreach my $key (sort {$b <=> $a} keys %bp){
	$val += $bp{$key}*$key;

	if (($val >= $n10) && (! exists($nx{10}))){
	    $nx{10}=$key;
	    push(@nx_values, $key);
	}
	if (($val >= $n20) && (! exists($nx{20}))){
	    $nx{20}=$key;
	    push(@nx_values, $key);
	}
	if (($val >= $n30) && (! exists($nx{30}))){
	    $nx{30}=$key;
	    push(@nx_values, $key);
	}
	if (($val >= $n40) && (! exists($nx{40}))){
	    $nx{40}=$key;
	    push(@nx_values, $key);
	}
	if (($val >= $n50) && (! exists($nx{50}))){
	    $nx{50}=$key;
	    push(@nx_values, $key);
	    last;
	}
    }
    close(FSA);
    my $str='';
    my $i=10;
    foreach my $id (@nx_values){
	$str .= "\tn$i= $id";
	$i +=10;
    }
    
    return($entries, $min, $max, $str);
}

sub makeFastaIndex{
    my ($fasta) = @_;

    my %idx=();
    open(FSA,"<$fasta") || (die "file not found: $fasta\n");

    my $pos=0;
    while (<FSA>){
	if (m/^>(\S+)/){
	    $pos=tell(FSA);
	    $idx{$1}=$pos-length($_);
	}
    }
    close(FSA);
    return(\%idx);
}

sub Q_subset_of_S{
    my $value=0;
    if ($Qname ne $Sname){
	if (($Qstart == 1) && ($Qend == $Qlen)){
	    $value=1;
	}
    }
    return($value);
}

sub S_subset_of_Q{
    my $value=0;
    if ($Qname ne $Sname){
	if ( ( ($Sstart == 1) && ($Send == $Slen)) || (($Sstart==$Slen) && ($Send==1))){
	    $value=1;
	}
    }
    return($value);
}

sub splitFSA{
    my ($fsa, $N) = @_;

    my $i=0;
    my $fh;
    my $fh_out;
    my %filehandle=();
    my $rest=0;
    my $max_lines=5000000;
    my $lines=0;
    my %buffer=();

    open($fh, "<", "$fsa") || die "cannot open $fsa: $!";
    
    for my $l (0..$N-1) {	
	my $out= "$fsa" . "." . "$l"; 	
	open(my $fh_out, ">", $out) || die "cannot open $out: $!";
	$filehandle{$l}=$fh_out;	
	if ($cleanup){
	    push(@rmFiles, $out);
	}
    }
    $fh_out = $filehandle{0};
    while (defined($_=<$fh>)){
	if (m/^>/){
	    $rest= $i%$N;
	    $fh_out=$filehandle{$rest};
	    $i++;
	}
	$buffer{$rest} .= $_;
	$lines++;
	if ($lines > $max_lines){
	    #
	    # empty buffer
	    #
	    for my $l (0..$N-1) {
		$fh_out=$filehandle{$l};
		print $fh_out $buffer{$l};		
	    }
	    %buffer=();
	    $lines=0;
	}
    }
    close($fh);

    for my $l (0..$N-1) {
	$fh_out=$filehandle{$l};
	print $fh_out $buffer{$l};
	close($fh_out);
    }
    %buffer=();
    return(0);
}

sub calcWordSize{
    my $wordSize=$wordSize_min;
    my $HSP_min=1;
    my $HSP=$HSP_min;
    my $ws=$wordSize_min;
    while ($HSP >= $HSP_min){
	$HSP = $overlap_min - $ws + 1 - (1-$percId_min/100) * $ws * $overlap_min;
	if ($HSP >= $HSP_min){
	    printf LOG "# word size= $ws\tHSP= %.0f\n",$HSP if ($Verbose);
	    $wordSize=$ws;
	}
	$ws++;
    }
    print LOG "# wordsize needed to find region with id>=$percId_min% with an overlap>= $overlap_min:\t$wordSize\n" if ($verbose || $Verbose);
    return($wordSize);
}

sub checkStatus{
    my ($file)=@_;

    my $fh;
    open($fh, '<' , $file) || die "Can not open $file for reading: $!\n";
    
    #
    # format of $file:
    # 3 tab separated columns
    # column 1: 'Extended'
    # column 2: number of extended sequences
    # column 3: iteration cycle
    #
    while (<$fh>){	
	chomp;
	@w=split(/\t/);
	if ($w[1] == 0){
	    $makePrelimFsa_flag=0;
	}
	$iteration=$w[2];
    }
    close($fh);
    return($makePrelimFsa_flag, $iteration);
}


sub prepareBlastCommands{
    my ($query, $db, $outDir, $alignment_flag) = @_;
    my @commands=();
    #
    # prepare to run blast
    #
    my $cmd;
    for (my $i=0;$i<$N_subsets;$i++){
	my $in= "$query". '.' . "$i";
	my $out= "$outDir/blast.$i.out";
	if ($cleanup){
	    push(@rmFiles, $out);
	}
	if (! -e $out){
	    if (! $alignment_flag){
		$cmd = "blastn -task blastn -word_size $wordSize -query $in -db $db -evalue $evalue_cut -perc_identity $percId_min -out $out -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send evalue pident'";
	    }
	    else{
		$cmd = "blastn -task blastn -word_size $wordSize -query $in -db $db -evalue $evalue_cut -perc_identity $percId_min -out $out -outfmt '6 qseqid sseqid qlen qstart qend slen sstart send evalue pident qseq sseq'";
	    }
	    push(@commands, $cmd);
	}
    }
    return(@commands);
}

sub getEntry{
    my ($filehandle, $rh, $name) = @_;

    seek($filehandle,$rh->{$name},0);
    $_=<$filehandle>;
    chomp;

    if (! m/^>(\S+)/){
	die "Error in fasta header - subroutine getEntry\n";
    }
	
    #
    # length by which the query sequence is extended
    #

    my $sequence='';
    while (<$filehandle>){	    
	chomp;
	if (m/^>/){
	    last;
	}
	$sequence .= $_
    }
    return($sequence);
}

