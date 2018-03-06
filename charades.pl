#!/usr/bin/perl
use warnings;
use strict;
use Cwd;
use Getopt::Long;

# This script looks at the first 100 000 sequences from a bisulfite split fastq file and then
# tries to find out which how the library was created. It will use sequence composition at different
# positions in the sequence.

my $debug;
my @files_to_analyse;
my $parent_directory = getcwd();
	
process_commandline();
files_to_analyse();
extract_100K_reads();
collect_base_composition_info();
guess_library();



###################
### SUBROUTINES ###
###################

sub process_commandline {

	my $help;
	my $output_dir;
	
	my $command_line =  GetOptions ('help'                => \$help,
									'output_dir=s'        => \$output_dir,
									);
	
	
	die "Please respecify command line options\n\n" unless ($command_line);
	die "\nPlease specifiy gzipped fastq files\n\n" unless (@ARGV);

    if ($help){
	print_helpfile();
	exit;
    }
}

sub print_helpfile{
	print "\n\nCharades will try and predict which bisulfite sequencing library preparation was used. It uses gzipped fastq files as input.\n\n";
	print ">>> USAGE: ./charades.pl [options] filename(s) <<<\n\n";
	print "--help\t\t\tprint this helpfile\n";
	print "--output_dir [path]\tOutput directory, either relative or absolute. Output is written to the current directory if not\n\t\t\tspecified explicitly.\n\n";
}


sub files_to_analyse {											

## Finding out which files to analyse

	my @files = @ARGV;
	my $length = @files;
	warn "\nScanning $length input file(s).\n";
	warn join("\n",@files)."\n" if $debug;


	foreach my $file(@files) {
	
		
	
		if ($file =~ /R1/) {
			push @files_to_analyse, $file;
		}
		elsif ($file !~ /R2|R3|R4/) {
			push @files_to_analyse, $file;
		}
	}
	warn "\nFiles selected for analysis are\n", join("\n",@files_to_analyse),"\n\n\n";
}




sub extract_100K_reads {

	my @subfolders;										
	system "mkdir charades_temp";
	
	warn "\nExtracting the first 100 000 reads from fastq file\n\n";
	foreach my $file_to_analyse(@files_to_analyse) { 
	
		my $file_name;											## Extracting the filenames if a path is given
		while ($file_to_analyse =~ /\//g) {						
			push @subfolders, pos$file_to_analyse;
			warn "The last slash is found at ", $subfolders[-1], "\n" if $debug;
		}
		
	
		if (scalar @subfolders > 0) {							## remove the path from the file name							
			$file_name = substr($file_to_analyse, ($subfolders[-1]));
		}
		else {$file_name = $file_to_analyse};					## if files are in the current directory
		system "zcat $file_to_analyse | head -400000 > charades_temp/$file_name.100K.fastq";
	}
}




sub collect_base_composition_info {
	
	open (my $out, '>', "sequence_composition_stats.csv") or die "Cannot write to file: $!";
	open (my $out_arff, '>', "sequence_composition_stats.arff") or die "Cannot write to file: $!";
	print $out "file_name,pos1_A,pos1_C,pos1_T,pos1_G,pos2_A,pos2_C,pos2_T,pos2_G,pos3_A,pos3_C,pos3_T,pos3_G,pos4_A,pos4_C,pos4_T,pos4_G,pos8_A,pos8_C,pos8_T,",
				"pos8_G,pos9_A,pos9_C,pos9_T,pos9_G,pos10_A,pos10_C,pos10_T,pos10_G,pos11_A,pos11_C,pos11_T,pos11_G,pos12_A,pos12_C,pos12_T,pos12_G,",
				"pos20_A,pos20_C,pos20_T,pos20_G,pos30_A,pos30_C,pos30_T,pos30_G,method\n";
	#print $out "file_name\tposition\tpercentA\tpercentC\tpercentT\tpercentG\tpercentN\n";
	
	print $out_arff "\@relation training_data_summary_sorted\n
\@attribute pos1_A numeric
\@attribute pos1_C numeric
\@attribute pos1_T numeric
\@attribute pos1_G numeric
\@attribute pos2_A numeric
\@attribute pos2_C numeric
\@attribute pos2_T numeric
\@attribute pos2_G numeric
\@attribute pos3_A numeric
\@attribute pos3_C numeric
\@attribute pos3_T numeric
\@attribute pos3_G numeric
\@attribute pos4_A numeric
\@attribute pos4_C numeric
\@attribute pos4_T numeric
\@attribute pos4_G numeric
\@attribute pos8_A numeric
\@attribute pos8_C numeric
\@attribute pos8_T numeric
\@attribute pos8_G numeric
\@attribute pos9_A numeric
\@attribute pos9_C numeric
\@attribute pos9_T numeric
\@attribute pos9_G numeric
\@attribute pos10_A numeric
\@attribute pos10_C numeric
\@attribute pos10_T numeric
\@attribute pos10_G numeric
\@attribute pos11_A numeric
\@attribute pos11_C numeric
\@attribute pos11_T numeric
\@attribute pos11_G numeric
\@attribute pos12_A numeric
\@attribute pos12_C numeric
\@attribute pos12_T numeric
\@attribute pos12_G numeric
\@attribute pos20_A numeric
\@attribute pos20_C numeric
\@attribute pos20_T numeric
\@attribute pos20_G numeric
\@attribute pos30_A numeric
\@attribute pos30_C numeric
\@attribute pos30_T numeric
\@attribute pos30_G numeric
\@attribute method {NOMEseq,9N_PBAT,6N_PBAT,6N_sc,UMI_RRBS,Swift,WGBS,Truseq,Amplicon,RRBS,scNOME}\n\n
\@data\n";

	chdir ("./charades_temp") or die "Cannot move to temporary folder: $!";
	
	my @file_heads = <*.100K.fastq>;
	my %info;			## This is the data structure that contains all the base composition info (hash of array of hashes)


	FILE: foreach my $file_head(@file_heads) {
		my $file_name = $file_head;
		my $remove = substr($file_name,-11);
		$file_name =~ s/$remove//;


		## Read sequences in the fastq file

		open (my $in, $file_head) or die "Cannot open fastq file: $! ";


		## Collecting all the reads from the fastq file.
		my @reads;
		warn "Summarising sequence composition for $file_name\n";

		while (<$in>) {
	
			my $header = $_;
			unless ($header =~ /^@/) {
				warn "$file_head does not look like it's in fastq format. Skipping.\n";
				next FILE;
			}
			my $sequence = <$in>;
			unless ($header =~ /A|T|C|G|N/) {
				warn "Sequence in $file_head contains bases other than A, T, C, G, N. Skipping.\n";
				next FILE;
			}
			my $plusline = <$in>;
			unless ($header =~ /^@/) {
				warn "$file_head does not look like it's in fastq format. Skipping.\n";
				next FILE;
			}
			my $quality_scores = <$in>;
			last unless ($quality_scores);
			chomp $sequence;
			$sequence = uc($sequence);
			if ($sequence =~ /A|T|C|G|N/) {
				push @reads, $sequence;
				#warn "$reads[$#reads]\n";
			}
		}

		# Extracting the base composition at different positions in the read.
		# Base composition is very different between different library preps.
	
		my @interesting_positions = (1,2,3,4,8,9,10,11,12,20,30);
		
		my $As = 0;
		my $Cs = 0;
		my $Ts = 0;
		my $Gs = 0;
		my $Ns = 0;
		
	
		foreach my $interesting_position(@interesting_positions) {
			foreach my $read(@reads) {
				$read = uc($read);
				my $base = substr($read, $interesting_position -1,1);
		
				unless ($base) {
					warn "A read in $file_head is shorter than expected with only ", length($read), " bases. 
					Is this a trimmed file? \nSkipping Sequence: ", $read, "\n";
					next;
				} # There shouldn't be any shorter reads unless the runs are really old.
			
				if ($base eq "A") {
					++ $As;
				}
				elsif ($base eq "C") {
					++$Cs;
				}
				elsif ($base eq "T") {
					++$Ts;
				}
				elsif ($base eq "G") {
					++$Gs;
				}
				elsif ($base eq "N") {
					++$Ns;
				}
				else {
					warn "Unexpected base found: $base Only A,C,T,G and N allowed\n";
				}
			}
			#warn "There were $As As, $Cs Cs, $Ts Ts, $Gs, Gs and $Ns other bases at position $interesting_position\n";
			
			if ($As+$Cs+$Ts+$Gs+$Ns == 0) {
				warn "File $file_head does not seem to contain a (long enough)sequence.\n";
				next FILE;
			}
	
			my $percentA = int($As / ($As+$Cs+$Ts+$Gs+$Ns) * 100);
			my $percentC = int($Cs / ($As+$Cs+$Ts+$Gs+$Ns) * 100);
			my $percentT = int($Ts / ($As+$Cs+$Ts+$Gs+$Ns) * 100);
			my $percentG = int($Gs / ($As+$Cs+$Ts+$Gs+$Ns) * 100);
			my $percentN = int($Ns / ($As+$Cs+$Ts+$Gs+$Ns) * 100);
	
			if ($percentN > 5) {
				warn "WARNING: There are more than 5% Ns at position $interesting_position.\n";
			}
			
			warn "Percentages for A,C,T,G at position $interesting_position are $percentA, $percentC, $percentG, $percentT\n" if $debug;
			#print $out "$file_name\t$interesting_position\t$percentA\t$percentC\t$percentT\t$percentG\t$percentN\n";
			
			## Generating a data structure containing the composition info
			my $position_info = 	{position => $interesting_position,
									 percentA => $percentA,
									 percentC => $percentC,
									 percentT => $percentT,
									 percentG => $percentG};
			
			push @{$info{$file_name}},$position_info;
			
			
			## resetting base counts
			$As = 0;
			$Cs = 0;
			$Ts = 0;
			$Gs = 0;
			$Ns = 0;	
		}
		close $in or die;
	}
	
	my @temp_files = (<*.100K.fastq>);
	foreach my $temp_file(@temp_files) {
		unlink $temp_file or die "Couldn't delete file: $!";
	}
	chdir ("..") or die "Cannot move away from temporary folder: $!";
	rmdir ("charades_temp") or die "Cannot delete temporary folder: $!";
	
	foreach my $file (keys %info) {
		print $out $file, ",", $info{$file} -> [0] -> {percentA}, ",",$info{$file} -> [0] -> {percentC},",",$info{$file} -> [0] -> {percentT},",",$info{$file} -> [0] -> {percentG},
					",",$info{$file} -> [1] -> {percentA},",",$info{$file} -> [1] -> {percentC},",",$info{$file} -> [1] -> {percentT},",",$info{$file} -> [1] -> {percentG},
					",",$info{$file} -> [2] -> {percentA},",",$info{$file} -> [2] -> {percentC},",",$info{$file} -> [2] -> {percentT},",",$info{$file} -> [2] -> {percentG},
					",",$info{$file} -> [3] -> {percentA},",",$info{$file} -> [3] -> {percentC},",",$info{$file} -> [3] -> {percentT},",",$info{$file} -> [3] -> {percentG},
					",",$info{$file} -> [4] -> {percentA},",",$info{$file} -> [4] -> {percentC},",",$info{$file} -> [4] -> {percentT},",",$info{$file} -> [4] -> {percentG},
					",",$info{$file} -> [5] -> {percentA},",",$info{$file} -> [5] -> {percentC},",",$info{$file} -> [5] -> {percentT},",",$info{$file} -> [5] -> {percentG},
					",",$info{$file} -> [6] -> {percentA},",",$info{$file} -> [6] -> {percentC},",",$info{$file} -> [6] -> {percentT},",",$info{$file} -> [6] -> {percentG},
					",",$info{$file} -> [7] -> {percentA},",",$info{$file} -> [7] -> {percentC},",",$info{$file} -> [7] -> {percentT},",",$info{$file} -> [7] -> {percentG},
					",",$info{$file} -> [8] -> {percentA},",",$info{$file} -> [8] -> {percentC},",",$info{$file} -> [8] -> {percentT},",",$info{$file} -> [8] -> {percentG},
					",",$info{$file} -> [9] -> {percentA},",",$info{$file} -> [9] -> {percentC},",",$info{$file} -> [9] -> {percentT},",",$info{$file} -> [9] -> {percentG},
					",",$info{$file} -> [10] -> {percentA},",",$info{$file} -> [10] -> {percentC},",",$info{$file} -> [10] -> {percentT},",",$info{$file} -> [10] -> {percentG}, ",?\n";
		print $out_arff $info{$file} -> [0] -> {percentA}, ",",$info{$file} -> [0] -> {percentC},",",$info{$file} -> [0] -> {percentT},",",$info{$file} -> [0] -> {percentG},
					",",$info{$file} -> [1] -> {percentA},",",$info{$file} -> [1] -> {percentC},",",$info{$file} -> [1] -> {percentT},",",$info{$file} -> [1] -> {percentG},
					",",$info{$file} -> [2] -> {percentA},",",$info{$file} -> [2] -> {percentC},",",$info{$file} -> [2] -> {percentT},",",$info{$file} -> [2] -> {percentG},
					",",$info{$file} -> [3] -> {percentA},",",$info{$file} -> [3] -> {percentC},",",$info{$file} -> [3] -> {percentT},",",$info{$file} -> [3] -> {percentG},
					",",$info{$file} -> [4] -> {percentA},",",$info{$file} -> [4] -> {percentC},",",$info{$file} -> [4] -> {percentT},",",$info{$file} -> [4] -> {percentG},
					",",$info{$file} -> [5] -> {percentA},",",$info{$file} -> [5] -> {percentC},",",$info{$file} -> [5] -> {percentT},",",$info{$file} -> [5] -> {percentG},
					",",$info{$file} -> [6] -> {percentA},",",$info{$file} -> [6] -> {percentC},",",$info{$file} -> [6] -> {percentT},",",$info{$file} -> [6] -> {percentG},
					",",$info{$file} -> [7] -> {percentA},",",$info{$file} -> [7] -> {percentC},",",$info{$file} -> [7] -> {percentT},",",$info{$file} -> [7] -> {percentG},
					",",$info{$file} -> [8] -> {percentA},",",$info{$file} -> [8] -> {percentC},",",$info{$file} -> [8] -> {percentT},",",$info{$file} -> [8] -> {percentG},
					",",$info{$file} -> [9] -> {percentA},",",$info{$file} -> [9] -> {percentC},",",$info{$file} -> [9] -> {percentT},",",$info{$file} -> [9] -> {percentG},
					",",$info{$file} -> [10] -> {percentA},",",$info{$file} -> [10] -> {percentC},",",$info{$file} -> [10] -> {percentT},",",$info{$file} -> [10] -> {percentG}, ",?\n";
	}
	
	close $out or die "Could not close output file:$!\n";
	close $out_arff or die "Could not close output file:$!\n";
	
}	





sub guess_library {
	
	warn "\n\nUsing Weka Logistic Regression classifier to predict library method\nTaining data from file training_data_summary_sorted.arff\n\n";
	#system "java -cp '/bi/apps/weka/3.8.2/weka.jar' weka.classifiers.trees.J48 -t training_data_summary_sorted.arff -T sequence_composition_stats.arff -distribution -p 0 > weka.output";
	system "java -cp '/bi/apps/weka/3.8.2/weka.jar' weka.classifiers.functions.Logistic -R 3 -t training_data_summary_sorted.arff -T sequence_composition_stats.arff -distribution -p 0 > weka.output";
	
	open (my $in,"weka.output") or die "Cannot open Weka output file: $!";
	
	while(<$in>) {
		if ($_ =~ /\s+inst#/) {
			my $header = $_;
			my $i = 0;
			while (<$in>) {
				my $prediction = $_;
				warn $prediction, "\n" if $debug;
				last unless ($prediction =~ /\S+/);			## removes the empty line at the end of the file
				my (undef,undef,undef, $class, $probabilities) = split(/\s+/,$prediction);
				$class =~ s/^\d://;
				my ($NOMEseq,$PBAT_9N,$PBAT_6N,$sc_6N,$UMI_RRBS,$Swift,$WGBS,$Truseq,$Amplicon,$RRBS,$scNOME) = split(/,/,$probabilities);
							
				warn "The predicted bisulfite library for $files_to_analyse[$i] is $class \n\n",
				"Probabilities for each potential class are\n",
				"NOMEseq = ",$NOMEseq,"\n",
				"9N_PBAT = ",$PBAT_9N,"\n",
				"6N_PBAT = ",$PBAT_6N,"\n",
				"6N_sc = ",$sc_6N,"\n",
				"UMI_RRBS = ",$UMI_RRBS,"\n",
				"Swift = ",$Swift,"\n",
				"WGBS = ",$WGBS,"\n",
				"Truseq = ",$Truseq,"\n",
				"Amplicon = ",$Amplicon,"\n",
				"RRBS = ",$RRBS,"\n",
				"scNOME = ",$scNOME,"\n\n\n";
				++$i;
				}
		}
	
	
	
	
	}
}