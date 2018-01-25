#!/usr/bin/perl
use warnings;
use strict;

# This script looks at the first 100 000 sequences from a bisulfite split fastq file and then
# tries to find out which how the library was created. It will use sequence composition at different
# positions in the sequence.

my $debug = 0;
my @files_to_analyse;

	

#files_to_analyse();
#extract_100K_reads();
collect_base_composition_info();

#guess_library();



###################
### SUBROUTINES ###
###################


sub files_to_analyse {											## Finding out which files to analyse
	
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

	my @subfolders;										## Extracting the filenames if a path is given
	system "mkdir charades_temp";
	
	warn "\nExtracting the first 100 000 reads\n";
	foreach my $file_to_analyse(@files_to_analyse) { 
	
		my $file_name;
		if ($file_to_analyse =~ /\//g) {						## remove the path from the file name
			push @subfolders, pos$file_to_analyse;
			warn "The last slash is found at ", $subfolders[-1], "\n" if $debug;
		}
		
	
		if (scalar @subfolders > 0) {							## if files are in the current directory						
			$file_name = substr($file_to_analyse, ($subfolders[-1]));
		}
		else {$file_name = $file_to_analyse};
		system "zcat $file_to_analyse | head -400000 > charades_temp/$file_name.100K.fastq";
	}
}




sub collect_base_composition_info {
	
	open (my $out, '>', "sequence_composition_stats.txt") or die "Cannot write to file: $!";
	#print $out "file_name\tpos1_A\tpos1_C\tpos1_T\tpos1_G\tpos2_A\tpos2_C\tpos2_T\tpos2_G\tpos3_A\tpos3_C\tpos3_T\tpos3_G\tpos4_A\tpos4_C\tpos4_T\tpos4_G\tpos8_A\tpos8_C\tpos8_T\tpos8_G\tpos9_A\tpos9_C\tpos9_T\tpos9_G\tpos10_A\tpos10_C\tpos10_T\tpos10_G	pos11_A	pos11_C\tpos11_T\tpos11_G\tpos12_A\tpos12_C\tpos12_T\tpos12_G\tpos20_A\tpos20_C\tpos20_T\tpos20_G\tpos30_A\tpos30_C\tpos30_T\tpos30_G\n";
	print $out "file_name\tposition\tpercentA\tpercentC\tpercentT\tpercentG\tpercentN\n";
	
	chdir ("./charades_temp") or die "Cannot move to temporary folder: $!";
	
	my @file_heads = <*.100K.fastq>;

	FILE: foreach my $file_head(@file_heads) {
		my $file_name = $file_head;
		my $remove = substr($file_name,-11);
		$file_name =~ s/$remove//;


		## Read sequences in the fastq file

		open (my $in, $file_head) or die "Cannot open fastq file: $! ";


		## Collecting all the reads from the fastq file.
		my @reads;
		warn "\nSummarising sequence composition for $file_name\n\n";

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
				my @bases = split(//,$read);
				my $base = $bases[$interesting_position -1];
				#warn "The base at position $interesting_position is a $base\n";
		
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
			
			warn "Percentages for A,C,T,G at position $interesting_position are $percentA, $percentC, $percentG, $percentT\n";
			print $out "$file_name\t$interesting_position\t$percentA\t$percentC\t$percentT\t$percentG\t$percentN\n";
			
			## Generating a data structure containing the composition info
			
			
			
	
			$As = 0;
			$Cs = 0;
			$Ts = 0;
			$Gs = 0;
			$Ns = 0;	
		}
		close $in or die;
	}
	close $out or die "Could not close output file:$!\n";
}




sub guess_library {
	
	open (my $in, "training_data_summary.txt") or die "Couldn't open summary file: $!";
	my $header_line = <$in>;
	my %info;
	
	## Create a data structure that contains all the info (hash of array of hashes)
	while (<$in>) {
		chomp;
		my @sequence_composition = split (/\t/);
		my ($name,$pos,$percentA,$percentC,$percentT,$percentG) = @sequence_composition[0,1,2,3,4,5];
		#warn "$name\t$pos\t$percentA\t$percentC\t$percentT\t$percentG\n";
		
		my $position_info = 		{position => $pos,
									 percentA => $percentA,
									 percentC => $percentC,
									 percentT => $percentT,
									 percentG => $percentG};
									 
		push @{$info{$name}},$position_info;							 
		
	}
	
	#foreach my $file (keys %info) {
		#warn "\n$file\n\tposition ", $info{$file} -> [0] -> {position}, "\tpercentA ", $info{$file} -> [0] -> {percentA}, "\n";
		#warn "\tposition ", $info{$file} -> [2] -> {position}, "\tpercentA ", $info{$file} -> [2] -> {percentA}, "\n";
		#warn "\tposition ", $info{$file} -> [5] -> {position}, "\tpercentA ", $info{$file} -> [5] -> {percentA}, "\n";
	#}
	
	## Guess which library from sequence composition
	
	foreach my $file (keys %info) {
	
		if ($info{$file}->[7]->{percentC} <= 8) {	   				# percent C at position 30 (this goes into the WGBS branch) 
			if ($info{$file}->[4]->{percentC} <= 0) {   			# percent C at position 11
				if ($info{$file}->[1]->{percentC} <= 5) {			# percent C at position 8
					if ($info{$file}->[3]->{percentG} <= 26) {		# percent G at position 10
						warn "$file seems to be WGBS\n";			# leaf WGBS
					}
					if ($info{$file}->[3]->{percentG} > 26) {		# percent G at position 10
						warn "$file seems to be Amplicon\n";		# leaf Amplicon
					}
				}
				if ($info{$file}->[1]->{percentC} > 5) {			# percent C at position 8
					warn "$file seems to be UMI-RRBS\n";			# leaf UMI-RRBS
				}
			}
			if ($info{$file}->[4]->{percentC} > 0) {				# percent C at position 11
				if ($info{$file}->[3]->{percentC} <= 2) {			# percent C at position 10
					if ($info{$file}->[0]->{percentG} <= 29) {		# percent G at position 4
						warn "$file seems to be Swift\n";			# leaft Swift
					}
					if ($info{$file}->[0]->{percentG} > 29) {		# percent G at position 4
						warn "$file seems to be Truseq\n";			# leaf Truseq
					}
				}
				if ($info{$file}->[3]->{percentC} > 2) {			# percent C at position 10
					warn "$file seems to be NOMEseq\n";				# leaf NOMEseq
				}
			}
		}
		if ($info{$file}->[7]->{percentC} > 8) {	   				# percent C at position 30 (this goes into the PBAT branch) 
			if ($info{$file}->[1]->{percentC} <= 18) {				# percent C at position 8
				if ($info{$file}->[1]->{percentG} <= 12) {			# percent G at position 8
					if ($info{$file}->[0]->{percentC} <= 29) {		# percent C at position 4
						warn "$file seems to be scNOME\n";			# leaf scNOME
					}
					if ($info{$file}->[0]->{percentC} > 29) {		# percent C at position 4
						if ($info{$file}->[0]->{percentT} <= 11) {	# percent T at position 4
							warn "$file seems to be scNOME\n";		# leaf scNOME
						}
						if ($info{$file}->[0]->{percentT} > 11) {	# percent T at position 4
							warn "$file seems to be 9N_sc\n";		# leaf 9N_sc
						}
					}
				}
				if ($info{$file}->[1]->{percentG} > 12) {			# percent G at position 8
					if ($info{$file}->[0]->{percentT} <= 17) {		# percent T at position 4
						warn "$file seems to be scNOME\n";			# leaf scNOME
					}
					if ($info{$file}->[0]->{percentT} > 17) {		# percent T at position 4
						warn "$file seems to be 6N_sc\n";			# leaf 6N_sc
					}
				}
			}
			if ($info{$file}->[1]->{percentC} > 18) {				# percent C at position 8
				if ($info{$file}->[1]->{percentG} <= 1) {			# percent G at position 8
					warn "$file seems to be 6N_PBAT\n";				# leaf 9N_PBAT
				}
				if ($info{$file}->[1]->{percentG} > 1) {			# percent G at position 8
					warn "$file seems to be 9N_PBAT\n";				# leaf 9N_PBAT
				}
			}
		}
	}
	
	
	
	
	
	
	
	
	
	
	close $in or die $!;
	warn "\n\n";
}