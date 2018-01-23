#!/usr/bin/perl
use warnings;
use strict;

# This script looks at the first 100 000 sequences from a bisulfite split fastq file and then
# tries to find out which how the library was created. It will use sequence composition at different
# positions in the sequence.



# Finding out which files to analyse

my @files = <training_datasets/file_heads/*.fastq>;
#warn join("\n",@files)."\n";

my @files_to_analyse;

foreach my $file(@files) {
	$file = substr($file, 29);
	if ($file =~ /R1/) {
		push @files_to_analyse, $file;
	}
	elsif ($file !~ /R2|R3|R4/) {
		push @files_to_analyse, $file;
	}
}

warn "\nFiles chosen for analysis are\n", join("\n",@files_to_analyse),"\n\n\n";

my $As = 0;
my $Cs = 0;
my $Ts = 0;
my $Gs = 0;
my $Ns = 0;	

collect_base_composition_info();
#guess_library();



###################
### SUBROUTINES ###
###################

sub collect_base_composition_info {

	open (my $out, '>', "charades_sequence_composition.txt") or die "Cannot write to file: $!";
	print $out "file_name\tposition\t%A\t%C\t%T\t%G\t%N\n";

	FILE: foreach my $file_to_analyse(@files_to_analyse) {


		# Read sequences in the fastq file

		open (my $in, "training_datasets/file_heads/$file_to_analyse") or die "Cannot open fastq file: $! ";


		# Collecting all the reads from the fastq file.
		my @reads;

		while (1) {
	
			my $header = <$in>;
			my $sequence = <$in>;
			my $plusline = <$in>;
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
		foreach my $interesting_position(@interesting_positions) {

			foreach my $read(@reads) {
				$read = uc($read);
				my @bases = split(//,$read);
				my $base = $bases[$interesting_position -1];
				#warn "The base at position $interesting_position is an $base\n";
		
				unless ($base) {
					warn "A read in $file_to_analyse is shorter than expected with only ", length($read), " bases. 
					Is this a trimmed file?\nSequence: ", $read, "\n";
					next;
				} # There really shouldn't be any shorter reads!
			
				count_base_occurance($base);
			}
			#warn "There were $As As, $Cs Cs, $Ts Ts, $Gs, Gs and $Ns other bases at position $interesting_position\n";
			
			if ($As+$Cs+$Ts+$Gs+$Ns == 0) {
				warn "File $file_to_analyse does not seem to contain a (long enough)sequence.\n";
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
	
			warn "In file $file_to_analyse, the percentages for A,C,T,G at position $interesting_position are $percentA, $percentC, $percentG, $percentT\n";
			print $out "$file_to_analyse\t$interesting_position\t$percentA\t$percentC\t$percentT\t$percentG\t$percentN\n";
	
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




sub count_base_occurance {
	my $base = shift @_;
	
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