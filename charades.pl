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
my ($output_dir, $parent_dir, $project,$r_path) = process_commandline();


files_to_analyse();
extract_100K_reads();
collect_base_composition_info();
guess_library();
draw_graph();



###################
### SUBROUTINES ###
###################

sub process_commandline {

	my $help;
	my $output_dir;
	my $project;
	my $r_path = "R";
	
	my $command_line =  GetOptions ('help'                => \$help,
									'o|output_dir=s'      => \$output_dir,
									'p|project=s'		  => \$project,
									'r_path=s'			  => \$r_path,
									);
	
	
	die "Please respecify command line options\n\n" unless ($command_line);
	

    if ($help){
	print_helpfile();
	exit;
    }
		
	die "\nPlease specifiy gzipped fastq files\n\n" unless (@ARGV);
	
	### PARENT DIRECTORY
    my $parent_dir = getcwd();
    unless ($parent_dir =~ /\/$/){    ## making parent directory end in a slash
		$parent_dir =~ s/$/\//;
    }

    ### OUTPUT DIRECTORY
    if (defined $output_dir){
		unless ($output_dir eq ''){
			unless ($output_dir =~ /\/$/){     ## making output directory end in a slash
			$output_dir =~ s/$/\//;
			}
	    
	    if (chdir $output_dir){
			$output_dir = getcwd(); #  making the path absolute
			unless ($output_dir =~ /\/$/){
				$output_dir =~ s/$/\//;
			}
	    }
	    
		else {
			mkdir $output_dir or die "Unable to create directory $output_dir $!\n";
			warn "\nCreating output directory $output_dir \n\n"; 
			chdir $output_dir or die "Failed to move to $output_dir\n";
			$output_dir = getcwd(); #  making the path absolute
			unless ($output_dir =~ /\/$/){
				$output_dir =~ s/$/\//;
			}
	    }
	    warn "\nOutput will be written into the directory: $output_dir\n";
		}
    }
    else {
		$output_dir = '';
    }

    # Changing back to parent directory
    chdir $parent_dir or die "Failed to move to $parent_dir\n";
	
	
	### PROJECT NAME
	if (defined $project) {
		$project = $project.".";
	}
	else {
		warn "\nNo project name defined, using generic output file names.\n";
		$project = '';
	}
	
	
	return ($output_dir, $parent_dir, $project,$r_path);
	
	
	
}

sub print_helpfile{
	print "\n\nCharades will try and predict which bisulfite sequencing library preparation was used. It uses gzipped fastq files as input.\n\n";
	print ">>> USAGE: ./charades.pl [options] filename(s) <<<\n\n";
	print "\nOptions:\n\n";
	print "--output_dir / -o [path]\tOutput directory, either relative or absolute. If not specified, output is written to the current directory.\n";
	print "--project / -p [project name]\tSpecifies a project name for all tested samples which will be added to the output file names.\n\n\n";
}


sub files_to_analyse {											

## Finding out which files to analyse

	my @files = @ARGV;
	unless (@files) {
		print "Please provide one or more fastq.gz files for analysis\n\n";
	}
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
	@files_to_analyse = sort(@files_to_analyse);
	warn "\nFiles selected for analysis are\n", join("\n",@files_to_analyse),"\n\n\n";
}




sub extract_100K_reads {

	my @subfolders;										
	mkdir "${output_dir}charades_temp" or die "Could not create temporary directory: $!";
	
	warn "\nExtracting the first 100 000 reads from fastq files\n\n";
	my $file_name;

	foreach my $file_to_analyse(@files_to_analyse) { 
	
		$file_name = $file_to_analyse;
		$file_name =~ s/^.*\///;
		system "zcat $file_to_analyse | head -400000 > ${output_dir}charades_temp/$file_name.100K.fastq";
	}
}




sub collect_base_composition_info {
	
	open (my $out, '>', "${output_dir}${project}sequence_composition_stats.csv") or die "Cannot write to file: $!";
	open (my $out_arff, '>', "${output_dir}${project}sequence_composition_stats.arff") or die "Cannot write to file: $!";
	print $out "file_name,pos1_A,pos1_C,pos1_T,pos1_G,pos2_A,pos2_C,pos2_T,pos2_G,pos3_A,pos3_C,pos3_T,pos3_G,pos4_A,pos4_C,pos4_T,pos4_G,pos8_A,pos8_C,pos8_T,",
				"pos8_G,pos9_A,pos9_C,pos9_T,pos9_G,pos10_A,pos10_C,pos10_T,pos10_G,pos11_A,pos11_C,pos11_T,pos11_G,pos12_A,pos12_C,pos12_T,pos12_G,",
				"pos20_A,pos20_C,pos20_T,pos20_G,pos30_A,pos30_C,pos30_T,pos30_G,method\n";
	#print $out "file_name\tposition\tpercentA\tpercentC\tpercentT\tpercentG\tpercentN\n";
	
	print $out_arff "\@relation base_signature\n
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
\@attribute method {Swift,6N_PBAT,6N_sc,9N_sc,WGBS,RRBS,scNOME,UMI_RRBS,9N_PBAT,Amplicon,NOMEseq,Truseq}\n\n
\@data\n";

	chdir "${output_dir}charades_temp" or die "Cannot move to temporary folder: $!";
	
	my @file_heads = <*.100K.fastq>;
	my %info;			## This is the data structure that contains all the base composition info (hash of array of hashes)

	unless (@file_heads) {
		die "Cannot find files for first 100 000 reads.Exiting.\n";
	}
	
	FILE: foreach my $file_head(@file_heads) {
		my $file_name = $file_head;
		warn $file_name, "\n" if $debug;
		my $remove = substr($file_name,-11);
		$file_name =~ s/$remove//;
		warn $file_name, "\n" if $debug;


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
			unless ($sequence =~ /^[A|T|C|G|N]*$/) {
				warn "Sequence in $file_head contains bases other than A, T, C, G, N. Skipping.\n";
				next FILE;
			}
			my $plusline = <$in>;
			unless ($plusline =~ /^\+/) {
				warn "$file_head does not look like it's in fastq format. Skipping.\n";
				next FILE;
			}
			my $quality_scores = <$in>;
			last unless ($quality_scores);
			chomp $sequence;
			$sequence = uc($sequence);
			if ($sequence =~ /^[A|T|C|G|N]*$/) {
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
	chdir $parent_dir or die "Cannot move away from temporary folder: $!";
	rmdir ("${output_dir}charades_temp") or die "Cannot delete temporary folder: $!";
	
	foreach my $file (sort keys %info) {
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
	
	warn "\n\nUsing Weka Logistic Regression classifier to predict library method\nTaining data from file training_data_20181211.arff\n\n";
	unless (-e "training_data_20181211.arff") {
		die "Couldn't find file >> training_data_20181211.arff <<. Please add it to the same directory as charades.pl.\n\n\n";
	}
	#system "java -cp '/bi/apps/weka/3.8.2/weka.jar' weka.classifiers.trees.J48 -t training_data_20181210.arff -T sequence_composition_stats.arff -distribution -p 0 > ${output_dir}weka.output";
	system "java -cp '/bi/apps/weka/3.8.2/weka.jar' weka.classifiers.functions.Logistic -R 1 -t training_data_20181211.arff -T ${output_dir}${project}sequence_composition_stats.arff -distribution -p 0 > ${output_dir}weka.output";
	
	open (my $in,"${output_dir}weka.output") or die "Cannot open Weka output file: $!";
	open (my $out, '>', "${output_dir}${project}class.probabilities.txt") or die "Cannot create class probability file: $!";
	
	## making the header for the class probability file which is then to be used for the graphical representation
	print $out "sample\tSwift\tPBAT_6N\tsc_6N\tsc_9N\tWGBS\tRRBS\tscNOME\tUMI_RRBS\tPBAT_9N\tAmplicon\tNOMEseq\tTruseq\n";
	
	while(<$in>) {
		unless ($_) {
			warn "No sequence composition data was collected. Exiting.\n";
		}
		if ($_ =~ /\s+inst#/) {
			my $header = $_;
			my $i = 0;
			my $sample;
			
			while (<$in>) {
				my $prediction = $_;
				warn $prediction, "\n" if $debug;
				last unless ($prediction =~ /\S+/);			## removes the empty line at the end of the file
				my (undef,undef,undef, $class, $probabilities) = split(/\s+/,$prediction);
				$class =~ s/^\d+://;
				my ($Swift,$PBAT_6N,$sc_6N,$sc_9N,$WGBS,$RRBS,$scNOME,$UMI_RRBS,$PBAT_9N,$Amplicon,$NOMEseq,$Truseq) = split(/,/,$probabilities);
				#{Swift,6N_PBAT,6N_sc,9N_sc,WGBS,RRBS,scNOME,UMI_RRBS,9N_PBAT,Amplicon,NOMEseq,Truseq}		
				# THIS NEEDS TO BE THE SAME ORDER AS THE CLASS ATTRIBUTE IN THE ARFF FILE!!
				
				# removing the star from the predicted class to make them all numeric
				foreach my $probability($Swift,$PBAT_6N,$sc_6N,$sc_9N,$WGBS,$RRBS,$scNOME,$UMI_RRBS,$PBAT_9N,$Amplicon,$NOMEseq,$Truseq) {
					$probability =~ s/\*//;
				}
				
				
				# remove the path from the file name
				$sample = $files_to_analyse[$i];
				$sample =~ s/^.*\///;
				
				warn "The predicted bisulfite library for $sample is $class \n";
				print $out "$sample\t$Swift\t$PBAT_6N\t$sc_6N\t$sc_9N\t$WGBS\t$RRBS\t$scNOME\t$UMI_RRBS\t$PBAT_9N\t$Amplicon\t$NOMEseq\t$Truseq\n";

				### The commented out sections print the probabilities for each class to STOUT, either as such or from a hash. I'll leave these in
				### for now to see what my preferred option is to process the information downstream
					
				#warn "The predicted bisulfite library for $files_to_analyse[$i] is $class \n\n",
				#"Probabilities for each potential class are\n",
				#"NOMEseq = ",$NOMEseq,"\n",
				#"9N_PBAT = ",$PBAT_9N,"\n",
				#"6N_PBAT = ",$PBAT_6N,"\n",
				#"6N_sc = ",$sc_6N,"\n",
				#"UMI_RRBS = ",$UMI_RRBS,"\n",
				#"Swift = ",$Swift,"\n",
				#"WGBS = ",$WGBS,"\n",
				#"Truseq = ",$Truseq,"\n",
				#"Amplicon = ",$Amplicon,"\n",
				#"RRBS = ",$RRBS,"\n",
				#"scNOME = ",$scNOME,"\n\n\n";
				++$i;
				
				## Collect top two class probabilities
				
				
				#my %probs = (
				#				NOMEseq => $NOMEseq,
				#				PBAT_9N => $PBAT_9N,
				#				PBAT_6N => $PBAT_6N,
				#				sc_6N => $sc_6N,
				#				UMI_RRBS => $UMI_RRBS,
				#				Swift => $Swift,
				#				WGBS => $WGBS,
				#				Truseq => $Truseq,
				#				Amplicon => $Amplicon,
				#				RRBS => $RRBS,
				#				scNOME => $scNOME,
				#				);
								
				
				#foreach my $method (sort {$probs{$b} <=> $probs{$a}} keys %probs) {
				#	print "$method\t $probs{$method}\n";
				#}
				#print "\n\n";				
					
			}
			
		}
	}
	close $in;
	close $out or die "Could not close filehandle for class probabilities";
}



sub draw_graph {

	my $infile = "${output_dir}${project}class.probabilities.txt";
	my $graph = "${project}probabilities.png";
	
	if ($output_dir eq "") {
		$output_dir = getcwd();
	}
	
	my $r_script = <<"END_SCRIPT";
	
		directory <- "$output_dir";
		file <- "$infile";
		graph <- "$graph";
		
		setwd(directory)

		read.delim(file,stringsAsFactors = FALSE) -> class.probabilities
		row.names(class.probabilities) <- class.probabilities[,1]
		class.probabilities[2:ncol(class.probabilities)] -> class.probabilities
		substr(row.names(class.probabilities),1,20) -> row.names(class.probabilities)
		as.matrix(t(class.probabilities)) -> class.probabilities

		library(RColorBrewer)
		cols = brewer.pal(12,"Paired")
		
		png(filename = "$graph", width = 800, height = max(400,60 + (20*ncol(class.probabilities))))

		par(oma = c(1,1,1,9), mar = c(5,12,4,3))

		barplot(class.probabilities,
			horiz = T,
			col = cols,
			main = "Probabilities for each method",
			las = 1,
			xlab = "probability"
			)

		par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
		plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

		legend ("topright",
			legend = row.names(class.probabilities),
			fill = cols,
			xpd = T,
			inset = c(0.08,0.2)
			)

		dev.off()
		
END_SCRIPT
		
		warn "The path to R is $r_path\n" if $debug;
		
		open (my $in, "| \"$r_path\" --no-save > /dev/null 2>&1") or die "Can't open pipe to '$r_path'  - do you need to set an rpath?";

		print $in $r_script;

		close $in or die "Can't write to '$r_path' pipe - Do your file names start with 20 common characters?";

		
}











