#!/usr/local/bin/perl 
#
#require "./masses_and_fragments.pl";

# reads each file in the directory - either the first argument on the command line, or the current directory
#for each fasta file in the directory, the name/sequence is read in and then the protein is digested with 
#trypsin to get the resulting peptides - 
#a single output file is created - "all_predigested.fasta" that contains, 
#for the description line: the peptide sequence, followed by a number representing the number of proteins that 
#  resulted in that peptide when digested with trypsin
#for the sequence line: the peptide sequence

use warnings;
use strict;
#use JSON; # imports encode_json, decode_json, to_json and from_json.

my $dir;
my $incompletes;
my $use_tail;

if ($ARGV[0] && $ARGV[0]=~/\w/) { $dir=$ARGV[0];}  
else { $dir = "/Users/sarahkeegan/fenyolab/data_and_results/Llama_MiSeq/MISEQ_DB_2017_06_12/2015-5094-H/trimmed/merged/fasta/removed_primers/translated/"; } 

if ($ARGV[1] && $ARGV[1]=~/\w/) { $incompletes=$ARGV[1];}
else { $incompletes=1; }

if ($ARGV[2] && $ARGV[2]=~/\w/) { $use_tail=$ARGV[2];}
else { $use_tail=0; }

my %protein_peptides = ();

my $total_pep_count=0;
my %PEP=();
my %PEP_proteins=();
my %PEP_proteins_count=();
my $files_count=0;
my $peptides_count=0;
my $proteins_count=0;
my %proteins_count=();
my $line="";
my @TAIL_SEQUENCES = ("SEPKIPQPQPKPQ", "SAHHSEDPSSKCP", "SEPKTPKPQPQPQPQ", "SGTNEVCKCPKCPAPEL", "EPKIPQPQPKPQ", "AHHSEDPSSKCP", "EPKTPKPQPQPQPQ", "GTNEVCKCPKCPAPEL");

if (!opendir(DIR,"$dir")) { print "Error reading $dir\n"; exit(1); }

my @allfiles=readdir DIR;
closedir DIR;
my $count_fasta = 0;

if (!open (IND_OUT, ">$dir/protein_peptides.fasta"))
{
	print "Error creating $dir/protein_peptides.fasta\n";
	exit(1);
	
}

foreach my $filename (@allfiles)
{
	if ($filename!~/\.fas?t?a?$/i) { next; } #can be *.fa, *.fas, *.fast, *.fasta, case insensitive
	
	if (!open (IN,"$dir/$filename")) { print "Error opening $dir/$filename.\n"; next; }
	$count_fasta++;
	#print qq!$filename\n!;
	my $name="";
	my $sequence="";
	while ($line=<IN>)
	{
		chomp($line);
		if ($line =~ s/^>//)
		{#if the current line starts with a '>', then it is the description line, remove the '>' and input the 'name'
			
			my $name_=$line;
			if ($name_ =~ s/^\"//)
			{#if the line begins with quotes (")
			 #remove quotes at beginning and end
			 #replace any non-word character (not letters or numbers) with '_'
				$name_ =~ s/\"$//;
				$name_ =~ s/[^\w]/_/g; 
			}
			else
			{#line does not begin with quotes ("), remove anything after (and including) the first whitespace character 
				$name_ =~ s/\s.*$//;
			}
			if ($name =~ /\w/ and $sequence =~ /\w/)
			{#both name and sequence were inputted, add the protein to the count
				
				
			
				#digest the protein with trypsin, the function fills the %PEP hash with the resulting peptides
				%PEP=();
				
				#add tail sequence to protein so that we don't miss c terminus peptides:
				if($use_tail)
				{
					my $new_sequence;
					foreach my $tail (@TAIL_SEQUENCES)
					{
						$new_sequence = $sequence . $tail;
						DigestTrypsin($name,$new_sequence,$incompletes);
					}
				}
				else
				{
					DigestTrypsin($name,$sequence,$incompletes);
				}
				
				if($proteins_count % 50000 == 0) { print "$proteins_count\n"; } #if($protein_count >= 10000) { last; } }
				#if ($proteins_count > 1000) {
				#	last;
				#}
				
				
				#strip primer info from name before outputting!
				$name =~ s/(_fr[012])_fwd_(\w*)_rev_(\w*)$//;
				$name = $name . $1;
				
				print IND_OUT ">$name\n";
				foreach my $peptide (keys %PEP)
				{#we want to keep a count of the number of proteins that contained a particular peptide when digested 
				 #with trypsin - a list of each 'name' of the proteins is stored in the hash %PEP_proteins and the number 
				 #of proteins for a given peptide is stored in the hash %PEP_proteins_count
					if (length($peptide)>6)
					{#disregard peptides that are too short
						#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
						#{
						#	$PEP_proteins_count{$peptide}++;
						#	$PEP_proteins{$peptide} .= qq!#$name#!;
						#}
						#else
						#{
						#	print "$name $peptide\n";
						#}
						$PEP_proteins_count{$peptide}++;
						print IND_OUT "#$peptide#";
					}
				}
				print IND_OUT "#\n";
				
				$proteins_count++;
				$proteins_count{$filename}++;
				#print "$proteins_count. $name\n";
				
			}
			$name=$name_;
			$sequence="";
		}
		else
		{#it is a line of the sequence
			$sequence.="$line";
		}
	}	
	
	#the last protein left in the file, do the same as above
	if ($name =~ /\w/ and $sequence =~ /\w/)
	{
		#digest the protein with trypsin, the function fills the %PEP hash with the resulting peptides
		%PEP=();
		
		#add tail sequence to protein so that we don't miss c terminus peptides:
		if($use_tail)
		{
			my $new_sequence;
			foreach my $tail (@TAIL_SEQUENCES)
			{
				$new_sequence = $sequence . $tail;
				DigestTrypsin($name,$new_sequence,$incompletes);
			}
		}
		else
		{
			DigestTrypsin($name,$sequence,$incompletes);
		}
		
		#strip primer info from name before outputting!
		$name =~ s/(_fr[012])_fwd_(\w*)_rev_(\w*)$//;
		$name = $name . $1;
		
		print IND_OUT ">$name\n";
		foreach my $peptide (keys %PEP)
		{
			#print qq!$peptide\n!;
			if (length($peptide)>6)
			{
				#if (!$PEP_proteins{$peptide} || ($PEP_proteins{$peptide} !~ /#$name#/))
				#{
				#	$PEP_proteins_count{$peptide}++;
				#	$PEP_proteins{$peptide} .= qq!#$name#!;
				#}
				$PEP_proteins_count{$peptide}++;
				print IND_OUT "#$peptide#";
			}
		}
		print IND_OUT "#\n";
		
		$proteins_count++;
		$proteins_count{$filename}++;
		#print "$proteins_count. $name\n";
	}
	close(IN);
	$files_count++;
	print qq!$files_count. $filename $proteins_count{$filename}\n!;
}
close(IND_OUT);

if ($count_fasta == 0)
{
	print "Warning: No fasta files processed!\n";
	exit(1);
}

if (open (OUT,">$dir/all_predigested.fasta"))
{
	foreach my $peptide (keys %PEP_proteins_count) #(keys %PEP_proteins)
	{
		print OUT qq!>$peptide $PEP_proteins_count{$peptide}\n$peptide\n!;
	}
	close(OUT);
}
else
{
	print "Error creating $dir/all_predigested.fasta\n";
	exit(1);
}

##save to JSON the protein peptides index
#my $utf8_encoded_json_text = encode_json \%protein_peptides;
#if (open (IND_OUT, ">$dir/protein_peptides.json"))
#{
#	print IND_OUT $utf8_encoded_json_text;
#	close(IND_OUT);
#}

exit(0);

sub DigestTrypsin
{#fills the %PEP hash with the peptides resulting from digesting the sequence with trypsin
 #arg1 - name, arg2 = sequence, arg3 = # of incompletes
	my $name = shift();
	my $seq = shift();
	my $incompletes = shift();

	my $temp=$seq;
	my @pep=();
	my @start=();
	my @end=();
	my $aa="";
	my $aa_="";
	my $i=0;

	for($i=0;$i<=$incompletes;$i++)
	{
		$start[$i]=0;
		$end[$i]=-1;
		#$pep[$i]="[";
	}
	my $aa_count=0;
	while ($temp =~ s/^\s*([A-Z\*])//)
	{
		$aa="\U$1";
		$aa=~s/I/L/g;
		if ( (($aa_=~/R/ or $aa_=~/K/) and $aa!~/P/) or $aa_=~/\*/)
		{
			for($i=0;$i<=$incompletes;$i++)
			{
				$PEP{"$pep[$i]"}=1;
				$pep[$i]=$pep[$i+1];
				$start[$i]=$start[$i+1];
				$end[$i]=$end[$i+1];
			}
			$start[$incompletes]=$aa_count;
			$end[$incompletes]=$aa_count-1;
		}
		for($i=0;$i<=$incompletes;$i++)
		{
			if ($aa!~/\*/) { $pep[$i].=$aa; }
			$end[$i]++;
		}
		$aa_=$aa;
		$aa_count++;
	}
	for($i=0;$i<=$incompletes;$i++)
	{
		$PEP{"$pep[$i]"}=1;
	}
}


