#!/usr/local/bin/perl 

use warnings;
use strict;
use JSON; # imports encode_json, decode_json, to_json and from_json.


my $dir="";
my $out_dir = "";

if ($ARGV[0] && $ARGV[0]=~/\w/) { $dir="$ARGV[0]";}
else { $dir = "."; } 

if ($ARGV[1] && $ARGV[1]=~/\w/) { $out_dir="$ARGV[1]";}
else { $out_dir = "."; }

#my %ctla4_hits_5094_7_22 = ('M00587_82_000000000_AA4PE_1_2118_15018_24150_1' => 'GGGAGGCTTGGTGCAGGCTGGGGGGTCTCTGAGACTCTCCTGTGCGGCCTCTGGAAGTGGTTTCTTTATTCATGCCATGGCCTGGTACCGCCAAACTCCAGGGAAGCAGCGCGAGTTGGTCGCAGAGCTTACTAGCGAGGATAGCACAACCTATGCAGACTCCGTGAAGGGCCGAGTCACCATCTCCAGAGACGGGAACTTCAAGAACACGGCAAATCTGCAAATGAACAACCTGAAACCAGAGGACACGGCCGTTTATTACTGTGCAGCAGCACGAACACTAGTTTTTAGTAGGCCCCGCTATGACTACTGGGGCCAGGGG',
#                        'M00587_82_000000000_AA4PE_1_2113_7222_14022_1' => 'GGGAGGATTGGTGCAGGCTGGGGGCTCTCTGAGAATTTCCTGTGCAGCCGCTGGATTCAGCTTCAGTACTGCTGGCACGGCCTGGTTCCGCCAGGCTCCAGGGAAGGAGCGTGAGTTTGTGGGAGGTATCAACTGGAGTGGTGATAGTACAAGGCTTGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGAGAGTACCAAGAACGTGGTGTATCTGCAAATGAACAGCCTGAAACCTGAGGACACGGCCGTTTACTACTGTGCAGCAAACCTGAATAGCGACTATACCCCGAGTAGATGGACTTACTGGGGCCAGGGG',
#                        'M00587_82_000000000_AA4PE_1_1119_14834_18252_1' => 'CCCTTTGCCCCAGTTATCCATGCCAGTCCTCAGGTCTGACGCACAGTAATAGCGGGCCGTGTCCTCAGGTTTCAGGCTGTTCATTTGCAGAAACACCGTGTTCTTGGAGTTGTCTCTGGAGATGGTGAATCGGCCCTTCACGGAGTCTGCATAGTATGTGGTGCCACCACCCGTACTAATACTAGAGACCCACTCGAGCCCCTTCCCTGGAGCCTGGCGGACCCAGGTCATGTAGTAGCTACCGAAGGGGAATCCAGAGGCTGTACAGGAGAGTCTTAGAGACCCCCCAGGCTGCACCGAGCCTCCC',
#                        'M00587_82_000000000_AA4PE_1_1118_21917_19811_1' => 'GGGAGTCTTGGTGCAGCCTGGGGGGTCTCTGAGACTCTCCTGTGCAGCCTCTGGTTTCAGTTTTGATGCTTTTGACATAAGCTGGTTCCGTCAGGCCCCAGGGAAGGAGCGCGAGGGGGTCTCATGTATTGATAGTAGTGGTACTTTCACATACTATACAGGCTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGGTGTATCTGCAAATGAACAGCCTGAAACCCGAAGACACGGCCGTTTATTACTGTGCAACACCTCGGGCGCGTTGCGATAGTGACCGGTACGAGTACTCTGACTATTCCAACCGGGGCCAGGGG',
#                        'M00587_82_000000000_AA4PE_1_1117_24349_23082_1' => 'GGGAGGCTTCGTGCAGACTGGGGGGTCTCTGAGACTCTCTTGTGCAGCATCTGGATTCGTGTTCAATGATTATACCATGGGCTGGTACCGCCAGGCTCCAGGGAAGCAGCGTGAGTTGGTCGCAACTATAAATACGCGTGGCGTCACATGGTATTCGAGCGCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAATACGGTGTATCTGCAATTGAACAGCCTGAAACCTGATGACACGGCCGTCTATTATTGTAACGCACCATTTAAGGACTTTGGTTCTTGGGGCCAGGGG',
#                        'M00587_82_000000000_AA4PE_1_2103_9372_1708_1' => 'GGGAGGCTTGGTGCAGGCTGGGGGGTCTCTGAGACTCTCCTGTGTGGCCTCTGGAAGCGGCTTCTTTATTTATGCCATGGGCTGGTACCGCCAAACTCCAGGGAAGCAGCGCGAGTTGGTCGCAGGTCTTACTAGCGATGATAGCACAAGCTATGCAGACTCCGTGAAGGGCCGATACACCGTCTCCAGAGAGAAGAACTTCAAGAACACGGTATATCTGCAAATGAACAACCTGAAACCTGAGGACACGGCCGTTTATTACTGTGCAGCGTCAGGAACACTCGAATTTGGCAGTCCCCGTTATGACTACTGGGGCCAGGGG');
#my %ctla4_hits_5094_7_17 = ('M00587_82_000000000_AA4PE_1_2102_16222_8804_1' => 'CCCCTGGCCCCGGTAAGTCCATCTACTCGGGGTATAGTCGCTATTTAGGTTTGCTGCACAAGAGTAAACGGCCGTGTCCTCAGGTTTCAGGCTGTTCATTTCAAGATACACCGTGTTCTTGGCACTGTCTCTGGAGAGGGTGAATCGGCCCTTCACGGAGTCTGCAACACTTGTAGTATCACCACTCCAACTGATACCTCCCACAAACTCACGCTCCTTCCCTGGAACCTGGCGGTACCAGGCCGTGCCAGCGGTACTTAAGGTGAATCCAGAGGCTGTACAGGAGAGTCTCAGAGAGCCCCCAGCCTGCACCAATCCTCCC');

my %redundancy_listing;
my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
		"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
		"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
		"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
		
		"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
		"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
		"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
		"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
		
		"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
		"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
		"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
		"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
		
		"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
		"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
		"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
		"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");

print "Executing dnatoprot_longest_nr.pl\n";

if (!opendir(DIR,"$dir")) { print "Error reading $dir ($!)\n"; exit(2); }

my %longest_orf=();
my %orf=();
my %orf_trim=();

my @allfiles=readdir DIR;
closedir DIR;
my $count_fasta = 0;
foreach my $filename (@allfiles)
{#for each fasta file  or fastq file in the directory:
	my $ftype;
	if ($filename =~ /\.fas?t?a?$/i) { $ftype = 'FASTA'; }
	elsif($filename =~ /\.fastq$/i || $filename =~ /\.fq$/i) { $ftype = 'FASTQ'; }
	else { next; }
	
	#print "$filename\n";
	if (!open (IN,"$dir/$filename")) { print "Error opening $dir/$filename ($!).\n"; next; }
	print "Opened $dir/$filename\n";
	$count_fasta++;
	my $name="";
	my $description="";
	my $sequence="";
	my $line="";
	my $entry_count = 0;
	while ($line=<IN>)
	{
		chomp($line);
		if ($line=~/^[>@](\S+)\s?(.*)$/)
		{#if the current line is the description line - begins with '>', and then 
		 #one or more non-whitespace, followed by 0 or more whitespace char's, then anything until the end
			$entry_count++;
			#if ($entry_count % 1000 == 0) { print "$entry_count\n"; }
			#if ($entry_count > 1000) { last; }
			
			my $name_=$1; 
			my $description_=$2;
			$name_ =~ s/[\:\-\\\/]/_/g;
			if ($name=~/\w/ and $sequence=~/\w/)
			{#the entire sequence has been read in, so do the conversion:
				my $size = length($sequence);
				my $longest_orf_direction="";
				my $longest_orf_frame="";
				my $longest_orf_length=0;
				my $longest_orf_seq="";
				foreach my $direction ("fwd","rev")
				{#read both forward and reverse
					my $seq="";
					if ($direction =~ /^fwd$/) { $seq=$sequence; } 
					else 
					{ 
						$seq = reverse $sequence;
						$seq =~ tr/ATCG/TAGC/;
					}
					for (my $k=0;$k<3;$k++)
					{#for each reading frame:
						my $protein="";
						for(my $n=$k;$n<$size;$n=$n+3)
						{
							my $triplet = substr($seq, $n, 3);
							if ($mapping{$triplet} and $mapping{$triplet} =~ /[\w\*]/) { $protein .= $mapping{$triplet}; } # '*' is stop codon
							else { $protein.="X"; } # X is unknown, doesn't code for anything must be error in sequence
						}
						#remove X's at beginning and end
						$protein =~ s/X+$//;
						$protein =~ s/^X+//;
				
						my $temp="$protein*";
						my $index=0;
						
						#find all the orf's - the sequence up to the next stop codon
						while ($temp =~ s/^([^\*]*)\*//) 
						{# $1 is 0 or more of any char's except '*', followed by '*', starting at the beginning
						 #so we get the AA sequence up to the next stop codon, and remove it
							my $orf_seq=$1;
							
							#remove X at beginning and end
							$orf_seq =~ s/X+$//;
							$orf_seq =~ s/^X+//;
							
							my $orf_length=length($orf_seq); 
							if ($orf_length>6)
							{
								#print out the current open reading frame to the file that will contain all the orf's for each
								#sequence found - $index will indicate the position of the orf in the sequence 
								
								#don't need if only keeping longest orf - saves memory 
								#$orf{$orf_seq} .= "$name\_$direction\_fr$k\_$index, ";
							}
							$index+=$orf_length+1;
							if ($longest_orf_length<$orf_length)
							{#store the info for the longest orf, will print only longest orf's to a separate file
								$longest_orf_direction=$direction;
								$longest_orf_frame=$k;
								$longest_orf_length=$orf_length;
								$longest_orf_seq=$orf_seq;
							}
						}
					}
				}
				if ($longest_orf_length>6) { $longest_orf{$longest_orf_seq} .= "$name\_$longest_orf_direction\_fr$longest_orf_frame, "; }
			}
			$name=$name_;
			$description=$description_;
			$sequence="";
		}
		else
		{
			$sequence .= "\U$line";
			if ($ftype eq 'FASTQ')
			{#skip next 2 lines
				$line=<IN>;
				$line=<IN>;
			}
			
		}
		
	}	
	if ($name=~/\w/ and $sequence=~/\w/)
	{#do the same as above for the last sequence in the file
				
		my $size = length($sequence);
		my $longest_orf_direction="";
		my $longest_orf_frame="";
		my $longest_orf_length=0;
		my $longest_orf_seq="";
		foreach my $direction ("fwd","rev")
		{
			my $seq="";
			if ($direction=~/^fwd$/) { $seq=$sequence; } 
			else 
			{ 
				$seq = reverse $sequence;
				$seq=~tr/ATCG/TAGC/;
			}
			for (my $k=0;$k<3;$k++)
			{
				my $protein="";
				for(my $n=$k;$n<$size;$n=$n+3)
				{
					my $triplet = substr($seq, $n, 3);
					if ($mapping{$triplet} and $mapping{$triplet}=~/[\w\*]/) { $protein.=$mapping{$triplet}; } else { $protein.="X"; }
				}
				$protein=~s/X+$//;
				$protein=~s/^X+//;
				
				my $temp="$protein*";
				my $index=0;
				while ($temp=~s/^([^\*]*)\*//)   
				{
					my $orf_seq=$1;
					$orf_seq=~s/X+$//;
					$orf_seq=~s/^X+//;
					my $orf_length=length($orf_seq); 
					if ($orf_length>6)
					{
						$orf{$orf_seq}.="$name\_$direction\_fr$k\_$index, ";
					}
					$index+=$orf_length+1;
					if ($longest_orf_length<$orf_length)
					{
						$longest_orf_direction=$direction;
						$longest_orf_frame=$k;
						$longest_orf_length=$orf_length;
						$longest_orf_seq=$orf_seq;
					}
				}
			}
		}
		if ($longest_orf_length>6) { $longest_orf{$longest_orf_seq} .= "$name\_$longest_orf_direction\_fr$longest_orf_frame, "; }
	}
	close(IN);		
}
print "Translation finished: $count_fasta files processed.\n";
if ($count_fasta == 0)
{
	print "Warning: No fasta files processed ($!)!\n";
	exit(3);
}
print "Opening $out_dir/longest_nr.fasta for writing.\n";
if (open (OUT,">$out_dir/longest_nr.fasta"))
{#this file will contain the set of the unique, longest length reading frames found with a count in their description of how many dna sequences 
 #resulted in this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
	print "Success opening $out_dir/longest_nr.fasta for writing.\n";
	foreach my $seq (sort keys %longest_orf)
	{
		my $temp=$longest_orf{$seq};
		if ($temp =~ s/^([^\,]+)\, //) #get characters up to next ', '
		{
			my $name=$1;
			my $count=0;
			
			my $found_name = '';
			while($temp =~ s/^([^\,]+)\, //) #count how many dna sequences result in current $seq and print result
			{
				$count++;
				my $next_name = $1;
				
				#my $next_name_ = $next_name;
				#$next_name_ =~ s/_fwd_fr\d$//;
				#$next_name_ =~ s/_rev_fr\d$//;
				#if (defined $ctla4_hits_5094_7_22{$next_name_})
				#{
				#	$found_name = $next_name;
				#}
			}
			
			#my $name_ = $name;
			#$name_ =~ s/_fwd_fr\d$//;
			#$name_ =~ s/_rev_fr\d$//;
			#if ($found_name || defined $ctla4_hits_5094_7_22{$name_})
			#{
			#	if ($found_name)
			#	{
			#		print "$found_name ($name): $count\n";
			#	}
			#	else
			#	{
			#		print "$name: $count\n";
			#	}	
			#}
			
			if ($count==0) { print OUT qq!>$name\n$seq\n!; }
			else
			{
				print OUT qq!>$name + $count other\n$seq\n!;
				#print OUT qq!>$longest_orf{$seq}\n$seq\n!;
			}
			
			#remember the IDs of all the DNA sequences that produce this AA sequence
			#will be saved to JSON file
			#if (defined $redundancy_listing{$name})
			#{
			#	print "Error: $name already found in JSON hash ($redundancy_listing{$name})";
			#}
			#$redundancy_listing{$name} = $longest_orf{$seq};
		}
	}
	close(OUT);
}
else
{
	print "Error creating $out_dir/longest_nr.fasta ($!)\n";
	exit(4);
}

#save JSON
#my $utf8_encoded_json_text = encode_json \%redundancy_listing;
#if (open (OUT,">$out_dir/repeat_ids.json"))
#{
#	print OUT $utf8_encoded_json_text;
#	close(OUT);
#}

print "Exiting program...\n";
exit(0);


