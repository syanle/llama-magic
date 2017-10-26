#!/usr/local/bin/perl

#for every sequence in each fasta file in the given directory (directory is the argument
#passed on the command line), converts the dna sequence to the possible protein sequences, 
#there are 3 reading frames for both forwards reading and backward reading of the dna sequence, 
#therefore 6 possible protein sequences from a dna sequence

#the $dir\_.fasta file contains all possible sequences
#the $dir\_split.fasta file contains the sequences split up based on the stop codons - all orf's
#the $dir\_longest.fasta file contains only the longest orf for each sequence
#the $dir\_split_nr.fasta file contains all unique orfs with the name of the first dna sequence that produced the orf, and a count of 
#  how many others produced the orf
#the $dir\_longest_nr.fasta file contains only the unique, longest orfs, with the name of the first dna sequence (that produced the orf) 
#  and a count of how many others produced the longest orf
#the $dir\_split_trim_nr.fasta file contains the unique orf's that are trimmed at each end (for what?)

use strict;

my $dir="";

if ($ARGV[0]=~/\w/) { $dir="$ARGV[0]";} else { $dir="/Users/sarahkeegan/Desktop/for_Brian/R2/"; }

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
		
if (opendir(dir,"$dir"))
{
	my %longest_orf=();
	my %orf=();
	my %orf_trim=();
	
	my @allfiles=readdir dir;
	closedir dir;
	
	open (OUT__,">$dir\_split.fasta");
	open (OUT_,">$dir\_longest.fasta");
	
	if (open (OUT,">$dir\_.fasta"))
	{
		foreach my $filename (@allfiles)
		{
			if ($filename =~ /.+\.fas?t?a?$/i) #can be *.fa, *.fas, *.fast, *.fasta, case insensitive
			{#for each fasta file in the directory:
				print qq!Processing: $filename\n!;
				if (open (IN,"$dir/$filename"))
				{
					my $name="";
					my $description="";
					my $sequence="";
					my $line="";
					while ($line=<IN>)
					{
						chomp($line);
						if ($line=~/^>(\S+)\s?(.*)$/)
						{#if the current line is the description line - begins with '>', and then 
						 #one or more non-whitespace, followed by 0 or more whitespace char's, then anything until the end
							my $name_=$1; 
							my $description_=$2;
							$name_ =~ s/[\:\-]/_/g;
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
											if ($mapping{$triplet} =~ /[\w\*]/) { $protein .= $mapping{$triplet}; } # '*' is stop codon
											else { $protein.="X"; } # X is unknown, doesn't code for anything must be error in sequence
										}
										#remove X's at beginning and end
										$protein =~ s/X+$//;
										$protein =~ s/^X+//;
										
										#print out the entire AA sequence for the current direction and reading frame to the output file
										#(directory-name.fasta)
										print OUT qq!>$name\_$direction\_fr$k $direction frame $k\n$protein\n!;
										
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
												print OUT__ qq!>$name\_$direction\_fr$k\_$index $direction frame $k $index\n$orf_seq\n!;
												$orf{$orf_seq} .= "$name\_$direction\_fr$k\_$index, ";
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
								if ($longest_orf_length>6)
								{
									#(if the longest orf is > 6), print the longest orf to a separate file which contains only these for each dna sequence
									print OUT_ qq!>$name\_$longest_orf_direction\_fr$longest_orf_frame $longest_orf_direction frame $longest_orf_frame\n$longest_orf_seq\n!;
									$longest_orf{$longest_orf_seq} .= "$name\_$longest_orf_direction\_fr$longest_orf_frame, ";
								}
							}
							$name=$name_;
							$description=$description_;
							$sequence="";
						}
						else
						{
							$sequence .= "\U$line";
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
											if ($mapping{$triplet}=~/[\w\*]/) { $protein.=$mapping{$triplet}; } else { $protein.="X"; }
										}
										$protein=~s/X+$//;
										$protein=~s/^X+//;
										print OUT qq!>$name\_$direction\_fr$k $direction frame $k\n$protein\n!;
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
												print OUT__ qq!>$name\_$direction\_fr$k\_$index $direction frame $k $index\n$orf_seq\n!;
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
								if ($longest_orf_length>6)
								{
									print OUT_ qq!>$name\_$longest_orf_direction\_fr$longest_orf_frame $longest_orf_direction frame $longest_orf_frame\n$longest_orf_seq\n!;
									$longest_orf{$longest_orf_seq} .= "$name\_$longest_orf_direction\_fr$longest_orf_frame, ";
								}
					}
					close(IN);
				}			
			}
		}
		close(OUT);
	}
	close(OUT_);
	close(OUT__);

	if (open (OUT,">$dir\_longest_nr.fasta"))
	{#this file will contain the set of the unique, longest length reading frames found with a count in their description of how many dna sequences 
	 #resulted in this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
		foreach my $seq (sort keys %longest_orf)
		{
			my $temp=$longest_orf{$seq};
			if ($temp =~ s/^([^\,]+)\, //) #get characters up to next ', '
			{
				my $name=$1;
				my $count=0;
				while($temp =~ s/^([^\,]+)\, //) { $count++; } #count how many dna sequences result in current $seq and print result
				if ($count==0)
				{#only 1 name in the list
					print OUT qq!>$name $name\n$seq\n!;
				}
				else
				{
					print OUT qq!>$name + $count other\n$seq\n!;
				}
			}
			else { print "Error\n"; }
		}
		close(OUT);
	}	
	if (open (OUT,">$dir\_split_nr.fasta"))
	{#this file will contain the set of all unique reading frames found with a count in their description of how many dna sequences 
	 #produced this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
		foreach my $seq (sort keys %orf)
		{
			my $ok=1;
			my $seq_=$seq;
			
			#cut off beginning up to K or R and the end back up to K or R (?)
			if ($seq_ !~ s/^[^KR]*[KR]//) { $ok=0; }
			if ($seq_ !~ s/([KR])[^KR]*$/$1/) { $ok=0; }
			
			my $orf_length=length($seq_); 
			if($orf_length>6 and $ok==1) { $orf_trim{$seq_} .= $orf{$seq}; }
			
			my $temp=$orf{$seq};
			if ($temp =~ s/^([^\,]+)\, //) #get characters up to next ', '
			{
				my $name=$1;
				my $count=0;
				while($temp =~ s/^([^\,]+)\, //) { $count++; } #count how many dna sequences result in current $seq and print result
				if ($count==0)
				{#only 1 name in the list
					print OUT qq!>$name $name\n$seq\n!;
				}
				else
				{
					print OUT qq!>$name + $count other\n$seq\n!;
				}
			}
			else { print "Error\n"; }
		}
		close(OUT);
	}
	
	if (open (OUT,">$dir\_split_trim_nr.fasta"))
	{#this file will contain the set of all unique (and trimmed) reading frames found with a count in their description of how many dna sequences 
	 #produced this AA sequence - only the name of the first dna sequence is outputted, following a number for the count
	 #the sequence is trimmed by cutting off the beginning and end up to the first encountered K or R
		foreach my $seq (sort keys %orf_trim)
		{
			my $temp=$orf_trim{$seq};
			if ($temp =~ s/^([^\,]+)\, //)
			{
				my $name=$1;
				my $count=0;
				while($temp =~ s/^([^\,]+)\, //) { $count++; }
				if ($count==0)
				{
					print OUT qq!>$name $name\n$seq\n!;
				}
				else
				{
					print OUT qq!>$name + $count other\n$seq\n!;
				}
			}
			else { print "Error\n"; }

		}
		close(OUT);
	}
}
