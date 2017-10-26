#!/usr/local/bin/perl
#filters a tab-delimited text file by extracting the columns based on column header, then outputting the columns
#of interest to another tab-delimited text file, (in the order given in the input)

#arg1 is a file or directory (where all files w/ extension xls will be read in as tab delimited text files)
#arg2 is a file name that contains a list of columns of interest in the order that is expected in the output file

use strict;

my $input_fname_dir; my $columns_file;
my @output_col_names;

if ($ARGV[0]=~/\w/) { $input_fname_dir=$ARGV[0];} else { $input_fname_dir="C:\\NCDIR\\Llama\\4761\\to_map\\For Mapping_Flag_61_manually_processed_011814.txt"; }
if ($ARGV[1]=~/\w/) { $columns_file=$ARGV[1];} else { $columns_file="columns.txt"; }

#input the column names/order from the columns file:
open(IN, "$columns_file") || die "Error: Could not open columns input file, \'$columns_file\'.\n";
my $cols = "";
while(<IN>)
{
	chomp();
	$cols .= " " . $_;
}
@output_col_names = split(' ', $cols);
close(IN);

#foreach my $col (@output_col_names) { print "$col\n"; }

my $files_converted = 0;

#is fasta_name_dir a file or a directory?
if(-f $input_fname_dir && ($input_fname_dir =~ /\.xls$/ || $input_fname_dir =~ /\.txt$/))
{
	filter_tandem_results($input_fname_dir);
	$files_converted++;
}
elsif(-d $input_fname_dir)
{
	opendir(DIR, "$input_fname_dir") || die "Error: Could not open directory $input_fname_dir.\n";
	my @input_files = readdir DIR;
	foreach my $file (@input_files)
	{
		#if($file =~ /\.xls$/i)
		if($file =~ /\.xls$/i || $file =~ /\.txt$/i)
		{
			filter_tandem_results("$input_fname_dir/$file");
			$files_converted++;	
		}
	}
}
else { print OUT "Error: $input_fname_dir is not an (xls/txt) file or directory.\n"; }

print "Done! $files_converted files converted.\n";

sub filter_tandem_results
{
	my $fname = shift;
	my $outfile = $fname;
	if(!($outfile =~ s/\.xls$/\.txt/))
	{#no xls extension, try txt extension
		$outfile =~ s/\.txt$/_\.txt/;
	}
	
	print "Filtering file $fname.\n";
	
	my %col_names;
	
	if(open(IN, $fname))
	{
		#open outfile
		if(!open(OUT, ">$outfile")) { close(IN); print "Error: Could not open $outfile for writing.\n"; return 0; }
		
		#read in the column headers
		my $line = <IN>;
		chomp($line);
		my @line = split(/\t/, $line);
		for(my $i = 0; $i <= $#line; $i++) { $col_names{$line[$i]} = $i; }
		
		#print out desired columns to out file
		my $first = 1;
		for(my $i = 0; $i <= $#output_col_names; $i++) 
		{ 
			if(exists $col_names{$output_col_names[$i]}) 
			{ 
				if($first) { print OUT "$output_col_names[$i]"; $first = 0; } 
				else { print OUT "\t$output_col_names[$i]"; }
			}
			else { print "Warning: Column name \'$output_col_names[$i]\' does not exist in file $fname.\n"; }
		}
		print OUT "\n";
			
		#read in the column data, printing out the headers that we want...
		#fix!  only want headers w/ expect >= 3!
		while(<IN>)
		{
			chomp();
			my @line = split(/\t/);
			$first = 1;
			for(my $i = 0; $i <= $#output_col_names; $i++) 
			{ #for the desired column names, in the desired order, go through and find the values in the line and print them
				 
				if(exists $col_names{$output_col_names[$i]}) 
				{ 
					if($col_names{$output_col_names[$i]} > $#line) 
					{ print "Warning: missing data in file $fname: $col_names{$output_col_names[$i]}, $#line\n     line = $_\n"; }
					else 
					{ 
						if($first) { print OUT "$line[$col_names{$output_col_names[$i]}]"; $first = 0; } 
						else { print OUT "\t$line[$col_names{$output_col_names[$i]}]"; }
					}
				}
			}
			print OUT "\n";
		}
		close(OUT);
		close(IN);
	}
	else { print "Error: Could not open $fname for reading.\n"; return 0; }
}

#################


