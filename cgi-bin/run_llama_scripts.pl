#!/usr/local/ActivePerl-5.24/bin/perl

#runs the scripts sent in by command line arguments for (CGI web script) llama_magic.pl
#options are:
#db_scripts (runs both dnatoprot_longest_nr.pl and digest_fasta.pl)
#search_and_map_scripts (runs xtandem (command line/xml interface), parse_xtandem.pl, filter_tandem_results.pl, and map_peptides_to_proteins.pl)

# "search_and_map_scripts" "C:/Code/NCDIR/Llama/results/64/31" "C:/Code/NCDIR/Llama/results/64" "10" "0.4" "ppm" "Daltons" "/results/64/31" 1 1 1 0

use warnings;
use strict;

my $include_extra_files = 0;
my $expect_threshold = '0.1';
my $function = $ARGV[0];

if ($include_extra_files)
{
    use File::Copy;
}

if ($function eq 'db_scripts')
{
    my $results_dir = $ARGV[1];
    my $use_primers = $ARGV[2];
    my $add_tails = $ARGV[3];
    
    my $error_fail = 0;
    
    #create status txt file
    open(OUT, ">$results_dir/protein/status.txt") or die "Failed to create status file: $results_dir/protein/status.txt";
    
    #copy vhh_nuc.fasta to dna folder
    if ($include_extra_files && !$use_primers)
    {
	if(!copy("../vhh_nuc.fasta", "$results_dir/dna/vhh_nuc.fasta"))
	{
	    print OUT "ERROR: Copying of '../vhh_nuc.fasta' to '$results_dir/dna/vhh_nuc.fasta' failed.\n";
	}
    }
    
    my $cmd_out;
    if ($use_primers)
    {
	unless(-d "$results_dir/dna/primers_removed")
	{
	    if(!mkdir "$results_dir/dna/primers_removed")
	    {
		print OUT "ERROR: Cannot create directory '$results_dir/dna/primers_removed': $!\n";
		$error_fail = 1;
	    }
	}
	
	print OUT qq!Calling: "../remove_primers.pl" "$results_dir/dna" "$results_dir/dna/primers_removed"\n!;
	$cmd_out = `"perl" "../remove_primers.pl" "$results_dir/dna" "$results_dir/dna/primers_removed" 2>&1`;
	if ( $? == -1 )
	{
	    print OUT "ERROR: command failed (remove_primers.pl): $!\n";
	    $error_fail = 1;
	}
	elsif($? >> 8 != 0) # exit value not 0, indicates error...
	{
	    printf OUT "ERROR: command (remove_primers.pl) exited with value %d\n", $? >> 8;
	    print OUT "$cmd_out\n";
	    $error_fail = 1;
	}
    }
    
    if (!$error_fail)
    {
	my $input_dir;
	if ($use_primers)
	{
	    print OUT "Success: $cmd_out\n";
	    $input_dir = "$results_dir/dna/primers_removed";
	}
	else { $input_dir = "$results_dir/dna"; }
	
	#run dnatoprot
	print OUT qq!Calling: "../dnatoprot_longest_nr_with_primers.pl" "$input_dir" "$results_dir/protein" "$use_primers"\n!;
	$cmd_out = `"perl" "../dnatoprot_longest_nr_with_primers.pl" "$input_dir" "$results_dir/protein" "$use_primers" 2>&1`;
	if ( $? == -1 )
	{
	    print OUT "ERROR: command failed (dnatoprot_longest_nr.pl): $!\n";
	    $error_fail = 1;
	}
	elsif($? >> 8 != 0) # exit value not 0, indicates error...
	{
	    printf OUT "ERROR: command (dnatoprot_longest_nr.pl) exited with value %d\n", $? >> 8;
	    print OUT "$cmd_out\n";
	    $error_fail = 1;
	}
    }
    
    if (!$error_fail)
    {
        print OUT "Success: $cmd_out\n";
        
        #Success! add vh_prot.fas sequences to the longest_nr file:
	if ($include_extra_files && !$use_primers)
	{
	    if(!open(PROTEINS_IN, "../vh_prot.fasta"))
	    {
		print OUT "ERROR: Opening of '../vh_prot.fasta' for reading failed: $!\n";
	    }
	    else
	    {
		if(!open(PROTEINS_IN_, "../crap.fasta"))
		{
		    print OUT "ERROR: Opening of '../crap.fasta' for reading failed: $!\n";
		}
		else
		{
		    if(!open(PROTEINS_OUT, ">>$results_dir/protein/longest_nr.fasta"))
		    {
			print OUT "ERROR: Opening of '$results_dir/protein/longest_nr.fasta' for appending failed: $!\n";
		    }
		    else
		    {
			print PROTEINS_OUT "\n"; #just in case
			while(<PROTEINS_IN>)
			{
			    print PROTEINS_OUT $_;
			}
			print PROTEINS_OUT "\n"; #just in case
			while(<PROTEINS_IN_>)
			{
			    print PROTEINS_OUT $_;
			}
			
			close(PROTEINS_IN);
			close(PROTEINS_OUT);
			print OUT "Success: Added vh_prot.fasta and crap.fasta to longest_nr.fasta.\n";
		    }
		}
		
	    }
	}
        
	#success! run digest_fasta
	print OUT qq!Calling: "../digest_fasta.pl" "$results_dir/protein" "1" "$add_tails"\n!;
	$cmd_out = `"perl" "../digest_fasta.pl" "$results_dir/protein" "1" "$add_tails" 2>&1`;
	if ( $? == -1 )
	{
	    print OUT "ERROR: command failed (digest_fasta.pl): $!\n";
	    $error_fail = 1;
	}
	elsif($? >> 8 != 0) # exit value not 0, indicates error...
	{
	    printf OUT "ERROR: command (digest_fasta.pl) exited with value %d", $? >> 8;
	    print OUT "$cmd_out\n";
	    $error_fail = 1;
	}
	else
	{
	    print OUT "Success: $cmd_out\n";
	    
	    #success!  rename all_predigested.fasta to longest_nr_predigested.fasta
	    if (!rename "$results_dir/protein/all_predigested.fasta", "$results_dir/protein/longest_nr_predigested.fasta")
	    {
		print OUT "ERROR: Could not rename '$results_dir/protein/all_predigested.fasta' to '$results_dir/protein/longest_nr_predigested.fasta'.\n";
		$error_fail = 1;
	    }
	}
    }
    #write to status file and close
    print OUT "DONE\n";
    close(OUT); 
}
elsif($function eq 'search_and_map_scripts')
{
    my $results_dir = $ARGV[1];
    my $taxon_dir = $ARGV[2];
    my $parent_err = $ARGV[3];
    my $fragment_err = $ARGV[4];
    my $parent_err_units = $ARGV[5];
    my $frag_err_units = $ARGV[6];
    my $cgi_results_dir = $ARGV[7];
    my $show_score = $ARGV[8];
    my $use_primers = $ARGV[9];
    my $add_tails = $ARGV[10];
    my $use_constants = $ARGV[11];
    my $new_primers = $ARGV[12];
    
    #create status txt file
    open(OUT, ">$results_dir/status.txt") or die "Failed to create status file: $results_dir/status.txt";
    
    #create xtandem input.xml
    if(open(XML_OUT, ">$results_dir/tandem/input.xml"))
    {
        print XML_OUT <<XMLTEXT;
<?xml version="1.0"?>
<bioml>
	<note type="input" label="spectrum, fragment monoisotopic mass error">$fragment_err</note>
	<note type="input" label="spectrum, parent monoisotopic mass error plus">$parent_err</note>
	<note type="input" label="spectrum, parent monoisotopic mass error minus">$parent_err</note>
	<note type="input" label="spectrum, fragment monoisotopic mass error units">$frag_err_units</note>
	<note type="input" label="spectrum, parent monoisotopic mass error units">$parent_err_units</note>

	<note type="input" label="list path, default parameters">../default_input.xml</note>
	<note type="input" label="list path, taxonomy information">$taxon_dir/taxonomy.xml</note>
XMLTEXT
	if($use_constants)
	{
	    print XML_OUT qq(\n\t<note type="input" label="protein, taxon">llama_constant</note>);
	}
	else
	{
	    print XML_OUT qq(\n\t<note type="input" label="protein, taxon">llama</note>);
	}
	
    #get list of mgf files and add a line in input.xml for each one
    my @file_names = <$results_dir/mgf/*.mgf>; 
    foreach my $cur_file (@file_names)
	{
        print XML_OUT qq(\n\t<note type="input" label="spectrum, path">$cur_file</note>);
    }
        
    print XML_OUT <<XMLTEXT;
	<note type="input" label="output, path">$results_dir/tandem/results/output.xml</note>
	<note type="input" label="output, results">valid</note>
</bioml>
XMLTEXT

        close(XML_OUT);
    
        #run xtandem with uploaded mgf files
        print OUT qq!Calling: "../tandem" "$results_dir/tandem/input.xml"\n!;
        my $cmd_out = `"../tandem" "$results_dir/tandem/input.xml" 2>&1`;
        if ( $? == -1 )
        {
            print OUT "ERROR: command failed (tandem.exe): $!\n";
        }
        elsif($? >> 8 != 0) # exit value not 0, indicates error...
        {
            printf OUT "ERROR: command (tandem.exe) exited with value %d\n", $? >> 8;
            print OUT "$cmd_out\n";
        }
        else
        {
            #parse xtandem to tab separated txt file, extract relevant columns
            print OUT qq!Calling: "../parse_xtandem_llama.pl" "$results_dir/tandem/results" "$expect_threshold"\n!;
            my $cmd_out = `"perl" "../parse_xtandem_llama.pl" "$results_dir/tandem/results" "$expect_threshold" 2>&1`;
            if ( $? == -1 )
            {
                print OUT "ERROR: command failed (parse_xtandem_llama.pl): $!\n";
            }
            elsif($? >> 8 != 0) # exit value not 0, indicates error...
            {
                printf OUT "ERROR: command (parse_xtandem_llama.pl) exited with value %d\n", $? >> 8;
                print OUT "$cmd_out\n";
            }
            else
            {
                #get name of txt file created above:
                #my @files = <$results_dir/tandem/results/*.txt>; 
                #my $pep_file = $files[0]; #should be exactly 1 txt file in directory
                
                my $pep_file = "output.xml.peptide_list.$expect_threshold.txt";
                
                #run map_peptides_to_proteins - output is candidate list html file
                
		#need to change this, need to add input of whether using primers / make changes in map_peptides_to_proteins.pl...
                print OUT qq!Calling: "../map_peptides_to_proteins.pl" "$results_dir/tandem/results/$pep_file" "$taxon_dir/protein/longest_nr.fasta" "$cgi_results_dir/tandem/results/output.xml" "$taxon_dir/protein/protein_peptides.fasta" "$show_score" "$use_primers" "$add_tails" "$new_primers" \n!;
                my $cmd_out = `"perl" "../map_peptides_to_proteins.pl" "$results_dir/tandem/results/$pep_file" "$taxon_dir/protein/longest_nr.fasta" "$cgi_results_dir/tandem/results/output.xml" "$taxon_dir/protein/protein_peptides.fasta" "$show_score" "$use_primers" "$add_tails" "$new_primers" 2>&1`;
                if ( $? == -1 )
                {
                    print OUT "ERROR: command failed (map_peptides_to_proteins.pl): $!\n";
                }
                elsif($? >> 8 != 0) # exit value not 0, indicates error...
                {
                    printf OUT "ERROR: command (map_peptides_to_proteins.pl) exited with value %d\n", $? >> 8;
                    print OUT "$cmd_out\n";
                }   
            }
            
        }
    }
    else { print OUT "ERROR: Failed to create xtandem input file: $results_dir/tandem/input.xml"; }     
    
    #write to status file and close
    print OUT "DONE\n";
    close(OUT); 
}


