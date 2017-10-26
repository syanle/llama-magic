#!/usr/local/bin/perl

use File::Copy;
use warnings;
use strict;

#input params are mgf_dir, root_dir

#open parameters file
# index     mgf_name    db_id   parent_error    parent_units    frag_error  frag_units

#for reach row, read settings, name title to match  mgf name with errors in name
#create folder structure under root_dir


my $params_file = "input_params.csv";
my $ms_list_file  = "ms_list.txt";

my $base_dir;
if ($ARGV[0] && $ARGV[0]=~/\w/) { $base_dir=$ARGV[0];}  
else { $base_dir = "/ifs/data/proteomics/projects/NCDIR/Llama/pipeline"; }
#{ $base_dir = "C:/Code/NCDIR/Llama"; } 

my $pipeline_dir = "pipeline";
my $results_dir = "results";

#print "$base_dir/$pipeline_dir/$params_file\n";
open(IN, "$base_dir/$pipeline_dir/$params_file") || die("Failed to open file: $!");
my @lines = <IN>;
close(IN);
for(my $i = 0; $i <= $#lines; $i++)
{
    if ($lines[$i] =~ /(\d+),([^,]+),(\d+),([\d\.]+),(\w+),([\d\.]+),(\w+)/)
    {
        my $i = $1;
        my $mgf_name = $2;
        my $db_id = $3;
        my $parent_err = $4;
        my $parent_err_units = $5;
        my $frag_err = $6;
        my $frag_err_units = $7;
        
        #read list of db searches, get next search id for the db id selected
        open(IN, "$base_dir/$results_dir/$ms_list_file") || die("Failed to open file: $!");
        my @lines = <IN>;
        close(IN);
        
        my $new_ms_id = 1;
        for(my $i = 0; $i <= $#lines; $i++)
        {
            if ($lines[$i] =~ /(\d+)\t(\d+)\t(.+)/)
            {
                if ($1 eq $db_id and $2 >= $new_ms_id) { $new_ms_id = $2+1; }
            }
        }
	
        #add 1000 since it's the pipeline - we don't want it to collide with any searches created manually
        if($new_ms_id < 1000) { $new_ms_id = $new_ms_id + 1000; }
            
        #create new dir structure for the search
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id") || die("Failed to make dir: $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/mgf") || die("Failed to make dir: $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/tandem") || die("Failed to make dir: $!");
        mkdir("$base_dir/$results_dir/$db_id/$new_ms_id/tandem/results") || die("Failed to make dir: $!");
            
        #copy mgf file:
        copy("$base_dir/$pipeline_dir/$mgf_name","$base_dir/$results_dir/$db_id/$new_ms_id/mgf/ms_1.mgf");
        
        my $new_ms_name = "$i.$mgf_name-$parent_err-$parent_err_units-$frag_err-$frag_err_units";
        open(OUT, ">>$base_dir/$results_dir/$ms_list_file") || die("Failed to open file: $!");
        print OUT "$db_id\t$new_ms_id\t$new_ms_name\n"; 
        close(OUT);
        
        #read db params to see if we should use primers and tail sequences
        my $USE_PRIMERS = 0; 
        my $ADD_TAIL_SEQUENCES = 0;
        my $NEW_PRIMERS = 0;
        if(open(IN, "<$base_dir/$results_dir/$db_id/db_params.txt"))
        { 
            while (<IN>)
            {
                if (/USE_PRIMERS=(\d)/) { $USE_PRIMERS = $1; }
                if (/ADD_TAIL_SEQUENCES=(\d)/) { $ADD_TAIL_SEQUENCES = $1; }
                if (/NEW_PRIMERS=(\d)/) { $NEW_PRIMERS = $1; }
            }
            close(IN);
        }
        else { die("Failed to open file: $!"); }
        
        #start qsub
        print "Running job: \n";
        print "qsub run_LM_script.sh \"$base_dir/$results_dir/$db_id/$new_ms_id\" \"$base_dir/$results_dir/$db_id\" \"$parent_err\" \"$frag_err\" \"$parent_err_units\" \"$frag_err_units\" \"/$results_dir/$db_id/$new_ms_id\" \"1\" \"$USE_PRIMERS\" \"$ADD_TAIL_SEQUENCES\" \"0\" \"$NEW_PRIMERS\"";
            print "\n";
        system("qsub run_LM_script.sh \"$base_dir/$results_dir/$db_id/$new_ms_id\" \"$base_dir/$results_dir/$db_id\" \"$parent_err\" \"$frag_err\" \"$parent_err_units\" \"$frag_err_units\" \"/$results_dir/$db_id/$new_ms_id\" \"1\" \"$USE_PRIMERS\" \"$ADD_TAIL_SEQUENCES\" \"0\" \"$NEW_PRIMERS\"")
                    
    }
}


