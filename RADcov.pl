#!/usr/bin/perl

# Script takes an input directory of bam files created by Stacks (from RADseq data)
# For each bam file runs the SAMtools depth function, then runs the RADcoverage shell script
# Samtools must be loaded before running

### Load Perl modules ###

use strict;
use warnings;
use Cwd;
use Parallel::ForkManager;
use Getopt::Std;

### Get Options ###

our $opt_d = getcwd;
# The directory with the bam data 
# If there is no "raw data" directory specified the script will default to the current working directory

our $opt_p = 15;
# The number of simultaneous processes to use in the steps managed by parallel fork manager (default is 15)

getopts("d:p:");

### Set variables and print settings ###

my $bamdir = $opt_d;
my $Depth_dir = $bamdir . "/Depth";
my $Cov_dir = $bamdir . "/Coverage";

print "~~~ Settings as specified: ~~~ \n";
print "bam directory: " . $bamdir . "\n";
print "depth directory: " . $Depth_dir . "\n";
print "coverage directory: " . $Cov_dir . "\n";
print "number of threads: " . $opt_p . "\n";
print "\n";

### Make Directories ###

mkdir $Depth_dir;
mkdir $Cov_dir;

### Open bam directory and list files ###

opendir DIR,$bamdir;
my @bamfiles = readdir DIR;
closedir DIR; 

### Initialise variables before loop ###

my $Min;
my $Depth;
my $Mout;

my $DepthFM = Parallel::ForkManager->new($opt_p);

foreach my $file (@bamfiles) {

	$DepthFM->start and next;
  if ($file =~ /bam/) {
  
    print "working on: " . $file . "\n";

		$Min = $bamdir . "/". $file;
    $Mout = $Cov_dir . "/" . $file . ".cov";
    $Depth = $Depth_dir . "/" . $file . ".depth";
    
		system("samtools depth $Min > $Depth");
		system("DATE=\$(date) \n echo \"Done running samtools for $file It is now \$DATE\" ");
    system("sh /home/francejm/JaFrance/RAD_Scripts/RADcoverage.sh $Depth > $Mout");
    
    print "done working on: " . $file . "\n"; 
     
		}
  $DepthFM->finish;
	}
$DepthFM->wait_all_children();


	

	


	







