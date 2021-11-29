#!/gsc/bin/perl

use warnings;
use strict;
use IO::File;
use POSIX;

if (@ARGV < 2){
    die "need more arguments
  arg 0 = sam output from bwa
  arg 1 = list of chrs and lengths
  arg 2 = window size
  arg 3 = chromosome we're working on
  arg 4 = read length
  arg 5 = ouput directory
";
}

my $windSize = $ARGV[2];

#first hash reference lengths (entrypoints)
my %chrLengths = ();
my $inFh = IO::File->new( $ARGV[1] );
while( my $line = $inFh->getline )
{
    chomp($line);
    my @fields = split("\t",$line);
    $chrLengths{$fields[0]} = $fields[1];
}
$inFh -> close;


my @posArr = ();
my @gcArr = ();
my @nArr = ();
my $chr = $ARGV[3];

$inFh = IO::File->new( $ARGV[0] );
while( my $line = $inFh->getline )
{
    chomp ($line);
    #skip header lines
    unless ($line =~ /^@/)
    {
        #assumes the file is named with suffix .fa or .fasta
	if ($line =~ /^([^\.]+)\.fa(sta)?:(\d+)-\d+/)
	{
	    my $start = $3;
	    my $thisChr = $1;
	    my @fields = split("\t",$line);
	    unless ($fields[2] eq "*"){
		
#		print STDERR "$chr == $fields[2]\n";
#		print STDERR "$start == $fields[3]\n";	       
#		print STDERR "$chr == $thisChr\n";

		if ( ($fields[2] eq $chr) && ($fields[3] == $start) 
		     && $chr eq $thisChr)
		{ 
		    #calculate appropriate window
		    my $key = floor($start/$windSize);
#		    print STDERR $key . "\n";
		    if (!(exists($posArr[$key]))){
			$posArr[$key] = 1;
			$gcArr[$key] = ($fields[9] =~ tr/GCgc//);
			$nArr[$key] = ($fields[9] =~ tr/nN//);
		    } else {
			$posArr[$key]++;
			$gcArr[$key] = $gcArr[$key] + ($fields[9] =~ tr/gcGC//);
			$nArr[$key] = $nArr[$key] +($fields[9] =~ tr/nN//);
		    }
		}	
	    }	
	}
    }
}

my $outFh = open (my $mapfile, ">$ARGV[5]/$ARGV[3].map") || die "Can't open output file.\n";
my $outFh2 = open (my $gcfile, ">$ARGV[5]/$ARGV[3].gc") || die "Can't open output file.\n";

my $len = $chrLengths{$chr};

for(my $i=0; $i<($len/$windSize); $i++){
    if (!(exists($posArr[$i]))){
	print $mapfile "0\n";
	print $gcfile "NA\n";
    } else {
	print $mapfile $posArr[$i]/$windSize . "\n";

	my $basecount =  ($posArr[$i] * $ARGV[4]);
	print $gcfile ($gcArr[$i] / ($basecount - $nArr[$i])) . "\n";
    }
}
