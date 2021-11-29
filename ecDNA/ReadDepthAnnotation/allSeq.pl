#!/gsc/bin/perl

use warnings;
use strict;
use IO::File;
use File::Basename;

if (@ARGV < 2){   
    die "need more arguments
  arg 0 = fasta file
  arg 1 = read length
";
}



my $seq = "";
my $inFh = IO::File->new( $ARGV[0] );
while( my $line = $inFh->getline )
{
    chomp($line);
    $seq = $seq . $line unless $line =~ /\>/;    
}


my $readSize = $ARGV[1];
my $stop = length($seq)-$readSize;

my $qual = "";
for(my $j=0;$j<$readSize;$j++)
{
    $qual = $qual . "~";
}


for(my $i=0;$i<$stop;$i++)
{
    my $read = substr($seq,$i,$readSize);
    #don't print lines that are all Ns
    unless ($read =~ /^[N|n]+$/)
    {
	print "@" . basename($ARGV[0]) . ":" . ($i+1) . "-" . ($i+$readSize) . "\n";
	print $read . "\n";
	print "+\n";
	print $qual . "\n";
    }
}
