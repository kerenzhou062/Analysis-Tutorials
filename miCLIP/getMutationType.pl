#!/usr/bin/perl -w
# modified from CTK
# add c2t mutation
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Carp;
use Bed;
use Sequence;

my $prog=basename($0);
my $cmdDir=dirname($0);

my $selectType = "";
my $t2c = 0;
my $c2t = 0;
my $nfrom;
my $nto;
my $summaryFile = "";
my $verbose = 0;

GetOptions ('t|mutationType=s'=>\$selectType,
            'nfrom:s'=>\$nfrom,
            'nto:s'=>\$nto,
			'summary:s'=>\$summaryFile,
			'v'=>\$verbose);


if (@ARGV != 2)
{
	print "Get specific types of mutations\n";
	print "usage: $prog [options] <tag.mutation.txt> <tag.mutation.type.bed>\n";
	print " -t [string]        : the type of mutations (del|ins|sub|t2c|c2t)\n";
    print " --nfrom [nucleotide]   : mutated nucleotide (a|t|c|g)\n";
    print " --nto [nucleotide]   : converted nucleotide (a|t|c|g)\n";
	print " --summary  [string]: print summary statistics to file\n";
	print " -v                 : verbose\n";
	exit (1);
}

my ($mutationFile, $outBedFile) = @ARGV;
my %summaryHash;

if ($selectType ne 'del' && $selectType ne 'ins' && $selectType ne 'sub' && $selectType ne 't2c' && $selectType ne 'c2t')
{
	Carp::croak "Wrong type of mutation \n";
}


if ($selectType eq 't2c')
{
	$t2c = 1;
	$selectType = 'sub';
}

if ($selectType eq 'c2t')
{
    $c2t = 1;
    $selectType = 'sub';
}

if (defined($nfrom) && defined($nto))
{
    if ($nfrom =~ m[ATCG] && $nto =~ m[ATCG]) {
        $nfrom = uc ($nfrom);
        $nto = uc ($nto);
    }
}
my ($fin, $fout);

open ($fin, "<$mutationFile") || Carp::croak "cannot open $mutationFile to read\n";
open ($fout, ">$outBedFile") || Carp::croak "cannot open $outBedFile to write\n";

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	print "$iter ...\n" if $verbose && $iter % 100000 == 0;

	my @cols = split (/\s+/, $line);
	
	my $strand = $cols[5];
	my ($from, $type, $to) = @cols [7..9];


	if ($type eq '-')
	{
		$summaryHash{'del'}++;
		print $fout join ("\t", @cols[0..5]), "\n" if $selectType eq 'del';
	}
	elsif ($type eq '+')
	{
		$summaryHash{'ins'}++;
		print $fout join ("\t", @cols[0..5]), "\n" if $selectType eq 'ins';
	}
	elsif ($type eq '>')
	{

		if ($strand eq '-')
		{
			$from = uc (revcom ($from));
			$to = uc (revcom ($to));
		}
		
		$summaryHash{'sub'}{"$from-$to"}++;
		$summaryHash{'sub'}{'total'}++;

		if ($selectType eq 'sub')
		{
			if ($t2c == 1)
			{
				next unless $from eq 'T' && $to eq 'C';
			}elsif ($c2t == 1)
            {
                next unless $from eq 'C' && $to eq 'T';
            }elsif (defined($nfrom) && defined($nto))
            {
                if ($nfrom =~ m[ATCG] && $nto =~ m[ATCG]) {
                    next unless $from eq $nfrom && $to eq $nto;
                }
            }
			print $fout join ("\t", @cols[0..5]), "\n";
		}
	}
}

close ($fin);
close ($fout);

if ($summaryFile ne '')
{
	open ($fout, ">$summaryFile") || Carp::croak "cannot open file $summaryFile to write\n";

	print $fout "Deletion: ", exists $summaryHash{'del'} ? $summaryHash{'del'} : 0, "\n";
	print $fout "Insertion: ", exists $summaryHash{'ins'} ? $summaryHash{'ins'} : 0, "\n";

	if (exists $summaryHash{'sub'})
	{
		print $fout "Substitution: ", $summaryHash{'sub'}{'total'}, "\n";
		foreach my $t (sort keys %{$summaryHash{'sub'}})
		{
			next if $t eq 'total';
			print $fout "$t: ", $summaryHash{'sub'}{$t}, "\n";
		}
	}

	close ($fout);
}


