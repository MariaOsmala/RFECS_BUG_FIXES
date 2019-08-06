#!usr/bin/perl

use strict;

MAIN:{
    my($chrom_file,$infile,$outfile)=@ARGV;
    my @chr_set=readpipe("cut -f1,1 $chrom_file");
    my @all_lens=readpipe("cut -f1,1 chr_file_lens");
    chomp @all_lens;
    chomp @chr_set;
    open(INF2,$infile);
    my $start=0;
    my $line_num=0;
    foreach my $chr(@chr_set)
    {
        my $outfile2=$outfile."/all_mods_".$chr;
        my $len=$all_lens[$line_num];
        print "chromosome file length:".$len."\n";
        open(my $out,">",$outfile2);
        my $i;
        for($i=$start;$i<$start+$len;$i++)
        {
            my $line=<INF2>;
            print $out $line;
        }
        $start=$i;
        $line_num=$line_num+1;
    }
}
