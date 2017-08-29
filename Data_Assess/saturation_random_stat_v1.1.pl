#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($Config,$fOut,$Q);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"dataconfig:s"=>\$Config,
				"Q:s"=>\$Q,
				) or &USAGE;
&USAGE unless ($Config and $fOut);
$Q||=33;
my %seq_hash;
#my $cfg=&ABSOLUTE_DIR($Config);
my $cfg=$Config;
#my $outdir=&ABSOLUTE_DIR($fOut);
my $outdir=$fOut;
unless (-d "$outdir/work_sh")
{
	`mkdir $outdir/work_sh`;
}
my $qc_stat="$Bin/fastq_qc_stat/v1.3/fastq_qc_stat";
my $base_quality="$Bin/bin/GCqual_Solexa_check_svg_final.pl";
my $GC_dis="$Bin/bin/GCcont_Solexa_check_svg_final_v1.pl";
open (IN,$cfg) or die $!;

my $sample;
while (<IN>) 
{
	chomp;
	next if ($_=~/\#/);
	next if ($_=~/^\s*$/);
	my @data=split /\t/,$_;
	if ($data[0]=~/Sample/i)
	{
		$sample=$data[1];
	}
	if ($data[0]=~/^fq1/)
	{
		$seq_hash{$sample}{1}=$data[1];
	}
	if ($data[0]=~/^fq2/)
	{
		$seq_hash{$sample}{2}=$data[1];
	}
}
open (QCSTAT,">$outdir/work_sh/qc_stat.sh") or die $!;
open (BQUALITY,">$outdir/work_sh/base_quality.sh") or die $!;
open (GCDIS,">$outdir/work_sh/gc_dis.sh") or die $!;
foreach my $sam (sort keys %seq_hash)
{
	unless (-d "$outdir/$sam")
	{
		`mkdir $outdir/$sam`;
	}
	print QCSTAT "cd $outdir/$sam && $qc_stat -a $seq_hash{$sam}{1} -b $seq_hash{$sam}{2} -f $sam -q 45 -Q $Q\n";
	#print BQUALITY "perl $base_quality -qu $outdir/$sam/$sam.quality -od $outdir/$sam\n"; 
	print BQUALITY "Rscript  $Bin/bin/quality_bar_v1.R --infile $outdir/$sam/$sam.quality --outfile $outdir/$sam/$sam.quality --col \"#377EB8\" -s $sam\n";#huangls 2015-7-31
	#print GCDIS "perl $GC_dis -gc $outdir/$sam/$sam.acgtn -od $outdir/$sam\n"; 
	print GCDIS "Rscript $Bin/bin/gc_plot.r --gc $outdir/$sam/$sam.acgtn --output $outdir/$sam --name $sam\n";#huangls 2015-7-31
}
close(QCSTAT);
close(BQUALITY);
close(GCDIS);
&qsub("$outdir/work_sh/qc_stat.sh");
&qsub("$outdir/work_sh/base_quality.sh");
&qsub("$outdir/work_sh/gc_dis.sh");
my @Sample=glob("$outdir/*/");
my @Paraf;
foreach my $sample (@Sample) 
{
	my @STAT=glob("$sample/*.stat");
	foreach my $stat (@STAT) 
	{
		push @Paraf,$stat unless $stat=~/.+cycleQ\.stat/;
	}
}
open (OUT,">$outdir/AllSample_GC_Q.stat");
print OUT "SampleID\tReadSum\tBaseSum\tGC(%)\tN(%)\tQ20(%)\tCycleQ20(%)\tQ30(%)\n";
foreach my $paraf (@Paraf) 
{
	my $now;
	my $file=basename($paraf);
	$file=~s/\.stat$//;
	open(IN,"$paraf")||die"can't open $paraf\n";
	<IN>;
	while (<IN>)
	{
		chomp;
		my @A=split/\s+/,$_,2;
		$now="$file\t$A[1]\n";
	}
	close(IN);
	print OUT $now;
}
close(OUT);
#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
#######################################################################################
sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
#############################################################
#&run_or_die($cmd);
sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}
######################################################################
## qsub
sub qsub()
{
	my ($shfile, $queue, $ass_maxproc) = @_ ;
	$queue ||= 'general.q' ;
	$ass_maxproc ||= 20 ;
	my $cmd = "sh $shfile" ;
	&run_or_die($cmd);
}
#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	Wujh <wujh\@biomarker.com.cn> 
Program Date:   2014.9.1
Description:	this program is used to ......
Usage:
  Options:
  -dataconfig <file>  input config file
	
  -Q <file>  quality shift for different plotform [optional] [default 33]
  
  -o <file>  outdir 
 
  -h         Help

USAGE
	print $usage;
	exit;
}
