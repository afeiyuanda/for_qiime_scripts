#!/usr/bin/perl -w
# Modifier:	Chen Xiang <chenx@biomarker.com.cn>
# Last Modified:	2014/5/29.
use strict;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use FindBin qw($Bin $Script);

my $programe_dir = basename($0);
my $path         = dirname($0);

my $ver    = "1.0";
my $Writer = "Chen Xiang <chenx\@biomarker.com.cn>";
my $Data   = "2014/5/29";
my $BEGIN  = time();

#######################################################################################
my ($quality,$od,$help);
GetOptions(
	'qu=s' => \$quality,
	'od=s' => \$od,
	'h' => \$help,
);
unless ($quality && $od && !$help) {
	&help;
}

mkdir $od if ($od);
$od=&ABSOLUTE_DIR($od);
$quality=&ABSOLUTE_DIR($quality);

my $Space=10;
my $frame_width=400;
my $frame_height=300;
my $X_len=$frame_width-$Space;
my $Y_len=$frame_height-$Space;
my $left_side=20;
my $right_side=20;
my $up_side=20;
my $down_side=20;

my $TitleSize=16;
my $YtextSize=14;
my $YnumSize=12;
my $YspaceLen=3;
my $XtextSize=14;
my $XnumSize=12;
my $XspaceLen=3;

my $IntextSize=10;

my $height_paper=$frame_height+$up_side+$TitleSize+$XtextSize+$XnumSize+$down_side+2*$XspaceLen;
my $width_paper=$frame_width+$left_side+$YtextSize+$YnumSize+$right_side+$YspaceLen;

########################################################################################
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart $programe_dir Time :[$Time_Start]\n\n";
########################################################################################

#### base quality file:
open IN,$quality or die $!;
my %quality;
my @title=split /\s+/,<IN>;
my $cycle_num=0;
while (<IN>) {
	chomp;
	next if (/^$/);
	$cycle_num++;
	my @units=split /\s+/,$_;
	for (my $i=1;$i<@units ;$i++) {
		$quality{$cycle_num}{$title[$i]}=$units[$i];
	}
}
close IN;

#### svg output:
my $quality_file=basename($quality);
open OUT,">$od/$quality_file.svg" or die $!;
print OUT &svg_paper($width_paper,$height_paper);
print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width/2,$up_side+$TitleSize/2,$TitleSize,"black","Quality Distribution along Reads");  ## Title
print OUT &svg_mid_txt($left_side+$YtextSize/2-$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height/2,$YtextSize,"black","Quality Score",3);  ## Y text
print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width/2,$up_side+$TitleSize+$XspaceLen+$frame_height+$XspaceLen+$XnumSize+$XtextSize+$XspaceLen,$XtextSize,"black","Reads Position,bp");   ## X text

## 矩形
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen,"black",0.5);   ## up line
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## left line
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$Y_len/8-1.5,"black",0.5);   ## left line (up)
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$Y_len/8+1.5,$left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## left line (down)

print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## right line
print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+$frame_height,$left_side+$YtextSize+$YnumSize+$XspaceLen+$frame_width,$up_side+$TitleSize+$XspaceLen+$frame_height,"black",0.5);   ## down line

## Y轴刻度
my $Ystep=10;
my $Y_base=$Y_len/$title[-1];
my $YstepNum=$title[-1]/$Ystep;
for (my $i=0;$i<=$YstepNum ;$i++) {
	print OUT &svg_mid_txt($left_side+$YtextSize,$up_side+$TitleSize+$XspaceLen+$YnumSize*0.4+($Y_len-$i*$Ystep*$Y_base),$YnumSize,"black",$i*$Ystep);
	print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen,$up_side+$TitleSize+$XspaceLen+($Y_len-$i*$Ystep*$Y_base),$left_side+$YtextSize+$YnumSize+$YspaceLen+$YspaceLen,$up_side+$TitleSize+$XspaceLen+($Y_len-$i*$Ystep*$Y_base),"black",0.5);
}

## Y轴不定坐标线
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen-3,$up_side+$TitleSize+$XspaceLen+$Y_len/8-3,$left_side+$YtextSize+$YnumSize+$YspaceLen+3,$up_side+$TitleSize+$XspaceLen+$Y_len/8,"black",0.5);
#print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen-3,$up_side+$TitleSize+$XspaceLen+$Y_len/8,$left_side+$YtextSize+$YnumSize+$YspaceLen+3,$up_side+$TitleSize+$XspaceLen+$Y_len/8+3,"black",0.5);

## X轴刻度
my $Xstep=50;
my $XstepNum=$cycle_num/$Xstep;
my $X_base=$X_len/$cycle_num;
for (my $i=0;$i<=$XstepNum ;$i++) {
	print OUT &svg_mid_txt($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Space+$Y_len+$XnumSize*4/3,$XnumSize,"black",$i*$Xstep);
	print OUT &svg_line($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Y_len+$Space,$left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$i*$Xstep*$X_base,$up_side+$TitleSize+$XspaceLen+$Y_len+$Space-$XspaceLen,"black",0.5);
}

## 底色
#print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space,$up_side+$TitleSize+$XspaceLen,$X_len,$frame_height-$Space-28*$Y_base,"palegreen");   ## good quality
#print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space,$up_side+$TitleSize+$XspaceLen+$frame_height-$Space-28*$Y_base,$X_len,(28-20)*$Y_base,"oldlace");   ## middle quality
#print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space,$up_side+$TitleSize+$XspaceLen+$frame_height-$Space-20*$Y_base,$X_len,(20+1)*$Y_base,"lavenderblush");   ## low quality

## reads中间虚线
print OUT &svg_dashed($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$X_len/2,$up_side+$TitleSize+$XspaceLen,$left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+$X_len/2,$up_side+$TitleSize+$XspaceLen+$frame_height,"gray","1",0.5);

## 碱基质量分布
foreach my $cycle (sort {$a<=>$b} keys %quality) {
	foreach my $qual (sort {$a<=>$b} keys %{$quality{$cycle}}) {
		my $opacity=sqrt($quality{$cycle}{$qual}/50);
		print OUT &svg_rect($left_side+$YtextSize+$YnumSize+$YspaceLen+$Space+($cycle-1)*$X_base,$up_side+$TitleSize+$XspaceLen+($Y_len-$qual*$Y_base),$X_base-1,$Y_base,"blue",$opacity);
	}
}

print OUT &svg_end();
close OUT;

#######################################################################################
my $svg_dir=dirname "$od/$quality_file.svg";
my $svg_name=basename "$od/$quality_file.svg";
`cd $svg_dir && perl /share/nas2/genome/bmksoft/tool/svg2xxx/v1.0/svg2xxx $svg_name`;
#######################################################################################

my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

#######################################################################################

sub sub_format_datetime #Time calculation subroutine
{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub svg_paper (){#&svg_paper(width,height,[color])
	#my $svg_drawer = getlogin()."@".(`hostname`);
	my $svg_drawer = "chenx"."@"."biomarker\.com\.cn";
	chomp $svg_drawer;
	my @svg_x=@_;
	my $line="";
	$line.="<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n";
	$line.="<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n";
	$line.="<svg width=\"$svg_x[0]\" height=\"$svg_x[1]\">\n";
	$line.="<Drawer>$svg_drawer</Drawer>\n";
	$line.="<Date>".(localtime())."</Date>\n";
	if (defined $svg_x[2]) {
		$line.="<rect x=\"0\" y=\"0\" width=\"$svg_x[0]\" height=\"$svg_x[1]\" fill=\"$svg_x[2]\"/>\n";
	}
	return $line;
}

sub max 
{
	my ($x1,$x2)=@_;
	my $max;
	if ($x1 > $x2) {
		$max=$x1;
	}
	else {
		$max=$x2;
	}
	return $max;
}

sub min 
{
	my ($x1,$x2)=@_;
	my $min;
	if ($x1 < $x2) {
		$min=$x1;
	}
	else {
		$min=$x2;
	}
	return $min;
}

sub color_gradient  #datanow,data1,data2,color1,color2
{
	my @svg_x=@_;
	my $out_color;
	if ($svg_x[0] >=$svg_x[2]) {
		$out_color=$svg_x[4];
	}
	elsif ($svg_x[0] <=$svg_x[1]) {
		$out_color=$svg_x[3];
	}
	else {
		my $tmp_red1=&hex2ten(substr($svg_x[3],1,2));
		my $tmp_gre1=&hex2ten(substr($svg_x[3],3,2));
		my $tmp_blu1=&hex2ten(substr($svg_x[3],5,2));
		my $tmp_red2=&hex2ten(substr($svg_x[4],1,2));
		my $tmp_gre2=&hex2ten(substr($svg_x[4],3,2));
		my $tmp_blu2=&hex2ten(substr($svg_x[4],5,2));
		my $new_red=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_red2-$tmp_red1)+$tmp_red1);
		my $new_gre=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_gre2-$tmp_gre1)+$tmp_gre1);
		my $new_blu=int(($svg_x[0]-$svg_x[1])/($svg_x[2]-$svg_x[1])*($tmp_blu2-$tmp_blu1)+$tmp_blu1);
		$new_red=&ten2hex($new_red);$new_red="0$new_red" if(length($new_red)==1);
		$new_gre=&ten2hex($new_gre);$new_gre="0$new_gre" if(length($new_gre)==1);
		$new_blu=&ten2hex($new_blu);$new_blu="0$new_blu" if(length($new_blu)==1);
		$out_color="#$new_red$new_gre$new_blu";
	}
	return $out_color;
}

sub ten2hex  #ten
{
	my $tmp_ten=$_[0];   
	my $hex_value=uc(sprintf("%lx",$tmp_ten));
	return $hex_value;
}

sub hex2ten  #hex
{
	my $tmp_hex=$_[0];
	my $ten_value=0;
	my $tmp_i=0;
	my $tmp_j=0;
	my @tmp_x=split(//,$tmp_hex);
	my %hash_hex=(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,'A',10,'B',11,'C',12,'D',13,'E',14,'F',15);
	for ($tmp_i=@tmp_x-1;$tmp_i>=0 ;$tmp_i--) {
		$ten_value+=$hash_hex{$tmp_x[$tmp_i]}*(16**$tmp_j);
		$tmp_j++;
	}
	return $ten_value;
}

sub svg_polygon  #colorfill,colorstroke,coloropacity,point1,point2,...
{
	my @svg_x=@_;
	my $svg_color=shift(@svg_x);
	my $svg_color2=shift(@svg_x);
	my $svg_trans=shift(@svg_x);
	my $svg_points=join(" ",@svg_x);
	my $line="<polygon fill=\"$svg_color\" stroke=\"$svg_color2\" opacity=\"$svg_trans\" points=\"$svg_points\"/>\n";
	return $line;
}

sub svg_circle  #&svg_circle(x,y,r,color,[info])
{
	my @svg_x=@_;
	my $line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" />\n";
	if (defined $svg_x[4]) {
		$line="<circle r=\"$svg_x[2]\" cx=\"$svg_x[0]\" cy=\"$svg_x[1]\" fill=\"$svg_x[3]\" onclick=\"alert('$svg_x[4]')\" onmousemove=\"window.status='$svg_x[4]'\" />\n";
	}
	return $line;
}

sub svg_txt  #&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_mid_txt #&svg_mid_txt(x,y,size,color,text,[vertical,0/1/2/3]);
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=0;
	}
	my $svg_matrix='';
	if ($svg_x[5]==0) {
		$svg_matrix="1 0 0 1";
	}
	if ($svg_x[5]==1) {
		$svg_matrix="0 1 -1 0";
	}
	if ($svg_x[5]==2) {
		$svg_matrix="-1 0 0 -1";
	}
	if ($svg_x[5]==3) {
		$svg_matrix="0 -1 1 0";
	}
	my $line="<text fill=\"$svg_x[3]\" transform=\"matrix($svg_matrix $svg_x[0] $svg_x[1])\" text-anchor=\"middle\" font-family=\"ArialNarrow-Bold\" font-size=\"$svg_x[2]\">$svg_x[4]</text>\n";
	return $line;
}

sub svg_dashed  #&svg_line(x1,y1,x2,y2,color,"10 5",[width])
{
	my @svg_x=@_;
	my $line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4]\"/>\n";
	if (defined $svg_x[6]) {
		$line="<line x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\" style=\"stroke-dasharray:$svg_x[5];fill:none;stroke:$svg_x[4];stroke-width:$svg_x[6]\"/>\n";
	}
	return $line;
}
sub svg_line  #&svg_line(x1,y1,x2,y2,color,[width])
{
	my @svg_x=@_;
	my $line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	if (defined $svg_x[5]) {
		$line="<line fill=\"$svg_x[4]\" stroke=\"$svg_x[4]\" stroke-width=\"$svg_x[5]\" x1=\"$svg_x[0]\" y1=\"$svg_x[1]\" x2=\"$svg_x[2]\" y2=\"$svg_x[3]\"/>\n";
	}
	return $line;
}

sub svg_rect  #&svg_rest(x,y,width,height,color,[opacity])
{
	my @svg_x=@_;
	if (!defined $svg_x[5]) {
		$svg_x[5]=1;
	}
	my $line="<rect x=\"$svg_x[0]\" y=\"$svg_x[1]\" width=\"$svg_x[2]\" height=\"$svg_x[3]\" fill=\"$svg_x[4]\" opacity=\"$svg_x[5]\"/>\n";
	return $line;
}

sub svg_end  #end
{
	return "</svg>\n";
}

sub help{
	print <<"End.";
	Description:
        Writer  : $Writer
        Data    : $Data
        Version : $ver
        function: base distribution
    Usage:
		-qu		quality file, must be given
		-od		outdir, must be given

		-h		Help document
End.
	exit;
}	
	
