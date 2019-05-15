#!perl
#by deduct the linear partion in the displacement vecotr, now this script support to uniformly distribute the images by turn the $uniform to True 20190515
#remove the check_d function, add opt_d function, this is a rather easy way to optimize the initial path based on hard sphere model ponychen 20190509
#add a sortatoms function ,sort the atoms if the atoms oder changes . but not works for stacking faults!!! ponychen 20190502.
#i like model systems in MS but it's awful to make neb path by using nebmake.pl under linux. sometimes you need further modulate the path.thus i
#wrote this script to create neb path in the MS. make life better.....
#usage: 1.drag your initial and final structure into MS
#       2.set the values $inifile $finfile $images_num $dmin
#       3.run it and click in the status.txt, it will tell you the atoms meeting too closely!! you may futher modurate these images by yourself!!
#caution: please make sure the cell and atom orders are exactly right between initial and final structures.
#author: ponychen  20190426
#email: 18709821294@outlook.com
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

#chage following values
my $inifile = "00"; #your initial xsd filename
my $finfile = "06"; #your final xsd filename
my $images_num = 5; #how much images would you create
my $dmin = 2.3;       #belowe which the atoms are thought to be close, unit angstrom
my $sort = "False"; # if True, sored the final structure by the initial structure
my $dd = 0.1; #the step size in optmize
my $dc = 0.0; #if the distance between atom in images and in initial atom is bigger than 1.0, this atom are thought to be active, not matrix or frozen 
my $uniform = "False";#if True, the images are set to uniformly distribute along the transition path. but this not always well.

#do not change below codes unless you konw what yuou are doing!!!
my @stepx;
my @stepy;
my @stepz;
my @unitstepx;
my @unitstepy;
my @unitstepz;

my $inidoc = $Documents{"$inifile.xsd"};
my $findoc = $Documents{"$finfile.xsd"};

my $iniatoms = $inidoc->UnitCell->Atoms;
my $finatoms = $findoc->UnitCell->Atoms;     

sub opt_d {
    #this function optimize the atoms that are too close based on hard sphere model
    my ($a, $b, @c, @d, @e) = @_;
    my $dmin2 = $dmin**2;
    for(my $m=0; $m<$a->Count; ++$m){
        my $de = 0;my $adx = 0; my $ady = 0; my $adz = 0; my $theta = 0; my $newadx = 0; my $newady = 0; my $newadz = 0; 
        for(my $n=0; $n<$a->Count; ++$n){
            if( $n == $m ) {next};
            my $dx = @$a[$m]->X - @$a[$n]->X;
            my $dy = @$a[$m]->Y - @$a[$n]->Y;
            my $dz = @$a[$m]->Z - @$a[$n]->Z;
            my $ds = $dx**2 + $dy**2 + $dz**2;
            if ($ds < $dmin2){
                $de += 1; 
                $adx += $dx;
                $ady += $dy;
                $adz += $dz;
                }};
        $newadx = $adx - $adx*$c[$m];
        $newady = $ady - $ady*$d[$m];
        $newadz = $adz - $adz*$e[$m];
        my $itr = 0;
        my $xx = @$a[$m]->X - @$b[$m]->X;
        my $yy = @$a[$m]->Y - @$b[$m]->Y;
        my $zz = @$a[$m]->Z - @$b[$m]->Z;
        my $dv = sqrt($xx**2+$yy**2+$zz**2);
        while( $dv > 1.0 and $de > 0 ){
            if($uniform eq "True"){
            @$a[$m]->X +=  $newadx*$dd;
            @$a[$m]->Y +=  $newady*$dd;
            @$a[$m]->Z +=  $newadz*$dd;}
            else{
            @$a[$m]->X +=  $adx*$dd;
            @$a[$m]->Y +=  $ady*$dd;
            @$a[$m]->Z +=  $adz*$dd;}
            $de  = 0; $adx = 0; $ady = 0; $adz = 0;
            for(my $n=0; $n<$a->Count; ++$n){ 
                if($n == $m){next};             
                my $dx = @$a[$m]->X - @$a[$n]->X;
                my $dy = @$a[$m]->Y - @$a[$n]->Y;
                my $dz = @$a[$m]->Z - @$a[$n]->Z;
                my $ds = $dx**2 + $dy**2 + $dz**2;
                if ($ds < $dmin2){
                   $de += 1; 
                   $adx += $dx;
                   $ady += $dy;
                   $adz += $dz;
                }};
            $itr += 1;
            if($itr > 100){
                last}
                }
                }}
                

sub sortatoms {
    my ($k, $l) = @_;
    my $sorteddoc = Documents->New("sorted.xsd");
    $sorteddoc->CopyFrom($inidoc);
    my $sortedatoms = $sorteddoc->UnitCell->Atoms;
    for(my $i=0; $i<$sortedatoms->Count; ++$i){
        my $rc = 5;
        my $atomnew = @$l[0];
        foreach my $atoml (@$l){
            my $r = (@$k[$i]->X - $atoml->X)**2 + (@$k[$i]->Y - $atoml->Y)**2 + (@$k[$i]->Z - $atoml->Z)**2;
            if ($r < $rc and @$k[$i]->ElementSymbol == $atoml->ElementSymbol){
                $rc = $r;
                $atomnew = $atoml;
                }
        }
        @$sortedatoms[$i]->X = $atomnew->X;
        @$sortedatoms[$i]->Y = $atomnew->Y;
        @$sortedatoms[$i]->Z = $atomnew->Z;
        
        }}
        
if ( $sort eq "True" ){

sortatoms($iniatoms, $finatoms);
$finatoms = $Documents{"sorted.xsd"}->UnitCell->Atoms;
$finfile = "sorted";
       }

#get the stepsize between images            
for(my $i=0; $i<$iniatoms->Count; ++$i){
    $stepx[$i] = (@$finatoms[$i]->X - @$iniatoms[$i]->X)/($images_num + 1);
    $stepy[$i] = (@$finatoms[$i]->Y - @$iniatoms[$i]->Y)/($images_num + 1);
    $stepz[$i] = (@$finatoms[$i]->Z - @$iniatoms[$i]->Z)/($images_num + 1);
    #get the unit step vector of each atom
    my $distance = sqrt($stepx[$i]**2 + $stepy[$i]**2 + $stepz[$i]**2);
    $unitstepx[$i] = $stepx[$i]/$distance;
    $unitstepy[$i] = $stepy[$i]/$distance;
    $unitstepz[$i] = $stepz[$i]/$distance;
    }

#create images
for(my $j=1; $j<=$images_num; ++$j){
    my $innerdoc = Documents->New("$j.xsd");
    $innerdoc->CopyFrom($inidoc);
    my $inneratoms = $innerdoc->UnitCell->Atoms;
    for(my $k=0; $k<$inneratoms->Count; ++$k){
        @$inneratoms[$k]->X = @$iniatoms[$k]->X + $j*$stepx[$k];
        @$inneratoms[$k]->Y = @$iniatoms[$k]->Y + $j*$stepy[$k];
        @$inneratoms[$k]->Z = @$iniatoms[$k]->Z + $j*$stepz[$k];
    }
    opt_d($inneratoms, $iniatoms, @unitstepx, @unitstepy, @unitstepz);
    }
