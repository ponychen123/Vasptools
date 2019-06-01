#!perl
#20190601 reconstruction all the code 
#20190601 adding weight function 2/d^-4, it looks wonderfoul
#20190529 fix bugs of 20190525 version
#20190525 adding pbc condition, this version may be the last one..
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
my $inifile = "00";   #your initial xsd filename
my $finfile = "06";   #your final xsd filename
my $images_num = 5;   #how much images would you create
my $dmin = 2.3;       #belowe which the atoms are thought to be close, unit angstrom
my $sort = "False";   #if True, sored the final structure by the initial structure, only expermentlly use
my $dd = 0.3;         #the step size in optmize
my $dc = 0.5;         #if the distance between atom in images and in initial atom is bigger than 0.5, this atom are thought to be active, not matrix or frozen 
my $uniform = "False";#if True, the images are set to uniformly distribute along the transition path. but this not always well.
my $maxitr = 100;     #the maximum iteration for each cycle in optimization.
my $A = 2;            #the amplitude of weight function A/d^-N
my $N = 4;            #the power of weight function A/d^-N

#do not change below codes unless you konw what yuou are doing!!!
my @stepx;
my @stepy;
my @stepz;
my @unitstepx;
my @unitstepy;
my @unitstepz;
my $uspx = \@unitstepx;
my $uspy = \@unitstepy;
my $uspz = \@unitstepz;

my $inidoc = $Documents{"$inifile.xsd"};
my $findoc = $Documents{"$finfile.xsd"};

my $iniatoms = $inidoc->UnitCell->Atoms;
my $finatoms = $findoc->UnitCell->Atoms;

#get the coordination matrix of three bias axis
my @vec;
    $vec[0] = $inidoc->SymmetryDefinition->VectorA->X;
    $vec[1] = $inidoc->SymmetryDefinition->VectorB->X;
    $vec[2] = $inidoc->SymmetryDefinition->VectorC->X;
    $vec[3] = $inidoc->SymmetryDefinition->VectorA->Y;
    $vec[4] = $inidoc->SymmetryDefinition->VectorB->Y;
    $vec[5] = $inidoc->SymmetryDefinition->VectorC->Y;
    $vec[6] = $inidoc->SymmetryDefinition->VectorA->Z;
    $vec[7] = $inidoc->SymmetryDefinition->VectorB->Z;
    $vec[8] = $inidoc->SymmetryDefinition->Vectorc->z;    

sub opt_d {
    #this function optimize the atoms that are too close based on hard sphere model
    my ($a, $b, $c, $d, $e) = @_;
    my $dmin2 = $dmin**2;
    my $de = 1;
    my $itr = 0;
    my @adx; my @ady; my @adz; my @newadx; my @newady; my @newadz; 
    while($de > 0){
        #initialize
        $de = 0;
        for(my $m=0; $m<$a->Count; ++$m){ 
            $adx[$m] = 0; $ady[$m] = 0; $adz[$m] = 0; $newadx[$m] = 0; $newady[$m] = 0; $newadz[$m] = 0;
            for(my $n=0; $n<$a->Count; ++$n){
                if( $n == $m ) {next};
                my $veca = @$a[$m]->FractionalXYZ->X - @$a[$n]->FractionalXYZ->X;
                my $vecb = @$a[$m]->FractionalXYZ->Y - @$a[$n]->FractionalXYZ->Y;
                my $vecc = @$a[$m]->FractionalXYZ->Z - @$a[$n]->FractionalXYZ->Z;
                #performing pbc condition
                if($veca>0.5){$veca-=1};
                if($veca<-0.5){$veca+=1};
                if($vecb>0.5){$vecb-=1};
                if($vecb<-0.5){$vecb+=1};
                if($vecc>0.5){$vecc-=1};
                if($vecc<-0.5){$vecc+=1};
                #transport the displacement vector to Cartersian coordination
                my $dx = $veca*$vec[0]+$vecb*$vec[1]+$vecc*$vec[2];
                my $dy = $veca*$vec[3]+$vecb*$vec[4]+$vecc*$vec[5];
                my $dz = $veca*$vec[6]+$vecb*$vec[7]+$vecc*$vec[8];
                my $ds = $dx**2 + $dy**2 + $dz**2;
                if ($ds < $dmin2){
                    #apply weighting function 2/d^-4
                    my $dss = sqrt($ds);
                    $de += 1; 
                    $adx[$m] += $A/$dss**$N*$dx;
                    $ady[$m] += $A/$dss**$N*$dy;
                    $adz[$m] += $A/$dss**$N*$dz;}
                    };
        if($uniform eq "True"){
             $newadx[$m] = $adx[$m] - $adx[$m]*@$c[$m];
             $newady[$m] = $ady[$m] - $ady[$m]*@$d[$m];
             $newadz[$m] = $adz[$m] - $adz[$m]*@$e[$m];}
             }
        #displace atoms by displacement vector
        for(my $m=0; $m<$a->Count; ++$m){
            my $xx = @$a[$m]->X - @$b[$m]->X;
            my $yy = @$a[$m]->Y - @$b[$m]->Y;
            my $zz = @$a[$m]->Z - @$b[$m]->Z;
            my $dv = sqrt($xx**2+$yy**2+$zz**2);
            #only when the atom is active and meet too closely with others,then this atom will be relaxed
            if($dv > $dc){
                if($uniform eq "True"){
                    @$a[$m]->X +=  $newadx[$m]*$dd;
                    @$a[$m]->Y +=  $newady[$m]*$dd;
                    @$a[$m]->Z +=  $newadz[$m]*$dd;}
                    else{
                         @$a[$m]->X +=  $adx[$m]*$dd;
                         @$a[$m]->Y +=  $ady[$m]*$dd;
                         @$a[$m]->Z +=  $adz[$m]*$dd;}
                         }
        $itr += 1;
        #if iretation steps reaching $maxitr, break this function and you should check relative parameters
        if($itr > $maxitr){
            last}
        }}    
        }        


sub sortatoms {
    #this is a silly function, most times it not works well...
    #it wasn't recommened to use this function at present
    my ($k, $l) = @_;
    my $sorteddoc = Documents->New("sorted.xsd");
    $sorteddoc->CopyFrom($inidoc);
    my $sortedatoms = $sorteddoc->UnitCell->Atoms;
    for(my $i=0; $i<$sortedatoms->Count; ++$i){
        my $rc = 100;
        my $atomnew = @$l[0];
        foreach my $atoml (@$l){
            my $veca = @$k[$i]->FractionalXYZ->X - $atoml->FractionalXYZ->X;
            my $vecb = @$k[$i]->FractionalXYZ->Y - $atoml->FractionalXYZ->Y;
            my $vecc = @$k[$i]->FractionalXYZ->Z - $atoml->FractionalXYZ->Z;
            #performing pbc condition
            if($veca>0.5){$veca-=1};
            if($veca<-0.5){$veca+=1};
            if($vecb>0.5){$vecb-=1};
            if($vecb<-0.5){$vecb+=1};
            if($vecc>0.5){$vecc-=1};
            if($vecc<-0.5){$vecc+=1};
            #transport the displacement vector to Cartersian coordination
            my $dx = $veca*$vec[0]+$vecb*$vec[1]+$vecc*$vec[2];
            my $dy = $veca*$vec[3]+$vecb*$vec[4]+$vecc*$vec[5];
            my $dz = $veca*$vec[6]+$vecb*$vec[7]+$vecc*$vec[8];
            my $r = $dx**2 + $dy**2 + $dz**2;
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
    opt_d($inneratoms, $iniatoms, $uspx, $uspy, $uspz);
    }
