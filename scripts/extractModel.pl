#!/usr/bin/perl
##############################################
##############################################
#This program extract the 
#first model that appears in
#the pdb file. If the record MODEL
#does not appear then it will read
#the first sequence of ATOM records.
#Only the CAs are stored in the output file
##############################################
##############################################


$inputName  = $ARGV[0];
$outputName = $ARGV[1]; 
chomp($inputName);
chomp($outputName);

unlink($outputName);

open(IN,$inputName) || die "ERROR - I can't open file:$inputName";
open(OUT,">$outputName") || die "ERROR - I can't open file:$outputName";

$r1  = "";
$r2  = "";
$r3  = "";
$r4  = "";
$r5  = "";
$r6  = "";
$r7  = "";
$r8  = "";
$r9  = "";
$r10 = "";
$r11 = "";
$r12 = "";
$flag = 1;
$X    = 0;
$Y    = 0;
$Z    = 0;
$absoluteResidueNumber = 0;

$_=<IN>;
@pdbHeader = split(/\s+/, $_);
print OUT "MOLECULE\t$pdbHeader[$#pdbHeader]\n";



while(<IN>)
{
 if($flag==1)
 {
     ($r1,$r2,$r3,$r4,$r5,$r6,$r7,$r8,$r9,$r10,$r11,$r12) = split(/\s+/, $_);

     if(!(index($r1,"ATOM")==-1) && !(index($r3,"CA")==-1) )
     {# we are in a model part of the pdb file 
	  $absoluteResidueNumber ++;
	  #print OUT "ATOM\t$r2\t$r3\t$r4\t";
	  if($r5=~ /[a-zA-Z]/)
	  {# There is a chain id in the file
	       $X = $r7;
	       $Y = $r8;
               $Z = $r9; 
	       #print OUT "$absoluteResidueNumber\t$r7\t$r8\t$r9\t$r10\t$r11\t$r12\n";
	   }
	  else
	  {
	      #print OUT "$absoluteResidueNumber\t$r6\t$r7\t$r8\t$r9\t$r10\t$r11\n";
	      $X = $r6;
	      $Y = $r7;
	      $Z = $r8;
	  }
	
	 # $completeLine="ATOM  %4i %4sA%3s A%4iA   %8.3f%8.3f%8.3f                         "; 
          
	  printf OUT "ATOM   %4i%4s  %3s A%4i    %8.3f%8.3f%8.3f                         \n",$r2,$r3,$r4,$absoluteResidueNumber,$X,$Y,$Z;
      }
     if(!(index($r1,"TER")==-1))
     {
	 $flag = 0;
     }
     #print $flag;

 }

}

print OUT "END\n";
close(IN);
close(OUT);

