#!/usr/bin/perl

$nst = 10;
$natom = 48;
$natomperfu = 3;
$angs = 0.529177;
$filepre = 'output';
$pattern = 'tote,fmax,tote/natom,omega';

for ($i = 0; $i <= $nst; $i++) { 
#
  $id = sprintf("%04d", $i);
  $file = $filepre.".".$id;
#  print "$file \n";
#
  open(IN, $file) || die "ecollect.pl : $!";
#  while(<IN>){
  while($l=<IN>){
#      print "$_";
      if ($l =~ /$pattern/){
	  $etot = substr($l,79,20);
	  $vol = substr($l,119,20);
	  if ($i==0){
#	      $e0=$etot;
	      $e0=0.0;
	  }    
	  $e=($etot-$e0)*27.2116*$natomperfu;
	  $v=$vol*$natomperfu;
#	  print "$v $e\n";
	  printf("%16.8f %16.8f\n", $v, $e);
      }
  }
  close(IN)
#
}

__END__
printf "$i \n";
print "$l\n";
print "$. : $l";
print "$. : $_";
