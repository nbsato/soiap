#!/bin/csh

set n = 1
set vmin = 1.0
set vmax = 1.0
set input = "input.stishovite"
set file1 = "vscale"
set file2 = "input"
set output = "output"

set dv = `echo "($vmax - $vmin)/$n" | bc -l`
#echo $dv

#@ i = 0
@ i = 1
while ($i <= $n)
    echo "i = $i"
    set i4 = `printf "%04d" $i`
    set v = `echo "$vmin + $dv*$i" | bc -l`
#    echo $v
    echo "volume_scale $v" > $file1
    set input_i = $file2.$i4
    set output_i = $output.$i4
    cat $input $file1 > $input_i
    ../../src/opt2_o $input_i >& $output_i
    @ i = $i + 1
end

rm -f $file1
