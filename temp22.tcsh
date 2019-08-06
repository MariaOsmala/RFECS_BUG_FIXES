#!/bin/bash 



for y in $prefix.*.rpkm
do 
    ./get_rpkm_sites 50 1 sites2_all $y $chrom_file temp1
    if ( $num_mods == 0 ) then
        cat temp1.list > all_mods_rpkm
    else
        paste all_mods_rpkm temp1.list > temp2
        mv temp2 all_mods_rpkm
    endif
    @ num_mods = $num_mods + 1
    echo $y
done

