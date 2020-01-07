
# create a series of regional maps
for pt in sp np q1 q2 q3 q4 ; do

ocnres=<ocnres>    #ie gx1PETM
date=<date>        #yymmdd

# space = 20 : fast, defines how many gridcells are plotted
ncl 'space=20' 'type="reg"' 'expt="'${ocnres}'"' 'plot_type="'${pt}'"' 'date="'${date}'"' plot_global_all.ncl

done


