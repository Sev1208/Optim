#!/bin/csh

##################################### cmd input ##########################
if (($#argv < 1)||($#argv > 5)) then
  echo "Invalid arguments ; correct syntax : SF_map_data.csh data1 data2 ..."
  echo "           data1,2: data type 1 to 5 (1-GPS, 2-InSAR, 3-levelling, 4-tilt, 5-gravi)"
  echo ""
  exit
endif

############################# initial cleaning ###########################
#\rm  ./Out/*.ps ./tmp/*
\rm  ./tmp/*
############################Variable Definitions##########################
set outdir=./Out
set tmpdir=./tmp
set sources=sources.txt
set j=-Jx0.0015
set r=-R0/4000/0/6000
set b=-B1000WeSn
set ok="-O -K"
set az=90

############################Plot the data#################################
foreach data ($argv)
  set out=$outdir/mogi_Plot_${data}_${az}.ps
  psbasemap $j $r -K $b -P > $out
  switch ($data)
  
  case 1:
    awk '{if (NR > 8) print $2,$3,$5,$6}' data_gps_mogi.dat > $tmpdir/tmp.gps
    gmt psvelo $tmpdir/tmp.gps $j $r -W2,black -Se250/0.9/5 $ok >> $out
	gmt psxy $sources $j $r -Gred -W1,black $ok >> $out
	gmt psvelo $j $r -W2,black -Se250/0.9/5 $ok -N <<fin >> $out
4500 1500 0.001 0
fin
    gmt pstext -F+f11p,Helvetica-Bold,black $r $j -O -N <<fin >> $out 
5000 1000 1 mm
fin
	gmt psconvert -A -Tj $out
    breaksw
	
  case 2:
    awk '{if (NR > 6) print $2,$3,$5}' data_insar_mogi.dat > $tmpdir/tmp.insar
	gmt xyz2grd $tmpdir/tmp.insar -G$tmpdir/tmp_insar.grd $r -I50
    set tmin=`gmtinfo -C $tmpdir/tmp.insar | awk '{print $5}'`
    set tmax=`gmtinfo -C $tmpdir/tmp.insar | awk '{print $6}'`
    set tstep=`gmtinfo -C $tmpdir/tmp.insar | awk '{if ((($6+$5)/10)<0) print -(($6+$5)/10); else print (($6+$5)/10)}'`
	echo 'tmin=' $tmin, 'tmax=' $tmax
	gmt makecpt -T$tmin/$tmax/$tstep -Z -I -Crainbow > $tmpdir/tmp_$data.cpt
	#gmt makecpt -T-2.2/0/0.0005  -Z -I -Crainbow > $tmpdir/tmp2_$data.cpt
	#gmt makecpt -T-2.2/0/0.2  -Z -I -Crainbow > $tmpdir/tmp3_$data.cpt
	gmt grdimage $tmpdir/tmp_insar.grd $j $ok -C$tmpdir/tmp_$data.cpt >> $out
	gmt psxy $sources $j $r -Gred -W1,black $ok >> $out
	gmt psscale -D7/6/6/0.3 -C$tmpdir/tmp_$data.cpt $ok -L >> $out
	gmt pstext -F+f11p,Helvetica-Bold,black $r $j -O -N <<fin >> $out 
5000 1500 mm
fin
	gmt psconvert -A -Tj $out
    breaksw
	
  case 3:
    awk '{if (NR > 6) print $2,$3,$5}' data_level_mogi.dat > $tmpdir/tmp.level
    #awk '{if (NR > 1) print $1,$2,$3}' mogi_data_out_uz > $tmpdir/tmp.uz
	gmt xyz2grd $tmpdir/tmp.level -G$tmpdir/tmp_level.grd $r -I50
    set tmin=`gmtinfo -C $tmpdir/tmp.level | awk '{print $5}'`
    set tmax=`gmtinfo -C $tmpdir/tmp.level | awk '{print $6}'`
    set tstep=`gmtinfo -C $tmpdir/tmp.level | awk '{if ((($6+$5)/10)<0) print -(($6+$5)/20); else print (($6+$5)/20)}'`
	#echo 'tmin=' $tmin, 'tmax=' $tmax
	gmt makecpt -T$tmin/$tmax/$tstep -Z -Cjet > $tmpdir/tmp_$data.cpt
	#gmt makecpt -T-2.2/0/0.01 -Z -Cseis > $tmpdir/tmp2_$data.cpt
	#gmt makecpt -T-2.2/0/0.2 -Z -Cseis > $tmpdir/tmp3_$data.cpt
	gmt grdimage $tmpdir/tmp_level.grd $j $ok -C$tmpdir/tmp_$data.cpt >> $out
	gmt psxy $sources $j $r -Gred -W1,black $ok >> $out
	gmt psscale -D7/6/6/0.3 -C$tmpdir/tmp_$data.cpt -O -L >> $out
    # gmt psvelo $tmpdir/tmp.level $j $r -Sn150 -W3,black $ok >> $out
    # gmt psvelo <<fin $j $r -Sn150 -W3,black $ok -N >> $out
# 5000 1500 0 0.002 
# fin
	# gmt pstext -F+f11p,Helvetica-Bold,black $r $j $ok -N <<fin >> $out 
# 5000 1000 2 mm
# fin
	# gmt pstext -F+f11p,Helvetica-Bold,black $r $j -O -N <<fin >> $out 
# 5000 6300 uz (mm)
# fin
	gmt psconvert -A -Tj $out
    breaksw
	
  case 4:
    awk '{if (NR > 7) print $2,$3,-$5,-$6}' data_tilt_mogi.dat > $tmpdir/tmp.tilt
    awk '{if (NR > 6) print $2,$3,$5}' data_level_mogi.dat > $tmpdir/tmp.level
	gmt xyz2grd $tmpdir/tmp.level -G$tmpdir/tmp_level.grd $r -I50
    set tmin=`gmtinfo -C $tmpdir/tmp.level | awk '{print $5}'`
    set tmax=`gmtinfo -C $tmpdir/tmp.level | awk '{print $6}'`
    set tstep=`gmtinfo -C $tmpdir/tmp.level | awk '{if ((($6+$5)/10)<0) print -(($6+$5)/20); else print (($6+$5)/20)}'`
	gmt makecpt -T$tmin/$tmax/$tstep -Z -Cjet > $tmpdir/tmp_level.cpt
	gmt grdimage $tmpdir/tmp_level.grd $j $ok -C$tmpdir/tmp_level.cpt >> $out
	gmt psscale -D7/6/6/0.3 -C$tmpdir/tmp_level.cpt $ok -L >> $out
	gmt psxy $sources $j $r -Gred -W1,black $ok >> $out
	gmt psxy $tmpdir/tmp.tilt $j $r $ok -Sx0.2 -W1 >> $out
    gmt psvelo $tmpdir/tmp.tilt $j $r -Se400000/0.9/5 -W2 $ok >> $out
	gmt psvelo <<fin $j $r -Se400000/0.9/5 -W2 $ok -N >> $out
4500 1500 0.000002 0
fin
	gmt pstext -F+f11p,Helvetica-Bold,black $r $j $ok -N <<fin >> $out 
5000 1000 2 microrad
fin
	gmt pstext -F+f11p,Helvetica-Bold,black $r $j -O -N <<fin >> $out 
5000 6300 uz (mm)
fin
	gmt psconvert -A -Tj $out
    breaksw
	
  case 5:
    awk '{if (NR > 1) print $2,$3,$5}' data_gravi_mogi.dat > $tmpdir/tmp.gravi
    set tmin=`gmtinfo -C $tmpdir/tmp.gravi | awk '{print $5}'`
    set tmax=`gmtinfo -C $tmpdir/tmp.gravi | awk '{print $6}'`
    set tstep=`gmtinfo -C $tmpdir/tmp.gravi | awk '{if ((($6+$5)/10)<0) print -(($6+$5)/10); else print (($6+$5)/10)}'`
	echo 'tmin=' $tmin, 'tmax=' $tmax, 'tstep=' $tstep
	gmt makecpt -T$tmin/$tmax/$tstep -Z -I -Crainbow > $tmpdir/tmp_$data.cpt
	#gmt makecpt -T-0.08/-0.010/0.005 -Z -I -Chot > $tmpdir/tmp_$data.cpt
	gmt psxy $tmpdir/tmp.gravi $r $j $ok -Sc0.5 -C$tmpdir/tmp_$data.cpt >> $out
	gmt psxy $sources $j $r -Gred -W1,black $ok >> $out
	gmt psscale -D7/6/6/0.3 -C$tmpdir/tmp_$data.cpt $ok -L >> $out
	gmt pstext -F+f11p,Helvetica-Bold,black $r $j -O -N <<fin >> $out 
5500 1500 microgals
fin
	gmt psconvert -A -Tj $out
    breaksw
  endsw
end

##########################################################################