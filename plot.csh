#!/bin/csh

set j=-Jx0.001
set r=-R659500/663500/12732000/12737100
set b=-B1000WeSn
set ok="-O -K"
set out=tmp_plot.ps

psbasemap $j $r -K $b -P > $out
gmt psvelo tmp.gps $j $r -W2,black -Se500/0.9/5 $ok >> $out