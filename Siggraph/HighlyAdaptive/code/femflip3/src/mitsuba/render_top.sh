#!/bin/sh
ulimit -c 0
i=$1
while [ $i -le $2 ]; do
f="$i"_particles.xml
if [ -s "$f" ]; then
echo "Starting scene $i..."
sleep 3
#/Applications/Mitsuba.app/Contents/MacOS/mitsuba
mitsuba -Dparticles="$i"_particles.xml -Dtet_filename="$i"_tet.ply -Dtarget="0.43, 0.43, 0.1" -Dorigin="-0.27, -0.27, 3" -Dup="0, 0, 1" -o "$i"_scene_top.exr particles.xml
exrtopng "$i"_scene_top.exr img/"$i"_scene_top.png
i=`expr $i + 1`
else
sleep 10
fi
done
