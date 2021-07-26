#!/bin/sh
ulimit -c 0
i=$1
while [ $i -le $2 ]; do
f="$i"_particles.xml
if [ -s "$f" ]; then
echo "Starting scene $i..."
sleep 3
mitsuba -Dparticles="$i"_particles.xml -Dtet_filename="$i"_tet.ply -Dtarget="0.4, 0.3, 0.25" -Dorigin="-0.7, -1.3, 1" -Dup="0, 0, 1" -o "$i"_scene_view.exr particles.xml
exrtopng "$i"_scene_view.exr img/"$i"_scene_view.png
i=`expr $i + 1`
else
sleep 10
fi
done
