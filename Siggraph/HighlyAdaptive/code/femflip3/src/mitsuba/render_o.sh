#!/bin/sh
ulimit -c 0
i=$1
while [ $i -le $2 ]; do
f="$i"_scene.serialized
if [ -s "$f" ]; then
echo "Starting scene (opaque) $i..."
sleep 3
exec 5< camera_"$i".txt
read target <&5
read origin <&5
read up <&5
mitsuba -Dfilename="$i"_scene.serialized -Dtet_filename="$i"_tet.ply -Dtarget="$target" -Dorigin="$origin" -Dup="$up" -o "$i"_scene_o.exr opaque.xml
exrtopng "$i"_scene_o.exr img/"$i"_scene_o.png
i=`expr $i + 1`
else
sleep 10
fi
done
