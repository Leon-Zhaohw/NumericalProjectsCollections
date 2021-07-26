#!/bin/sh
i=$1
while [ $i -le $2 ]; do
	f="$i"_scene_f.exr
	if [ -s "$f" ]; then
		dummy=1
	else
		echo "Starting scene (opaque) $i..."
		exec 5< camera_"$i".txt
		read target <&5
		read origin <&5
		read up <&5
		mitsuba -Dfilename="$i"_scene.serialized -Dtet_filename="$i"_tet.ply -Dtarget="$target" -Dorigin="$origin" -Dup="$up" -o "$i"_scene_f.exr fancy.xml
		exrtopng "$i"_scene_f.exr img/"$i"_scene_f.png
	fi
	i=`expr $i + 1`
done
