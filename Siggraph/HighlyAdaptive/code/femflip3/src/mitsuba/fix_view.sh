#!/bin/sh
i=$1
while [ $i -le $2 ]; do
	f="$i"_scene_view.exr
	if [ -s "$f" ]; then
		dummy=1
	else
		echo "Starting scene (particles) $i..."
		mitsuba -Dparticles="$i"_particles.xml -Dtet_filename="$i"_tet.ply -Dtarget="0.5, 0.5, 0.1" -Dorigin="0.2, -1.0, 1.5" -Dup="0, 0, 1" -o "$i"_scene_view.exr particles.xml
		exrtopng "$i"_scene_view.exr img/"$i"_scene_view.png
	fi
	i=`expr $i + 1`
done
