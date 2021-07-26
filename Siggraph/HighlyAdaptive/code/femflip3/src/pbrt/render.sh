#!/bin/sh
i=$1
while [ $i -le $2 ]; do
f="$i"_scene.pbrt
if [ -s "$f" ]; then
echo "Starting scene $i..."
pbrt "$i"_scene.pbrt
exrtopng "$i"_scene.exr img/"$i"_scene.png
#convert -gamma 4.0 "$i"_scene.exr img/"$i"_scene.png
i=`expr $i + 1`
else
sleep 10
fi
done
