#!/bin/sh
i=$1
while [ $i -le $2 ]; do
f="$i"_scene.png
if [ -s "$f" ]; then
dummy=1
else
echo "Copying..."
j=`expr $i - 1`
echo cp "$j"_scene.png "$i"_scene.png
cp "$j"_scene.png "$i"_scene.png
fi
i=`expr $i + 1`
done
