#!/bin/bash

PBRT=/opt/pbrt/bin
export PBRT_SEARCHPATH=$PBRT

if [ $# -lt 1 ] ; then
  echo "usage : ./render.sh [scene]"
  echo "  where scene is 'sta' or 'dyn'"
else
  cd render
  for f in $1*.pbrt.gz ; do
    fn=`basename $f .gz`
    if [ ! -e $fn.exr ]; then
      gunzip $f
      cp $fn volume_$1.pbrt
      gzip $fn
      echo Rendering $fn
      $PBRT/pbrt $1_scene.pbrt
      mv render_$1.exr $fn.exr
    fi
  done
  cd ..
fi