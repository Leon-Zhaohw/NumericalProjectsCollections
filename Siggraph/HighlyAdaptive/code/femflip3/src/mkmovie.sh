#!/bin/bash
rm -rf mitsuba_*.mp4
ffmpeg -qscale 1 -r 60 -b 9600 -i mitsuba/img/%d_scene_o.png mitsuba_o.mp4
ffmpeg -qscale 1 -r 60 -b 9600 -i mitsuba/img/%d_scene_w.png mitsuba_w.mp4
ffmpeg -qscale 1 -r 60 -b 9600 -i mitsuba/img/%d_scene_p.png mitsuba_p.mp4