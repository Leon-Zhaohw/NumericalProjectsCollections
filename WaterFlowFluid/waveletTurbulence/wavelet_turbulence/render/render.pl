#!/usr/bin/perl -w

$XRES = 75 * 6;
$YRES = 100 * 6;

$SAMPLES = 1; # preview
#$SAMPLES = 24; # high quality

$RENDER_SCRIPT = "scene.pbrt";

############################################################################
# rendering subroutine
############################################################################
sub renderFrame {
  my($i,$size) = @_;

  # get a zero padded version of the filename
  my $number = $i;
  if ($i < 1000) {
    $number = "0" . $number;
  }
  if ($i < 100) {
    $number = "0" . $number;
  }
  if ($i < 10) {
    $number = "0" . $number;
  }

  # make a copy of the pbrt script unique to this instance
  $thisPBRT = "render." . $number . ".pbrt";
  $copy = "cp ".$RENDER_SCRIPT." ".$thisPBRT;
  system($copy);
  system("perl -pi -e 's/##FRAMENUMBER##/" . $number . "/g' ". $thisPBRT);
  system("perl -pi -e 's/##XRES##/" . $XRES . "/g' ". $thisPBRT);
  system("perl -pi -e 's/##YRES##/" . $YRES . "/g' ". $thisPBRT);
  system("perl -pi -e 's/##SAMPLES##/" . $SAMPLES. "/g' ". $thisPBRT);
  system("perl -pi -e 's/##SIZE##/" . $size . "/g' ". $thisPBRT);

  # gunzip the volume file
  system("gunzip ../pbrt/density_".$size."_".$number.".pbrt.gz");

  system("chmod +rxw ".$thisPBRT);
  system("chmod +rxw ".$thisPBRT.".bak");

  # run pbrt
  system("/Users/sinithue/devel/pbrt/pbrt_vn_080106/bin/pbrt ".$thisPBRT);

  # rezip the volume file
  system("gzip ../pbrt/density_".$size."_".$number.".pbrt");

  # delete the pbrt file
  system("rm ".$thisPBRT);
  system("rm ".$thisPBRT.".bak");
}

############################################################################
# main
############################################################################
$argcnt = $#ARGV;
if ($argcnt != 1)
{
  print "usage: render.pl <start_frame #> <end_frame #>\n";
  exit(0);
}

# get the range
$startFrame = $ARGV[0];
$endFrame = $ARGV[1];

for ($frame = $startFrame; $frame <= $endFrame; $frame++)
{
  renderFrame($frame, "big");
  #renderFrame($frame, "small");
}
