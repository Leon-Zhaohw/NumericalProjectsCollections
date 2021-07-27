#!/usr/bin/perl -w

$RENDER_SCRIPT = "./pbrt_scenes/houdini.pbrt";
$outputPath = "";
@velocities = (3.02,3.02,3.02,3.02,3.02,2.6,2.125,2.09,2.01,1.955,1.915,1.875,1.855,1.83,1.81,1.79,1.775,1.765,1.75,1.74,1.73,1.725,1.715,1.715,1.71,1.705,1.705,1.705,1.705,1.705,1.71,1.715,1.72,1.73,1.735,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.815,1.83,1.84,1.85,1.865,1.88,1.89,1.905,1.92,1.935,1.945,1.96,1.97,1.985,1.995,2.005,2.015,2.025,2.035,2.05,2.06,2.065,2.075,2.085,2.095,2.105,2.11,2.12,2.13,2.135,2.145,2.15,2.16,2.17,2.175,2.185,2.19,2.2,2.205,2.215,2.22,2.23,2.235,2.245,2.25,2.255,2.265,2.27,2.275,2.285,2.29,2.3,2.305,2.31,2.315,2.325,2.33,2.335,2.34,2.35,2.355,2.36,2.365,2.37,2.38,2.385,2.39,2.395,2.4,2.405,2.415,2.42,2.425,2.43,2.435,2.445,2.45,2.455,2.46,2.465,2.47,2.475,2.48,2.485,2.49,2.495,2.505,2.51,2.515,2.52,2.525,2.5325,2.5375,2.5425,2.5475,2.5525,2.5575,2.5625,2.5675,2.5725,2.5775,2.5825,2.59,2.595,2.6,2.605,2.61,2.615,2.62,2.625,2.63,2.635,2.64,2.645,2.65,2.6525,2.6575,2.6625,2.6675,2.6725,2.675,2.68,2.685,2.69,2.695,2.7,2.7025,2.7075,2.7125,2.7175,2.72,2.725,2.73,2.735,2.7375,2.7425,2.7475,2.7525,2.755,2.76,2.7625,2.765,2.77,2.7725,2.7775,2.7825,2.785,2.79,2.7925,2.7975,2.8,2.805,2.81,2.8125,2.8175,2.8225,2.83,2.835,2.84,2.845,2.8525,2.8575,2.8625,2.8675,2.875,2.88,2.885,2.89,2.895,2.9,2.905,2.91,2.9125,2.9175,2.92,2.9225,2.925,2.9275,2.93,2.93,2.9325,2.9325,2.9325,2.9325,2.93,2.9275,2.925,2.9225,2.9175,2.915,2.91,2.9025,2.8975,2.89,2.8825,2.8725,2.8625,2.8525,2.8425);

############################################################################
# rendering subroutine
############################################################################
sub renderFrame {
    my($i) = @_;
    
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
    
    # compute the needed angle
    $angle = -43.55 + ($i - 3) * $velocities[$i];

    # make a copy of the pbrt script unique to this instance
    $thisPBRT = "render_temp." . $number . ".pbrt";
    $copy = "cp ".$RENDER_SCRIPT." ".$thisPBRT;
    system($copy);
    
    # build the final filename
    #$output = "./renders/frame.".$number.".exr"
    
    $output = $outputPath."\\/frame.".$number.".tga";
    
    $input  = "..\\/data\\/houdini_output\\/textures\\/height.".$number.".field3d.gz";
    $geometry = "..\\/data\\/houdini_output\\/surfaces\\/surface.".$number.".pbrt";
    
    $thisGeom = "surface_temp." . $number . ".pbrt";
    $copyGeom = "cp ".$geometry." ".$thisGeom;
    
    system($copyGeom);
    system("perl -pi -e 's/##OUTPUT##/" . $output. "/g' ". $thisPBRT);
    system("perl -pi -e 's/##INPUT##/" . $input. "/g' ". $thisPBRT);    
    system("perl -pi -e 's/##GEOMETRY##/" . $thisGeom. "/g' ". $thisPBRT);
    system("perl -pi -e 's/##ANGLE##/" . $angle . "/g' ". $thisPBRT);
    
    # run pbrt
    #system("pbrt --ncores 1 ".$thisPBRT);
    system("pbrt ".$thisPBRT);
    
    # delete the pbrt file
    system("rm ".$thisPBRT);
    system("rm ".$thisGeom);
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

$outputPath = "houdini_renders";
system("mkdir ".$outputPath);

for ($frame = $startFrame; $frame <= $endFrame; $frame++)
{
    renderFrame($frame);
}
