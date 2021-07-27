#!/usr/bin/perl -w

$RENDER_SCRIPT = "./pbrt_scenes/PhysBAM.pbrt";
$outputPath = "./PhysBAM_renders/";

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
    
    # make a copy of the pbrt script unique to this instance
    $thisPBRT = "render_temp." . $number . ".pbrt";
    $copy = "cp ".$RENDER_SCRIPT." ".$thisPBRT;
    system($copy);
    
    # build the final filename
    $output = ".\\/PhysBAM_renders\\/full.".$number.".tga";
    $input  = "..\\/data\\/PhysBAM_output\\/textures\\/height.".$number.".field3d.gz";
    $geometry = "..\\/data\\/PhysBAM_output\\/surfaces\\/surface.".$number.".pbrt";

    $thisGeom = "surface_temp." . $number . ".pbrt";
    $copyGeom = "cp ".$geometry." ".$thisGeom;
    system($copyGeom);
    
    system("perl -pi -e 's/##OUTPUT##/" . $output. "/g' ". $thisPBRT);
    system("perl -pi -e 's/##INPUT##/" . $input. "/g' ". $thisPBRT);    
    system("perl -pi -e 's/##GEOMETRY##/" . $thisGeom. "/g' ". $thisPBRT);
    
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

system("mkdir ".$outputPath);

for ($frame = $startFrame; $frame <= $endFrame; $frame++)
{
    renderFrame($frame);
}
