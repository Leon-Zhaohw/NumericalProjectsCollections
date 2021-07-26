This directory contains the source code for the paper

  D. Heck and T. Schloemer and O. Deussen, "Blue Noise Sampling with
  Controlled Aliasing", ACM Trans. Graph., 2013, in press

If you have questions about the paper or the included source code, please 
email

	Daniel Heck <dheck@gmx.de>

The main program is called "targetrdf", and it attempts to construct a point
set matching a given "target RDF" (or a given target spectrum). As a simple
example, here is how to generate a step blue noise point set, as presented in
the paper:

	./targetrdf --steprp --crit 0.606 --out step.txt


  --steprp      

        Chooses a radial power spectrum with a step characteristic

  --crit 0.606  

        Sets the position of the step; 0.606 is the numerical value for v_max,
        the maximal realizable step position

  --out step.txt

        Writes the final point set to 'step.txt'


There are many additional parameters that can be used to influence the
optimization process. A full list is shown by "targetrdf --help".


  --npts n
  --seed s

        Specifies the number of points to use. The points are initialized to
        random positions, and the seed for the RNG can be set to different
        values to generate a family of similar point sets.

  --in file

        Load initial point set from 'file'

  --reference file

        Loads a point set from 'file' and tries to construct a new point set
        with the same RDF.

  --rdf file
  --spectrum file
  --rp-as-rdf file

        Loads an RDF or power spectrum from 'file' and tries to construct a
        matching point set. 'rp-as-rdf' can be used to generate "dual" point
        sets by reinterpreting a power spectrum as an RDF.

  --peakrp
  --peak power
  --peaksmooth sigma
  --crit

        These paramters are used to generate the single-peak blue noise point
        sets discussed in the paper.

  --nbins n
  --smoothing sigma

        Specifies the number of bins and the amount of smoothing to use when
        calculating RDFs.

FILE FORMATS
------------

The default file formats used by 'targetrdf' are very simple.

Point sets:

  - the first line contains the number of points
  - the remaining lines contain the x/y coordinates of the points 
    in the unit square

        4096
        0.162776 0.74937
        0.0798951 0.880201
        0.58569 0.782718
        0.697048 0.380065
        0.883548 0.747492
        0.432829 0.344747
        ...

RDF/Power spectrum:

  - the file contains ordinate/abscissa pairs:

        6.10352e-05 0.00088529
        0.000183105 0.00107669
        0.000305176 0.00145942
        ...
        0.499695 1.0007
        0.499817 1.00067
        0.499939 1.00064

  - it is assumed that ordinate spacing is uniform

