README_istac.txt
J Pillow.
Sept 20, 2007 (version 1).

DESCRIPTION OF THE CODE CONTAINED IN THE ARCHIVE: code_iSTAC_v1.tgz

Demonstrates iSTAC filter-estimation/dimensionality-reduction using
simulated LNP data.  


INSTALLATION
============

Unnpack the archive using a compression utility such as WinZip
(or, from the command line in linux: 'tar -xvzf FILENAME.tgz').

Launch matlab and cd into the directory containing the code
(e.g. '/code_iSTAC/').



USE
===

Examine the script 'test_iSTAC_script.m' for a line-by-line tutorial
on how to use the code contained in this package, which goes through
several simulated examples.

The primary function used for estimating the filters is 'compiSTAC.m'



FUTURE PLANS
============

Still *not* included in this release:

1.  code for estimating the filters under a space-time separability
constraint (i.e., as in the last section of the JOV paper).  The
optimization under this constraint takes a slightly different form,
and I still need to comment and integrate the relevant functions.

2.  Bootstrap code for estimating statistical significance of the
number of filters (i.e., for deciding how many dimensions/filters are
truly present, and not the result of undersampling).

