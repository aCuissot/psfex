# PSFEx

[![Build Status](https://travis-ci.org/astromatic/psfex.svg?branch=master)](https://travis-ci.org/astromatic/psfex)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/astromatic-psfex/badge.svg)](https://scan.coverity.com/projects/astromatic-psfex "Coverity Badge")
[![Documentation Status](https://readthedocs.org/projects/psfex/badge/?version=latest)](http://psfex.readthedocs.io/en/latest/?badge=latest)

a utility that makes PSF models for use with the [SExtractor] program.

Check out the on-line [documentation], the [official web page], and the [user forum].

[SExtractor]: http://astromatic.net/software/sextractor
[documentation]: http://psfex.readthedocs.org
[official web page]: http://astromatic.net/software/psfex
[user forum]: http://astromatic.net/forum/forumdisplay.php?fid=21

Working on:
- basic refactor in 'master'
- Python translation, testing python libs and trying too implement a deep learning based method in branch 'pythonFroScratch'
- refactor the code, the includes, ... in branch 'working_on_opti'
- create a cython wrapper to use psfex in python in branch 'cython'
- use C++ to improve performances in branch 'cpp'
