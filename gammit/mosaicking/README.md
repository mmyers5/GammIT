# Mosaicking
## Dependencies
### MOPEX
See install instructions from the [MOPEX Download page].
  - Use `tcsh` to start a cshell and run any mopex commands
  - Use `chmod +x` on any patched files downloaded
  
## Usage
- See [mopex_listmaker.py] to create lists to be used in mopex. Outputs will be created in the data directory under a subdirectory labeled after the GRB and channel used.
- See [mopex_nl.py] to create the namelist to be used by mopex. Outputs will be created in the data directory under the cdf subdirectory and labeled according to the GRB and channel used.
- See [utils.py] in root directory to define constants (such as the data directory and mopex directory) used throughout these scripts.

To use the command-line mopex, run  `mosaic.pl -n namelist.nl`, where `namelist.nl` is a namelist located in the `cdf` directory of the current working directory.

[MOPEX Download page]: https://irsa.ipac.caltech.edu/data/SPITZER/docs/dataanalysistools/tools/mopex/mopexdownload/
[mopex_listmaker.py]: mopex_listmaker.py
[mopex_nl.py]: mopex_nl.py
[utils.py]: https://github.com/mmyers5/GammIT/blob/master/utils.py
