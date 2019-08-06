# Version 1.0.3 (2015-08-21 dev)

* change position parameters lag in acc.py, k in kmer.py, lamada and w in pse.py to optional parameter, add their default value.


# Version 1.0.2 (2015-08-16 release)

* Use function make_kmer_list to replace repDNA.util.make_kmer_list in util.py
* Add user_defined_property2 to make the user_indices.txt format more clear
* Fix acc.py used time is negative bug
* Fix acc.py and pes.py argparse bool type bug for -a command (It lead -a command always be true)
* Fix the inconsistent command argument Protein between kmer.py and acc.py

# Version 1.0.1 (2015-06-02)

* Add libSVM output file label.
* Optimize check_args and add const FILENAME and corresponding method name.
* Add cmd prompt.
* Fix default_e warning and pse algs.alphabet == 'PROTEIN' bug.
* Add cmd output used time and the help info about -a parameter.
* Add default indices argument.
* Fix PseKNC default indices bug.
* Change the result to 8 decimal.

# Version 1.0.0 (2015-05-09)

* Pse-in-One package release.