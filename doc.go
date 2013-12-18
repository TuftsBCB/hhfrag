/*
This package provides a Go implementation of Kalev and Habeck's HHfrag
algorithm from their 2011 paper in Bioinformatics. Notably, this implementation
uses hhsuite 2.0 with either hhsearch or hhblits from Johann Soeding in lieu of
the older hhsuite 1.5/1.6 version (which requires PSI-BLAST).

Note that this package will likely not be maintained, but it is provided for
anyone who would like to pursue using FragBag fragment libraries with HHfrag.
(In particular, the FragmentMap type implements the fragbag/bow.StructureBower
interface by computing BOWs on each segment and returning their sum as a single
BOW.)

The PDB database used in this package is not the regular PDB normally seen. In
particular, it corresponds to a database that contains both an HHblitz and a
set of PDB structures corresponding to entries in the HHblitz database. This
database can be created with the "build-pdb-hhm-db" shell script included in
this repository (although it has not been maintained). The shell script
requires the pdb2fasta tool (use "go get github.com/TuftsBCB/tools/pdb2fasta")
along with a properly configured HHsuite environment (see Johann Soeding's
handbook).

It is worth reading Kalev and Habeck's paper for additional insight in how this
package works. (In lieu of better documentation that simply does not exist.)
*/
package hhfrag
