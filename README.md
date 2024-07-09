Welcome to the Arboretum.
This package was developed to allow mass-spectrometry instruments
the ability to store and search PSM's (Peptide Sequence Matches)
in real time, during and between mass-spec (M2/M2) acquisitions. 

The forest.py file utilizes several existing python packages
and modules to design data storage "trees" of various formats;
some of which are not truly trees at all. The forest retains
access to each of these data storage forms regardless of which
class is utilized or recommended. 
The forest relies on a single AbstractPSMTree class to outline
the methods a tree has or must custom implement.
The binary search trees have implemented a function to allow
a search matching a range of values, not just a specific value.

The Arborist is our special class designed to handle creating,
shaping, naming, searching, and keeping track of the trees
created during the mass-spec run. The Arborist will handle the
forest, not the client -- the client should never have to
interact directly with the forest or the trees. 

Within the ArboretumTests folder are a few files to test the 
forest and the arborist, and to plot the speed of each tree.
The testers can be reconfigured to test all trees in order,
save the results in seconds, which the performance plotter 
can then graph and display against each other.
