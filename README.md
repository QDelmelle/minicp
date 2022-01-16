This work compares 2 implementations of the LNS-FFPA search algorithm to solve the DARP in minicp.
One method uses Insertion Sequence Variables (minicp/minicpISV/src/main/java/minicp/examples/DARP_LNSFFPA_Seqvar), 
the other (minicp/minicpISV/src/main/java/minicp/examples/DARP_RK) doesn't.
minicp/minicpISV/src/main/java/minicp/examples/DARPtest.java contains scripts to compare the 2 methods.

The ISV itself in implemented in minicp/minicpISV/src/main/java/minicp/engine/core/
Some ISV constraints have been implemented in minicp/minicpISV/src/main/java/minicp/engine/constraints/: First, Last, Precedence, TransitionTimes.
Unit tests for these contraints can be found in minicp/minicpISV/src/test/java/minicp/engine/constraints/

An approach to the dynamic DARP using ISVs has also been explored in minicp/minicpISV/src/main/java/minicp/examples/DynamicDARP.java,
but this is still a work in progress.

