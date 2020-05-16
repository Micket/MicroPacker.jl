# MicroPacker

Micropacker is a small program that packs inclusions into a periodic structure while attemptimg to minimize overlap.
It does this by simply discretizing the inclusions, and place them into a periodic voxel grid which minimizes overlap (in one pass).

The brute-force approach was shown to be highly effective at generating realistic WC-Co hardmetal microstructures.
This approach also manages to avoid introducing any unwanted anisotropy and minimal domain size bias.

It is ported from a Python code called CCBuilder, partly because I wanted to try out Julia.
Some functionality of CCBuilder has been left out, and a lot has been greatly generalized to accomodate more types of microstructures.
For more details see "CCBuilder: a software that produces synthetic microstructures of WC-Co cemented carbides" <https://www.sciencedirect.com/science/article/pii/S0263436818303457>

# Usage

Code under heave development, and things might still radically change
