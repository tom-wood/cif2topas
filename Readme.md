# cif2topas v0.1.1
`cif2topas` is a command line tool that reads a .cif file and outputs a text file with the relevant structural information formatted in a way that may be used in a TOPAS Academic structure.
## Installation
The only requirement is a working installation of Python (>3.6)
## Usage
The command line syntax is
```python fname_in [fname_out]```
where `fname_in` is the name of the cif file; `fname_out` is the name of the TOPAS-consistent output file (if not defined then it takes a similar name to `fname_in` and writes to the same location).
## Issues
Report any issues on github.
## Changes
### v0.1.1
- Fixed lack of atomic information after occ keyword