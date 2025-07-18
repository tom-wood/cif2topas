# cif2topas v0.1.2
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
### v0.1.2 (16/07/2025)
- Fixed non-supported escape character in regex
- Fixed when U values are given as a single full stop character
- Fixed when string values are preceded by only empty whitespace on a new line
### v0.1.3 (17/07/2025)
- Now enabled to deal with authors that have single quotation marks in their names