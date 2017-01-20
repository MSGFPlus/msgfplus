Overview
======

MS-GF+ (aka MSGF+ or MSGFPlus) performs peptide identification by scoring 
MS/MS spectra against peptides derived from a protein sequence database. 
It supports the HUPO PSI standard input file (mzML) and saves results in 
the mzIdentML format, though results can easily be transformed to TSV. 
ProteomeXchange supports Complete data submissions using MS-GF+ search results.

MS-GF+ is optimized for a variety of spectral types, i.e., combinations 
of fragmentation method, instrument, enzyme, and experimental protocols. 
It supports a variety of input file formats, including mzML, mzXML, 
Mascot Generic File (mgf), MS2 files, Micromass Peak List files (pkl), 
and Concatenated DTA files (_dta.txt).

Requirements
======

JRE 1.6 or greater
Main memory 2GB or greater (recommended 4GB)

Downloads / Updates
======

* https://github.com/sangtaekim/msgfplus/releases
* http://omics.pnl.gov/software/ms-gf

*Version number notes*

As of [January 20, 2016 (commit 375d462)](https://github.com/sangtaekim/msgfplus/commit/375d462e30cbe460b699091a7d6ba52bc192aba1) the version numbering scheme changed.
Previously the version number was the SVN commit number; git does not have simple commit numbers, so MSGFPlus was changed to a date-based version numbering scheme.

An example: v10282 became v2016.01.20

Installation
======

Unzip MSGFPlus.zip
Place MSGFPlus.jar in any folder

Usage Information
======

Type `java -jar MSGFPlus.jar` for command line arguments
To convert an mzid output file into a tsv file, run `java -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv`
For detailed documentation, see the "doc" subfolder, or visit:
* [GitHub repo HTML help pages](http://htmlpreview.github.io/?https://github.com/sangtaekim/msgfplus/blob/master/doc/index.html)
* https://omics.pnl.gov/software/ms-gf
* https://bix-lab.ucsd.edu/pages/viewpage.action?pageId=13533355

Contact Information
======

PNNL Proteomics [proteomics@pnnl.gov]
Sangtae Kim [sangtae.kim (at) gmail.com]

Publications
======

MS-GF+: Universal Database Search Tool for Mass Spectrometry, Sangtae Kim, Pavel A. Pevzner, 
Nat Commun. 2014 Oct 31;5:5277. doi: 10.1038/ncomms6277.
http://www.ncbi.nlm.nih.gov/pubmed/?term=25358478

Spectral Probabilities and Generating Functions of Tandem Mass Spectra: A Strike against Decoy Databases, Sangtae Kim, Nitin Gupta and Pavel Pevzner,
J Proteome Res. 2008 Aug;7(8):3354-63. doi: 10.1021/pr8001244.
http://www.ncbi.nlm.nih.gov/pubmed/?term=18597511

Source
======

https://github.com/sangtaekim/msgfplus
