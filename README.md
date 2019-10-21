# NEBbaschanger
Automate retriveal of primer sequences from NEB website
by NEB, they were not interested in giving me this information (big surprise). And researchers tend to want to stick to
what works, thus the implementation of Selenium on a Java driven website.

Inputs: 1. DNA fasta file that encodes the protein to be mutated.
                2. List of Amino acid changes using three letter codes (Sorry, I havn't added dictionary for one letter codes)
                The list should be a two column CSV with no headers as follows:
                18, Gln
                221, Tyr
                331, Arg

The first column contains the amino acid position to be mutated. The second position contains the three letter code you
want to mutate to. Again, sorry you can't use one letter codes.
The nucleotides to be changed to are codon optimized for expression in E.coli. I can add this as an alternative input if
there is any interest. The amino acids to change each position to are also hard coded as follows:
Leu', 'Val', 'Ile', 'Met', 'Phe', 'Trp', 'Ser', 'Thr', 'Tyr', 'Asn', 'Gln', 'Asp', 'Glu', 'Lys', 'Arg', 'His'
This can also be an optional input if anyone is interested.
The output contains temporary files stored in a newly created "Temp" directory. This is to ensure that we get an output in
case of getting booted off the website, though this has not been a problem so far. Final output is a list of all primers as well
as two seperate lists for forward and reverse primers. The primers are listed by the amino acid position that is being changed followed
by the amino acid it is changed to.
