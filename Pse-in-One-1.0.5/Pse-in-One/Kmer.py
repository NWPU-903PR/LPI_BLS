if __name__ == '__main__':
    import sys
    sys.dont_write_bytecode = True
    import argparse
    import subprocess
    from argparse import RawTextHelpFormatter

    parse = argparse.ArgumentParser(description="This is kmer module for generate kmer vector.",
                                    formatter_class=RawTextHelpFormatter)
    parse.add_argument('inputfile',
                       help="The input file in FASTA format.")
    parse.add_argument('outputfile',
                       help="The output file stored results.")
    parse.add_argument('alphabet', choices=['DNA', 'RNA', 'Protein'],
                       help="The sequence type.")

    parse.add_argument('-k', type=int, default=2,
                       help="The k value of kmer.")
    parse.add_argument('-r', default=0, type=int, choices=[1, 0],
                       help="Whether consider the reverse complement or not.\n"
                            "1 means True, 0 means False. (default = 0)")
    parse.add_argument('-f', default='tab', choices=['tab', 'svm', 'csv'],
                       help="The output format (default = tab).\n"
                            "tab -- Simple format, delimited by TAB.\n"
                            "svm -- The libSVM training data format.\n"
                            "csv -- The format that can be loaded into a spreadsheet program.")
    parse.add_argument('-l', default='+1', choices=['+1', '-1'],
                       help="The libSVM output file label.")

    args = parse.parse_args()

    cmder = ' '.join(sys.argv[1:])
    cmd = 'python Kmerr.pyc' + ' ' + cmder
    subprocess.call(cmd, shell=True)