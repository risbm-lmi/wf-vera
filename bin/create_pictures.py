import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-in', '--input',
                        type=str, required=True,
                        help="Input kraken report table")
    
    parser.add_argument('-out', '--output',
                        type=str, required=False,
                        help="Pattern for output")
    
    parser.add_argument('-cutoff', '--cutoff',
                        type=int, default=20,
                        help="Minimal cutoff for kraken2 counts")
    
    return parser.parse_args()
