'''
Collect all dNdS values.

usage: 
    python get_dNdS.py 
                    <-i | --input roh_input_path> 
                    <-o | --outdir output_path>

arguments:
    -i | --input roh_input_path
        The path to store fasta files
    -o | --output output_path
        Where to write the output
'''

import re
import sys
import getopt

from pathlib import Path
from typing import Dict, List, Tuple

def get_information(file_path:Path) -> Dict:

    dN_pattern = re.compile(r'tree length for dN:\W*(?P<dN>\d*\.\d*)')
    dS_pattern = re.compile(r'tree length for dS:\W*(?P<dS>\d*\.\d*)')
    output = dict()
    with file_path.open() as f:
        for line in f:
            if dN_pattern.search(line):
                output['dN'] = dN_pattern.search(line).groupdict()['dN']
            elif dS_pattern.search(line):
                output['dS'] = dS_pattern.search(line).groupdict()['dS']
    return output

def arg_parser():

    config = {
        "mlc_path":"",
        "output":"",
    }

    opts, _ = getopt.getopt(
                sys.argv[1:],
                "i:o:h",
                ["input=","outdir=","help"]
            )

    if not opts:
        print(__doc__)
        sys.exit(0)
    for k,v in opts:
        if k in {"-i","--input"}:
            config['mlc_path'] = Path(v)
        elif k in {"-o","--output"}:
            config['output'] = Path(v)
        elif k in {'-h', '--help'}:
            print(__doc__)
            sys.exit(0)
    return config

def main(config:Dict):

    with config['output'].open('w') as f:
        cnt = 0
        for i in config['mlc_path'].iterdir():
            if not i.suffix ==  '.mlc':
                continue

            cnt += 1
            output = {
                'gene_name':i.stem.split('_')[1]
            }
            output.update(get_information(i))
            
            if len(output) != 3:
                print(i.stem)
                continue

            f.write("{},{},{}\n".format(
                output['gene_name'],output['dN'], output['dS'])
            )

    if cnt == 0:
        print('No mlc file found')
    else:
        print('Processed {} files'.format(cnt))

if __name__ == '__main__':
    main(arg_parser())
