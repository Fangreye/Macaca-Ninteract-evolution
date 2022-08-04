'''
This scripts accept the output from selscan and finding outlier sliding windows.

Usage:
    python analyse_iHS.py -i path -a path -o path -w 100000 -t 0.99 -w 100000
                <-i | --input ihs_input_path> 
                <-a | --annotation annotation_path>
                <-o | --output output_path>
                [-t | --threshold float]
                [-w | --winsize int]
                [-h | --help]

arguments:
    -i | --input ihs_input_path
        The path of ihs output file
    -a | --annotation annotation_path
        The path of annotation file
    -o | --output output_path
        Where to write the output of sliding windows
    -t | --threshold float
        The threshold to determine if a snp is iHS extreme
    -w | --winsize
        The size of sliding windows
    -s | --selnrom
        Set if the input file has already normalized by selscan
    -h | --help :
        To show this message

'''

import sys
import logging
import pandas
import getopt

from scipy import stats
from pathlib import Path

def read_file(input_path:Path) -> pandas.DataFrame:
    '''
    Read in input
    '''
    my_data = pandas.read_csv(input_path, header=None,sep='\t')
    my_data.columns = ["chr","pos","1_freq","iHH_1","iHH_2","uiHS"]
    my_data.iloc[:,1] = my_data.iloc[:,1].astype('int')
    my_data.iloc[:,2:] = my_data.iloc[:,2:].astype('float')
    return my_data

def read_sel_norm(input_path):
    '''
    Read input from selscan norm windows
    '''
    my_data = pandas.read_csv(input_path, header=None,sep='\t')
    my_data.columns = ["chr","win_start","win_end","num_snp","prop_ext","rank"]
    my_data['win_end'] = my_data['win_end']-1

    return my_data

def read_anno(input_path:Path) -> pandas.DataFrame:
    '''
    Read in annotation
    '''
    my_data = pandas.read_csv(input_path, header=None, sep='\t')
    my_data.columns = ["chr","source","type","start","end",
                        "qual","strand","n","anno","name","is_Ninteract"]
    return my_data

def process_by_strand(my_data:pandas.DataFrame):
    '''
    Rebuild the start position by its strand
    '''
    for i in range(0,my_data.count().chr):
        if my_data['strand'][i] == "-":
            my_data.iloc[i,3], my_data.iloc[i,4] = \
                my_data.iloc[i,4], my_data.iloc[i,3]
    return my_data.sort_values(['chr','start'])

def get_chromosome_length(my_data) -> dict:
    '''
    Get the position of the last snp, used for sliding windows
    '''
    chromosomes = my_data.chr.unique()
    output_dict = dict.fromkeys(chromosomes,0)

    for chr in chromosomes:
        temp_data = my_data.loc[my_data['chr'] == chr]
        output_dict[chr] = temp_data['pos'].iloc[-1]
    
    return output_dict

def create_windows(chr_length_dict, win_size=100000) -> pandas.DataFrame:
    '''
    Ceate sliding window table
    '''
    total_frame = pandas.DataFrame({
            "chr":pandas.Series(dtype=int),
            "win_start":pandas.Series(dtype=int),
            "win_end":pandas.Series(dtype=int)
        })
    for chr in chr_length_dict:
        this_frame = pandas.DataFrame(
                    list(zip(
                            [chr for i in range(1, chr_length_dict[chr], win_size)],
                            range(1, chr_length_dict[chr], win_size),
                            range(win_size, chr_length_dict[chr], win_size)
                            )
                        ),
                    columns= ["chr","win_start","win_end"]
                )
        total_frame = pandas.concat([total_frame, this_frame], ignore_index=True)
    
    return total_frame

def get_normalized(my_data) -> pandas.Series:
    '''
    Calculate z-score based on chromosomes
    '''
    chromosomes = my_data.chr.unique()
    output_series = pandas.Series(dtype='float64')

    for chr in chromosomes:
        temp_data = my_data.loc[my_data['chr'] == chr]
        output_series = output_series.append(
                        stats.zscore(temp_data["uiHS"]), ignore_index=True
                    )

    return output_series

def get_windows_genes(win_data, anno_data):
    '''
    Detect if there is gene inside sliding windows
    '''
    num_genes = list()
    num_Ninteract = list()
    name_Ninteract = list()

    for _ ,windows in win_data.iterrows():
        temp_data = anno_data[anno_data['chr'] == windows['chr']].start.\
                        between(windows['win_start'],windows['win_end'])
        exact_dat = anno_data[anno_data['chr'] == windows['chr']][temp_data]
        num_genes.append(temp_data.sum())

        num_ni = exact_dat.is_Ninteract.sum()
        num_Ninteract.append(num_ni)

        if num_ni == 0:
            name_Ninteract.append('')
            continue

        name_Ninteract.append(
                ','.join(exact_dat[exact_dat['is_Ninteract'] == 1].name)
            )

    my_info = pandas.DataFrame(
            list(zip(num_genes, num_Ninteract, name_Ninteract)),
            columns=['num_genes', 'num_Ninteract', 'name_Ninteract']
        )

    return my_info

def get_windows_ihs(my_data, win_data, chromosome):
    '''
    Count the total number of snps in sliding window and the number of snps
    with extreme iHS
    '''
    total_snp_list = list()
    total_extreme_list = list()

    for chr in chromosome:
        temp_dat_snp = my_data.loc[my_data['chr'] == chr]
        temp_dat_win = win_data.loc[win_data['chr'] == chr]

        for _ , row in  temp_dat_win.iterrows():
            total_snps = temp_dat_snp[
                            temp_dat_snp['pos'].between(
                                row['win_start'], row['win_end']
                            ) 
                        ].pos.count()
            n_extreme_snps = temp_dat_snp[
                                temp_dat_snp['pos'].between(
                                    row['win_start'], row['win_end']
                                ) 
                        ].is_extreme.sum()
            
            total_snp_list.append(total_snps)
            total_extreme_list.append(n_extreme_snps)
        
    my_info = pandas.DataFrame(
            list(zip(total_snp_list, total_extreme_list)),
            columns=['number_snp', 'number_extreme_snp']
        )

    return my_info

def arg_parser():
    config = {
        'input_file': Path(r'F:\Work\Now_work\2021_05_18.ben_\19.iHS\arctoides_arctoides.total.ihs.out'),
        'anno_path' : Path(r'F:\Work\Now_work\2021_05_18.ben_\00.annotation\rheMac8.ncbiRefSeq.CDS.marked.corrected.out'),
        'log_path': Path(r'F:\Work\Now_work\2021_05_18.ben_\19.iHS\mylog.txt'),
        'output': r'',
        'threshold': stats.norm.ppf(0.99),
        'win_size' : 100000,
        'is_selscan_norm':0
        }

    opts, _ = getopt.getopt(
                sys.argv[1:],
                "i:a:o:t:w:hs",
                ["annotation=", "input=", "output=","threshold=","winsize=","help","selnrom"]
        )

    if not opts:
        print(__doc__)
        sys.exit(0)
    for k,v in opts:
        if k in {"-h","--help"}:
            print(__doc__)
            sys.exit(0)
        if k in {"-i","--input"}:
            config['input_file'] = Path(v)
        if k in {"-a","--annotation"}:
            config['anno_path'] = Path(v)
        if k in {"-o","--output"}:
            config['output'] = v
        if k in {"-t","--threshold"}:
            config['threshold'] = stats.norm.ppf(v)
        if k in {"-w","--winsize"}:
            config['win_size'] = int(v)
        if k in {"-s","--selnrom"}:
            config['is_selscan_norm'] = 1

    return config

def main(config):

    if config['is_selscan_norm'] == 0:
        print("Normal")
        logging.info("Reading {}".format(config['input_file']))
        my_data = read_file(config["input_file"])
        my_data = my_data.assign(z_score_iHS=get_normalized(my_data))
        threshold = config['threshold']
        is_extreme = my_data['z_score_iHS'].apply(lambda x: 1 if abs(x) > threshold else 0)
        my_data = my_data.assign(is_extreme=is_extreme)

        anno_data = read_anno(config["anno_path"])
        anno_data = process_by_strand(anno_data)

        win_size = config['win_size']
        chr_length_dict = get_chromosome_length(my_data)
        windows_table = create_windows(chr_length_dict, win_size)

        iHS_prop = get_windows_ihs(my_data, windows_table, chr_length_dict)
        gene_info = get_windows_genes(windows_table, anno_data)

        total_out = pandas.concat([windows_table, iHS_prop, gene_info], axis=1)

        my_data.to_csv(config['output']+"_snp_ihs.tsv", index=False, sep="\t")
        total_out.to_csv(config['output']+"_windows_ihs.tsv", index=False, sep="\t")
    else:
        print("Alter")
        my_data = read_sel_norm(config["input_file"])
        anno_data = read_anno(config["anno_path"])
        anno_data = process_by_strand(anno_data)
        gene_info = get_windows_genes(my_data, anno_data)

        total_out = pandas.concat([my_data, gene_info], axis=1)

        total_out.to_csv(config['output']+"_windows_ihs.tsv", index=False, sep="\t")

if __name__ == "__main__":
    main(arg_parser())
