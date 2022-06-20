'''
This script is used for analysing the roh blocks. 
It writes with python 3.7.2

usage: 
    python3 ROH_permutation.py 
                    <-i | --input roh_input_path> 
                    <-a | --annotation annotation_path >
                    [-o | --output output_path] [-s| --stats stats_output_path]
                    [-h | --help]

arguments:
    -i | --input roh_input_path
        The path of roh blocks file
    -a | --annotation annotation_path
        The path of annotation file
    -o | --output output_path
        Where to write the density file
    -s | --stats stats_output_path
        Where to write the stats file
    -h | --help :
        To show this message
'''

import getopt
import random
import sys
import math

from pathlib import Path
from typing import Tuple, Dict, List

stat_output_sequence = [
    'total gene number','contain othergene length','contain othergene blocks',
    'contain Ngene length','contain Ngene blocks','contain none length',
    'contain none blocks','Average length Ngene','Average length OtherGene',
    'Average length None','gene density nonNinteract','gene density Ninteract',
    'Proportion of N interact blocks with all blocks that have genes',
    'Proportion of N interact genes'
]


class Ndict(dict):
    '''
    A modifed dict class with new method
    '''
    def __init__(self):
        super(Ndict,self).__init__()
    
    def judgein(self, keys, value):
        '''
        judge whether a key is in a diction, if not, use that key to index a 
        new list, else append the value to the list.
        '''
        if keys in self:
            self[keys].append(value)
        else:
            self[keys] = [value]

def fast_count(file_path: Path) -> int:
    '''
    count the amount of lines
    '''
    with file_path.open('rb') as f:
        count = 0
        while True:
            data = f.read(0x400000)
            if not data:
                break
            count += data.count(b'\n')
    return count

def binarySearch(arr:List, x:int, left:int, right:int) -> Tuple[int,int]:
    '''
    The binary search implementation
    '''
    mid = int((left + right)/2)

    if right - left > 1:
        if arr[mid] == x:
            return (0,mid)
        elif arr[mid] > x:
            return binarySearch(arr, x, left, mid)
        else:
            return binarySearch(arr, x, mid, right)
    else:
        if x < arr[mid]:
            return (-1,mid)
        if x > arr[mid]:
            return (-2,mid)
        if x == arr[mid]:
            return (0,mid)

def get_annotation_hash(
        annotation_path: Path, exclude_chrs: List =["X","Y","M"] 
    ) -> Tuple[Dict[str,Dict], str]: 
    '''
    Get the coordiantes of the gene annotations
    '''
    anno_hash_dict = dict()
    # anno_hash_length = 0
    permu_str = ""

    with annotation_path.open() as f:
        for line in f:
            info_list = line.replace('\n','').split('\t')

            # filter unwanted annotations
            if "NW" in info_list[0]:
                continue
            if info_list[0] in exclude_chrs:
                continue

            # the information needed for downstream pipeline
            try:
                anno_single_dict = {
                    'start':int(info_list[3]),
                    'end':int(info_list[4]),
                    'is_N':int(info_list[10]),
                    'name':info_list[9]
                }
            except:
                print(line)
                print(info_list)
                sys.exit(1)

            if info_list[6] == "-":
                anno_single_dict['start'] = int(info_list[4])
                anno_single_dict['end'] = int(info_list[3])
            # the unique identifier for single annotation record
            hash_id = "{}_{}".format(
                            anno_single_dict['start'], anno_single_dict['end']
                        )
            # anno_hash_length = anno_hash_length + 1

            # build the quick search table, first index is the chromosome
            if not info_list[0] in anno_hash_dict:
                anno_hash_dict[info_list[0]] = {hash_id :anno_single_dict}
            else:
                anno_hash_dict[info_list[0]].update({hash_id :anno_single_dict})

            # 
            permu_str = permu_str + str(anno_single_dict['is_N'])

    return (anno_hash_dict,permu_str)

def get_supple_dict(anno_hash_dict:Dict) -> Dict:
    '''
    Build the auxiliary list for binary search, return 3 components:
        hash_list: the original hash in anno_hash_dict
        start_sort_list: sorted start positions
        index_list: the index of the start position inside the original 
                    hash dict
    '''
    supple_dict = dict()
    # keys in the first level of anno_hash_dict are chromosomes
    for key in anno_hash_dict.keys():
        # index the hash keys
        anno_hash_dict_key = list(anno_hash_dict[key].keys())
        start_list_org = [
                anno_hash_dict[key][x]['start'] for x in anno_hash_dict_key
            ]

        temp_list = sorted(
                        [p for p in enumerate(start_list_org)] ,
                        key = lambda x:int(x[1])
                    )

        # check if there are genes which have the same start position 
        if len(temp_list) != 1:
            # create Ndict to handle with duplication
            dup_Ndict = Ndict()
            for i in range(0,len(temp_list)-1):
                # temp_list is sorted, so only need to check the neighbours
                if temp_list[i][1] == temp_list[i+1][1] or \
                        temp_list[i][1] == temp_list[i-1][1] :
                    dup_Ndict.judgein(temp_list[i][1],temp_list[i][0])
            # judge the last element, as it can't be handled in the loop
            if temp_list[i+1][1] == temp_list[i][1]:
                dup_Ndict.judgein(temp_list[i+1][1],temp_list[i+1][0])

            # handle with those genes
            if dup_Ndict:
                for i in dup_Ndict:
                    # Order is not necessary here as we only need to tell how  
                    #   many genes does the roh have
                    for index in range(len(dup_Ndict[i])):
                        start_list_org[dup_Ndict[i][index]] = \
                            start_list_org[dup_Ndict[i][index]] + 0.01 * (index +1)
                # recreate the start list
                temp_list = sorted(
                        [p for p in enumerate(start_list_org)] ,
                        key = lambda x:x[1]
                    )

        supple_dict[key] = {
            'hash_list' : anno_hash_dict_key,
            'start_sort_list':[x[1] for x in temp_list],
            'index_list':{y:x for x,y in temp_list}
        }
    
    return supple_dict

def record_analyser(
        roh_path:Path, supple_dict:Dict, anno_hash_dict: Dict, 
        exclude_chrs: List =["X","Y","M"] 
    ) -> List[Dict]:
    '''
    Read in roh records and count their stats.
    '''
    record_collector = list()

    total_lines = fast_count(roh_path)
    sample_counter = list()
    # This is a counter of how many records are processed, just for debugging
    ## cnt = 0
    timer = 0
    with roh_path.open() as f:
        for line in f:
            timer = timer + 1
            sys.stdout.write("Processed: {}/{}\r".format(timer,total_lines))
            sys.stdout.flush()

            # exclude information line and 
            if "#" in line :
                total_lines = total_lines -1
                timer = timer -1
                continue
            temp = line.replace("\n","").split("\t")
            if temp[0] != "RG":
                continue
            if temp[2] in exclude_chrs:
                continue
            ## cnt = cnt+1
            # Initialize the container
            record = {
                'sample':temp[1],
                'chr': temp[2],
                'start':temp[3],
                'end':temp[4],
                'length':temp[5],
                'pos':temp[3],
                'containsgenes':0,
                'containsNinteractgenez':0,
                'num_genes':0,
                'num_ninteractgenes':0
            }

            # sample counter
            sample_counter.append(record['sample'])

            # make alias 
            supple_of_record = supple_dict[record['chr']]

            start_info = binarySearch(
                            supple_of_record["start_sort_list"],
                            int(record['start']),
                            left=0,
                            right=len(supple_of_record["start_sort_list"])
                        )
            end_info = binarySearch(
                            supple_of_record["start_sort_list"],
                            int(record['end']),
                            left=0,
                            right=len(supple_of_record["start_sort_list"])
                        )

            # handle with head and the tail to get the position of real block
            #   binary implementation here always in the right side
            start_pos = start_info[1]+1 if start_info[0] == -2 else start_info[1]

            if end_info[1] == len(supple_of_record["start_sort_list"]) - 1:
                end_pos = len(supple_of_record["start_sort_list"])
            elif end_info[0] == -1:
                end_pos = end_info[1]
            else:
                end_pos = end_info[1] + 1

            # find all genes situated in that block, check whether they are N-
            #   interact
            for index in range(start_pos,end_pos):
                # 
                info_index = anno_hash_dict[record['chr']][
                    supple_of_record['hash_list'][
                        supple_of_record["index_list"][
                            supple_of_record["start_sort_list"][index]
                        ]
                    ]
                ]

                if info_index['is_N'] == 1:
                    record["num_ninteractgenes"] = record["num_ninteractgenes"] + 1
                record["num_genes"] = record["num_genes"] + 1

            if record["num_ninteractgenes"] > 0 :
                record["containsNinteractgenez"] = 1
            if record["num_genes"] > 0:
                record['containsgenes'] = 1

            record_collector.append(record)
    sys.stdout.write("\n")
    return (record_collector,set(sample_counter))

def get_stats(record_collector: List[Dict]) -> Dict:
    '''
    Calculate the stats of this roh analysis
    '''
    total_gene_number = int()
    n_interact_number = int()
    contain_othergene_length = int()
    contain_Ngene_blocks = int()
    contain_Ngene_length = int()
    gene_density_Ninteract = float()
    contain_othergene_blocks = int()
    contain_othergene_length = int()
    gene_density_nonNinteract = float()
    contain_none_blocks = int()
    contain_none_length = int()

    for item in record_collector:
        total_gene_number = total_gene_number + int(item['num_genes'])
        if item['containsNinteractgenez'] == 1:
            n_interact_number = \
                n_interact_number + int(item['num_ninteractgenes'])
            contain_Ngene_blocks = contain_Ngene_blocks + 1
            contain_Ngene_length = contain_Ngene_length + int(item['length'])
            gene_density_Ninteract = gene_density_Ninteract + \
                                    int(item['num_genes'])/int(item['length'])
        elif item['containsgenes'] == 1:
            contain_othergene_blocks = contain_othergene_blocks + 1
            contain_othergene_length = contain_othergene_length + int(item['length'])
            gene_density_nonNinteract = gene_density_nonNinteract + \
                                    int(item['num_genes'])/int(item['length'])
        else:
            contain_none_blocks = contain_none_blocks +1
            contain_none_length = contain_none_length + int(item['length'])


    avg_contain_N = contain_Ngene_length/contain_Ngene_blocks
    avg_contain_o = contain_othergene_length/contain_othergene_blocks
    avg_without = contain_none_length/contain_none_blocks

    density_N = gene_density_Ninteract/contain_Ngene_blocks
    density_None = gene_density_nonNinteract/contain_othergene_blocks

    stat_dict = {
        'total gene number' : total_gene_number,
        'contain othergene length' : contain_othergene_length,
        'contain othergene blocks' : contain_othergene_blocks,
        'contain Ngene length' : contain_Ngene_length,
        'contain Ngene blocks' : contain_Ngene_blocks,
        'contain none length' : contain_none_length,
        'contain none blocks' : contain_none_blocks,
        'Average length Ngene':avg_contain_N,
        'Average length OtherGene':avg_contain_o,
        'Average length None':avg_without,
        'gene density nonNinteract' : "{:E}".format(density_None),
        'gene density Ninteract' : "{:E}".format(density_N),
        'Proportion of N interact blocks with all blocks that have genes' : 
            "{:.3%}".format(
                contain_Ngene_blocks/(
                    contain_othergene_blocks + contain_Ngene_blocks
                )
            ),
        'Proportion of N interact genes':
            "{:.3%}".format(n_interact_number/total_gene_number)
        } 

    return stat_dict

def permutation(
        test_stat:float,
        permu_str:str,
        record_collector:List,
        sample_list:List,
        perm_times:int = 1000
    ) -> float:
    # The permutation part, adapted from the perl script
    results = list()

    for p in range(perm_times):
        sys.stdout.write("Processed: {}/{}\r".format(p+1,perm_times))
        sys.stdout.flush()

        # other_len = list() ######3
        # N_len = list()     #############
        length_blocks_with_interacting=0
        n_with_interacting = 0
        length_blocks_with_other=0
        n_with_other = 0
        length_blocks_without=0

        for sample in sample_list:

            counter=0
            record_list = sample_list[sample]

            permu_list = [x for x in permu_str]
            random.shuffle(permu_list)

            for r_index in record_list:
                switch = 0
                record = record_collector[r_index]
                if record["num_genes"] == 0:
                    length_blocks_without = length_blocks_without + \
                                                int(record["length"])
                    continue

                for i in range(record["num_genes"]):
                    if permu_list[counter] == "1":
                        switch = 1  # found N_interact gene in permutated list
                    counter = counter + 1

                if switch == 1:
                    length_blocks_with_interacting = \
                            length_blocks_with_interacting + int(record['length'])
                    n_with_interacting = n_with_interacting + 1
                    # N_len.append(record['length']) ###########
                else:
                    length_blocks_with_other = \
                            length_blocks_with_other + int(record['length'])
                    n_with_other = n_with_other + 1
                    # other_len.append(record['length']) ##########

        stat = math.log(length_blocks_with_interacting/n_with_interacting) - \
                math.log(length_blocks_with_other/n_with_other)
        results.append(stat)

    results.sort()

    sign, index = binarySearch(results, test_stat, 0, perm_times)
    if sign == -1:
        true_index = 0
    elif index == perm_times -1 and sign == -2:
        true_index = perm_times
    else:
        true_index = index + 1

    quick_p = 1 - true_index/perm_times
        
    sys.stdout.write("\n")
    # For debugging, print the permuation stats list
    sys.stdout.write("Permutation stats: {}\n".format(str(results)))
    return quick_p

def arg_parser() -> Dict:
    '''
    parse arguments
    '''
    config = {
        'input_path':'',
        'annotation_path':'',
        'roh_output_path':Path('./roh_density.txt'),
        'stat_output_path':Path('./roh_stats.txt'),
    }
    opts, _ = getopt.getopt(
                sys.argv[1:],
                "i:a:o:s:h",
                ["annotation=", "input=", "output=","stats=","help"]
            )
    
    if not opts:
        print(__doc__)
        sys.exit(0)
    for k,v in opts:
        if k in {"-h","--help"}:
            print(__doc__)
            sys.exit(0)
        if k in {"-i","--input"}:
            config["input_path"] = Path(v)
            continue
        if k in {"-o","--output"}:
            config["roh_output_path"] = Path(v)
            continue
        if k in {"-a","--annotation"}:
            config["annotation_path"] = Path(v)
            continue
        if k in {"-s","--stats"}:
            config["stat_output_path"] = Path(v)
            continue
        print("Unrecognized argument: {}".format(k))
        sys.exit(1)
    
    return config

def main(config:Dict) -> None:
    sys.stdout.write("Analysis starts...\n")
    annotation_path = config["annotation_path"]
    roh_path = config["input_path"]
    out_path = config["roh_output_path"]
    stat_path = config["stat_output_path"]

    anno_hash_dict,permu_str = get_annotation_hash(annotation_path)
    supple_dict = get_supple_dict(anno_hash_dict)
    sys.stdout.write("Counting roh blocks... \n")
    record_collector, sample_list = record_analyser(
                                        roh_path,supple_dict,anno_hash_dict
                                    )
    stat_dict = get_stats(record_collector)

    test_stat = \
        math.log(stat_dict['contain Ngene length']/stat_dict['contain Ngene blocks']) - \
        math.log(stat_dict['contain othergene length']/stat_dict['contain othergene blocks'])
    sys.stdout.write("Performing permutation tests...\n")

    sample_Ndict = Ndict()
    for i in range(len(record_collector)):
        sample_Ndict.judgein(record_collector[i]['sample'],i)

    quick_p = permutation(test_stat,permu_str,record_collector,sample_Ndict)


    sys.stdout.write("Writing density files...\n")
    out_list = ["chr","pos","length","containsgenes",
    "containsNinteractgenez","num_genes","num_ninteractgenes"]
    with out_path.open('w') as f:
        _ = f.write("{}\n".format('\t'.join(out_list)))
        for item in record_collector:
            _ = f.write("{}\n".format('\t'.join([str(item[x]) for x in out_list])))

    global stat_output_sequence
    sys.stdout.write("Writing stats files...\n")
    with stat_path.open('w') as f:
        for i in stat_output_sequence:
            _ =f.write("{}: {}\n".format(i,stat_dict[i]))
        _ = f.write("Test Stat: {}\n".format(test_stat))
        _ = f.write("Quick P: {}".format(quick_p))

    sys.stdout.write("Completed!\n")



if __name__ == "__main__":
    main(arg_parser())


####### For Debugging, Do not enable them!
# new_stat = Ndict()
# for item in record_collector:
#     key = item['chr'] +"_"+ item['start'] + "_" + item['end']
#     new_stat.judgein(key,item)

# cnt = 0

# for i in new_stat:
#     if len(new_stat[i]) > 1 and new_stat[i][0]['containsNinteractgenez'] >0:
#         cnt = cnt + 1
config = {
#    'input_path':Path(r"G:\Now_work\2021_05_18.ben_\07.ROH_permutation\Data\roh_brown.txt"),
    'input_path':Path(r"F:\Work\Now_work\2021_05_18.ben_\06.corrected\roh_source\aureus_new\roh_thibetana.txt"),
    'annotation_path':Path(r"F:\Work\Now_work\2021_05_18.ben_\00.annotation\rheMac10.ncbiRefSeq.CDS.marked.corrected.out"),
    'roh_output_path':Path('./roh_density.txt'),
    'stat_output_path':Path('./roh_stats.txt'),
}

