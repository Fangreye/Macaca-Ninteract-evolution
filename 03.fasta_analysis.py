'''
Build input file for paml analysis

usage: 
    python fasta_analysis.py 
                    <-i | --input fasta_input_path> 
                    <-o | --outdir output_path>
                    <-b | --bed bed_file>
                    <-p | --prefix prefix>
                    [-h | --help]

arguments:
    -i | --input roh_input_path
        The path to store fasta files
    -o | --output output_path
        Where to write the output
    -b | --bed bed_file
        Bed Annotation, need to be 1-based coordinates
    -p | --prefix prefix
        The prefix added to the output
    -g | --genes
        File contains a set of genes to analysis
    -h | --help :
        To show this message
'''


import re
import sys
import getopt
import random

from pathlib import Path
from typing import Tuple, Dict, List

oxphos_match = 'Cbp6|TMEM177|HEL-S-271|CK820_G0026760|C7H14orf2|CI-MLRQ|DYT8|PIG47|ATP5MF|MC1DN22|MC4DN15|MTND4L|ECGF1|CI-ESSS|hATP1|AGGG|CK820_G0026028|LP3663|KIPP1184|MC5DN6|CK820_G0024747|UQPC|NADHDH2|CK820_G0014643|CI-AGGG|CI-B14.5a|MC1DN7|UQOR1|APLCC|C7orf44|UQCRQ|MTND5|CK820_G0040381|MC4DN4|COX10-AS|GP130|HSPC009|NDUFV2|CK820_G0039915|C3orf60|TMEM220|RIS1|CI-KFYI|MC1DN25|COX7A2|B17.2|UQCC3|CI-B22|NDUFB2|C4orf52|TTC19|ATP5G|CI-75k|UQCR5|PET100|CI-B16.6|NDUFAF4|COX8|TP|PLPM|BRP17|HT007|CTB-32O4.3|TYMP|MT-ND2|COX6A1|LSDMCA2|ATPSB|C7H7orf44|ATP5K|MT-ND5|UQOR22|TMEM126B|SMIM20|SCO1|FRTS5|MC1DN11|UQCR11|PCAG1|ISP|QCR7|ATP5A|ATP5S|ATP5SL|C6orf66|0710008D09Rik|P17.3|CK820_G0043666|ASHI|NDUFA8|VIII|COXVIb1|EGK_02698|MR-1|LYRM3|CI-B9|CI-ASHI|FLNMS|EGK_20031|Hs.6719|ACP|COXVIIAL|CI-42k|MT-ND3|SCOD1|UQCC|CF6|ATPIF1|MC3DN1|MC1DN19|EGK_19443|MOM2|COII|CK820_G0000612|h-BCS|NDUFV3|CUNH14orf2|NDUFA13|SCO2|MITRAC15|ND3|ND4|2P1|HIGD1A|NDUFAF6|COX7AL1|COA3|FP634|ATP5MPL|CI-AQDQ|SHY1|B17.2L|MTCO2L|B15|CBP3|NDUFAF2|mimitin|ATP5PF|bA792D24.4|MYP6|RCF1a|CR201_G0024819|MC4DN22|COX10AS|COX1|TIMMDC1|B22|QtsA-18546|CR201_G0031684|ATP12p|COX18|NDUFAF3|CI51KD|PKND11|CI|MWFE|BJS|CK820_G0029448|ND5|C8H8orf38|CR201_G0051441|ATPQ|PDSW|VIII-L|ATP5CL1|CK820_G0037941|COX4I1|CR201_G0035326|CK820_G0037171|MC1DN31|HLC-1|CK820_G0038953|NDUFS4|MZM1L|CDA016|C6orf203|HRPAP20|NDUFA10|MC4DN2|ATP5ME|MTND4|MC1DN9|CI-10k|MTCO2|M19|HMC08D05|CSRP2BP|QCR8|BCS|C6H6orf203|NDUFS5|COX11|CTB-32O4.2|COXII|MP68|SGDH|COXPD22|OMR|MC1DN23|ATP5PD|MC1DN4|CK820_G0020966|COX4-1|NDUFB6|UQCR6|C20orf44|MLRQ|bA22L21.1|CR201_G0053834|UQCRFS1|DAPIT|C6orf125|E3-3|ATP5F1D|CI-13kD-A|CI-B15|MC4DN7|ATP5G1|LRPPRC|NP17.3|ATP5J2|UQCC2|CK820_G0018229|NDUFB7|PTD|MC4DN21|CI-B17.2|CISGDH|C9orf123|FPD1|ND4L|MC1DN37|ATP12|C2H3orf60|CK820_G0043085|VA|NDUFB5|AQDQ|SCO1L|OSCP|CR201_G0039082|CK820_G0015730|ATP5AL2|MT-ND4|CK820_G0004038|CK820_G0047453|ND6|MC1DN32|NDUFB10|ATP5PO|MC4DN11|ATP5F1E|COXVIIa-L|NDUFS1|HIG1|MC1DN20|C8orf38|MC1DN17|COX6B|GRIM19|ATP5MG|COX2|COX19|MC4DN12|GRIM-19|ATP5D|ATPI|CMT4K|CI-B17|CI-75kD|CD14|NDUFA3|MR1|MITRAC12|CCDS41658.1|NDUFA7|CIB8|APT5H|RIP1|HSPC230|HSPC203|ATP5C1|CI-SGDH|ATP5MD|CYC1|KCTD14|MC3DN8|CI-42KD|MTRES1|MC1DN35|ATP5H|CEMCOX2|RNF113A|COX25|R1|NDUFA12L|MC1DN28|CI-18|COX6A|MTCYB|MC1DN29|NDUFAB1|CI13KDA|NDUFA1|CI15K|BCS1|CK820_G0045826|FKSG19|MC4DN6|LSFC|ATPO|CR201_G0028559|hCOX16|CI-B18|ATPMB|CI-MNLL|LYRM7|CK820_G0046631|C3H14orf2|bA119F7.2|CI-B14.5b|COX5B|MC1DN15|CI-13kA|MC1DN30|ATP11|BCS1-like|COX5A|BFZB|COXG|TMEM186|ATP5I|UQCR4|MC1DN12|C14H14orf2|PDC|C14orf112|UQCR|B14.5a|CI-MWFE|COX14|CK820_G0011622|ZNF183|HEL-S-123m|EGK_02244|CI-PGIV|DMAC1|NDUFA2|MT-ND4L|CMC1|FAM36A|EGK_13350|CTB-32O4.4|ATP5MC1|NDUFB3|GRACILE|NDUFC1|ATP5A1|COX11P|CR201_G0027796|UQCR7|QP-C|NDUFB1|CI-19KD|CI-18kDa|UNQ655|FASN2A|CK820_G0051292|COX7B|MC4DN8|CI-B12|CEMCOX4|C14orf2|MC4DN5|UQCC1|NDUFB4|CR201_G0012181|MC1DN18|PET117|MC1DN13|MT-ND6|EGK_09658|NDUFC2|CR201_G0040922|C12orf62|CGI65|MC4DN16|NDUFA4|COX20|ATP5J|B14.5b|MAMU-B18|H17|EGK_00460|MTND2|COX6C|ECSIT|CCNI1|UQBP|COX15|MMTN|MC5DN1|ATP5MK|MT-CYB|CYI|MT-CO2|MC1DN1|DMAC2|ATP5O|MC4DN1|MC1DN36|UQBC|TACO1|MITRAC7|6.8PL|COX17|ATP5|ATP5IF1|Grim19|CK820_G0023634|ATP5F1B|MC3DN10|TMEM70|MC5DN3|ND2|B16.6|ORM|COX-VA|CLONE-23970|MC1DN10|CK820_G0050806|CMTRID|C1H1orf31|ATP5JL|ATP5L|C1orf31|B17|ptr-mir-1307|MC4DN13|CK820_G0021909|B8|ATP5JG|CI-9k|bA6B20.2|COA1|ATP5E|CK820_G0051323|C11orf83|COX7AL|COX8-2|MTND6|COXFA4|VIIAL|COX16|COXIV|ATP5B|MC1DN24|NDUFV1|C5orf31|NDUFB11|TdRPase|CI-24k|CCDC44|CIB17.2|CK820_G0034816|DAP13|TAHCCP2|CGI-39|MC3DN2|MC3DN9|CI-75Kd|CGI-65|COX|h-BCS1|COX8A|dJ342K12.1|PRO1304|CR201_G0035861|SDAP|C19orf79|MC5DN5|MISTR1|CK820_G0041111|MLQ|COI|CK820_G0038292|ACAD9|C3orf1|MC3DN3|NDUFAF1|EGK_09856|CYTB|MC3DN7|TMEM261|ATP5F1A|PGIV|ATPM|QCR10|CK820_G0032202|NDUFB9|NDUFS6|MC3DN4|B12|MNLL|RISP|OXA1L2|CR201_G0044510|ATPIP|COXIV-1|ESSS-P17.3|ATPAF2|KEG83_p11|2010204O13Rik|CI-9KD|IP|COX6AL|PNKD|CIA30|HSPC125|BCS1L|ATP11p|MR-1S|NPD002|Np15|ACP1|MC4DN14|LRP130|MTND3|AGP|EGK_11651|CEMCOX1|COX6B1|MC3DN6|ATPAF1|Gliostatin|COX7C|PD-ECGF|CR201_G0027438|ATP5PB|B9|CI-15k|C6orf126|ESSS|COA6|ATPE|CK820_G0011089|SURF1|NDUFA12|COX4|HCVFTP2|CK820_G0041943|COX10|Zfp469|PRED31|My013|CI-51K|CK820_G0030520|MC4DN20|NDUFB8|FOXRED1|UQCRB|PNX|MRCAF1|FAM36B|COX8L|CCDC56|USMG5|MC5DN4|MNF1|F6|MC1DN5|ATP5F1|ATP5MJ|MC4DN19|MC4DN10|ATP5C|QPC|LOC722212|LOC716161|ATP5F1C|DMAC2L|ATP5S'
oxphos_match = re.compile(oxphos_match)
rep_match = 'POLRMT|TFAM|TFB1M|TFB2M'
rep_match = re.compile(rep_match)
ars_match = re.compile(r'ARS2')
mrp_match = re.compile(r'MRP|ICT1')

match_key = {
    'oxphos': oxphos_match,
    'rep' : rep_match,
    'ars' : ars_match,
    'mrp' : mrp_match
}

comple_dict = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C'
}

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

def read_fasta(infile:Path) -> Dict:
    '''
    Read in fasta file.
    '''
    patt = re.compile(r">(?P<chr>chr\d*):(?P<start>\d*)-(?P<end>\d*)")
    with infile.open() as f:
        temp_holder = dict()
        for line in f:
            if ">" in line:
                kw = patt.match(line).groups()
                temp_holder['_'.join(kw)] = list()
            else:
                temp_holder['_'.join(kw)].append(line.replace('\n',''))
    
    for item in temp_holder:
        temp_holder[item] =(item.split('_')[1] ,''.join(temp_holder[item]))

    return temp_holder

def read_bed(infile:Path) -> Dict:
    '''
    Read in bed file
    '''
    with infile.open() as f:
        temp_holder = dict()
        for line in f:
            kw = line.split()
            temp_holder['_'.join(kw[0:3])] = (kw[3], kw[5])
    return temp_holder

def check_length(info_dict:Dict) -> List:
    '''
    check if the seuqence length is correspondant to the annotation
    '''
    report_list = list()
    for key in info_dict:
        used = key.split('_')
        length = int(used[2]) - int(used[1]) + 1
        seq_length = len(info_dict[key][1])
        if seq_length != length:
            report_list.append(key)
        
    return report_list

def get_name(location:str, bed_dict:Dict):
    return bed_dict[location]

def get_sequence(sequence:str, strand:int=0) -> str:
    '''
    1: negative strand; 0: positive strand
    '''
    global comple_dict

    temp_list = list(sequence.upper().replace('*','N'))
    if strand == 1:
        temp_list.reverse()
        temp_list = [comple_dict[x] if x in comple_dict else x for x in temp_list ]
    
    return ''.join(temp_list).replace('*','N')

def get_sampled_gene(infile) -> List:
    with infile.open() as f:
        temp_holder = [x.replace('\n','') for x in f.readlines()]
    return temp_holder

def convert_legal(instring:str) -> str:
    return instring.replace("*","").replace("-","_")

def get_n_interact(instring: str) -> int:
    '''
    decide whether a gene is N-interact
    '''
    global match_key
    if match_key['oxphos'].fullmatch(instring):
        return 1
    elif match_key['rep'].fullmatch(instring):
        return 1
    elif match_key['ars'].search(instring):
        return 1
    elif match_key['mrp'].search(instring):
        return 1
    else:
        return 0


def write_outfile(out_dir,gene_name, sample_name,gene_sequence ,
                    prefix, is_interact) -> None:

    if is_interact == 1:
        cate = "Ninteract"
    elif is_interact == 0:
        cate = "Normal"

    with out_dir.joinpath(cate).joinpath(
                "{}{}.fasta".format(prefix, convert_legal(gene_name))
            ).open("a",encoding="UTF-8") as f:
        _ = f.write("{}     {}\n".format(sample_name,gene_sequence))

def arg_parser() -> Dict:
    '''
    parse arguments
    '''

    config = {
        "input":"",
        "prefix":"",
        "out_dir":"",
        "flag":0
    }

    opts, _ = getopt.getopt(
                sys.argv[1:],
                "i:o:p:b:g:h",
                ["input=","outdir=","preifx=","bed=","genes=","help"]
            )

    if not opts:
        print(__doc__)
        sys.exit(0)
    for k,v in opts:
        if k in {"-h", "--help"}:
            print(__doc__)
            sys.exit(0)
        if k in {"-i", "--input"}:
            config["input"] = Path(v)
            if not config["input"].exists():
                raise FileNotFoundError("No file inside!")
        if k in {"-o", "--outdir"}:
            config["out_dir"] = Path(v)
            if not config["out_dir"].exists():
                sys.stdout.write("Create output path: {}".format(v))
                config["out_dir"].mkdir()
            if not config["out_dir"].joinpath("Ninteract").exists():
                config["out_dir"].joinpath("Ninteract").mkdir()
            if not config["out_dir"].joinpath("Normal").exists():
                config["out_dir"].joinpath("Normal").mkdir()
        if k in {"-p", "--prefix"}:
            config["prefix"] = "{}_".format(v)
        if k in {"-b", "--bed"}:
            config["bed"] = Path(v)
        if k in {"-g", "--genes"}:
            config["genes"] = Path(v)
            config["flag"] = 1

    return config

def main(config:Dict):
    patt2 = re.compile(r"(?P<name>.*)\.fasta")

    bed_into = read_bed(config["bed"])

    for infile in config["input"].iterdir():
        if not patt2.match(infile.name):
            continue
    
        sample_name = patt2.match(infile.name).groupdict()["name"]
        fasta_info = read_fasta(infile)
        
        gene_seq_dict = Ndict()
        for item in fasta_info:
            gene_seq_dict.judgein(get_name(item, bed_into), fasta_info[item])

        gene_dict = Ndict()
        for i in bed_into:
            gene_dict.judgein(bed_into[i][0], i)

        whole_seq = dict()
        for item in gene_seq_dict:
            name, strand = item
            if strand == "+":
                whole_seq[name] = ''.join(
                                        [get_sequence(x[1]) for x in sorted(
                                            gene_seq_dict[item], key = lambda x:x[0]
                                            )
                                        ]
                                    )
            elif strand == "-":
                whole_seq[name] = ''.join(
                                        [get_sequence(x[1],1) for x in sorted(
                                            gene_seq_dict[item], 
                                            key = lambda x:x[0], reverse= True
                                            )
                                        ]
                                    )

        if config['flag'] == 0:
            for gene in whole_seq:
                is_iteract = get_n_interact(gene)
                write_outfile(config["out_dir"], gene, sample_name, whole_seq[gene],
                        config["prefix"], is_iteract)
        elif config['flag'] == 1:
            sampled_gene = get_sampled_gene(config["genes"])
            for gene in sampled_gene:
                is_iteract = get_n_interact(gene)
                write_outfile(config["out_dir"], gene, sample_name, whole_seq[gene],
                        config["prefix"], is_iteract)

if __name__ == "__main__":
    main(arg_parser())

###
# something = list()
# for gene in whole_seq:
#     if len(whole_seq[gene]) % 3 != 0:
#         something.append(gene)
# f1 = Path(r"F:\Work\Now_work\2021_05_18.ben_\10.dNdS\08.analysis_fasta\v10_Normal_qualified.txt").open('w')
# f2 = Path(r"F:\Work\Now_work\2021_05_18.ben_\10.dNdS\08.analysis_fasta\v10_Ninteract_qualified.txt").open('w')

# ninteract = list()
# normal = list()
# for gene in whole_seq:
#     flag = get_n_interact(gene)
#     if flag == 1:
#         ninteract.append(gene)
#         _ = f2.write("{}\n".format(gene))
#     else:
#         normal.append(gene)
#         _ = f1.write("{}\n".format(gene))
# f1.close()
# f2.close()

# import random

# with Path(r"F:\Work\Now_work\2021_05_18.ben_\10.dNdS\08.analysis_fasta\sample_Aureus.txt").open('w') as f:
#     for gene in ninteract:
#         _ = f.write("{}\n".format(gene))
#     for gene in random.sample(normal,3000):
#         _ = f.write("{}\n".format(gene))
