#!/bin/sh
#SBATCH --job-name=SRR1927139_MMUL.IN-28474_mito
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --time=4:00:00
#SBATCH --mem=65gb
#SBATCH --output=popgenWindows.%J.out
#SBATCH --error=popgenWindows.%J.err
#SBATCH --account=def-ben

perl /home/zhu46/softwares/NOVOPlasty-master/NOVOPlasty4.3.1.pl -c /home/zhu46/scratch/maca_800_config/SRR1927139.txt