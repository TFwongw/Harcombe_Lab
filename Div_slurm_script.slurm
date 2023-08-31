#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=16
#SBATCH --mem=100g
#SBATCH --tmp=10g
#SBATCH --mail-type=NONE


module load python3
source venv/bin/activate
unset PYTHONPATH
unset PYTHONHOME
unset PYTHONSTARTUP

cd Run_Comets/
python Div_COMETS.py
