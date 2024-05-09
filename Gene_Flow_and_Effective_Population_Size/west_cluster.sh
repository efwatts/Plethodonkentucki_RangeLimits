#!/bin/bash
#SBATCH --account=PHS0328
#SBATCH --time=48:00:00
#SBATCH --ntasks=40




# change to denovo folder 
cd /fs/scratch/PHS0328/kentucki_EFW/migrate


# run migrate
/users/PHS0328/efwatts/migrate-5.0.4/src/migrate-n parmfile_west
