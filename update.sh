rsync -avp ../nana/gdm/scripts/ ./scripts
rsync -avp ../nana/code/ ./code/
rsync -avp ../nana/gdm/data/*.csv ./data
rsync -avp ../nana/README.md .
rm -r scripts/scratch/
