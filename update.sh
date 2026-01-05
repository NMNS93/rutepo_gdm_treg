rm -r ./scripts ./code ./data
rsync -avp ../nana/gdm/scripts/ ./scripts
rsync -avp ../nana/code/ ./code/
rsync -avp ../nana/gdm/data/*.csv ./data
