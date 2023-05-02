# Update the remote host with files and directories from the local host ####
rsync -chavzP -e 'ssh -p 732' --exclude-from .gitignore . gerversk@files.cgrb.oregonstate.edu:spatafora/projects/taxonomy/
