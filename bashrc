alias tma='tmux a -t'
alias tml='tmux ls'
alias tmk='tmux kill-session -t'
alias tmn='tmux new -s'
#source phenix
#source /usr/local/phenix-1.12-2829/phenix_env.sh

#make ls command show type specific colors
#alias ls='ls -G'
alias ls='ls --color=auto'
alias chimera='/Applications/Chimera.app/Contents/MacOS/chimera'
export PHENIX_TRUST_OTHER_ENV=1

#launch pymol through script
alias pymol='~/scripts/pymol_launcher.py'
alias chimera='~/scripts/chimera_launcher.py'

#add script directory to pythonpath
export PYTHONPATH=$PYTHONPATH:/Users/brandonfrenz/scripts/:/Users/brandonfrenz/boruta_py:/Users/brandonfrenz/PyRosetta4.Release.python36.mac.release-221

# The next line updates PATH for the Google Cloud SDK.
if [ -f '/Users/brandonfrenz/google-cloud-sdk/path.bash.inc' ]; then source '/Users/brandonfrenz/google-cloud-sdk/path.bash.inc'; fi

# The next line enables shell command completion for gcloud.
if [ -f '/Users/brandonfrenz/google-cloud-sdk/completion.bash.inc' ]; then source '/Users/brandonfrenz/google-cloud-sdk/completion.bash.inc'; fi
alias ARGO=/Users/brandonfrenz/argo-rest/
alias al='argo list | head'
alias ad='argo delete'
alias ag='argo get'
alias kl='kubectl logs'
alias kubectlcp='kubectl cp --container=main'
#[ -f /usr/local/etc/bash_completion ] && . /usr/local/etc/bash_completion
#. "/Users/brandonfrenz/.conch/__init__.sh"
export ARTI_NAME=brandon
export ARTI_PASS=u8#19o5N
export OE_LICENSE='/Users/brandonfrenz/oe_license.txt'
export ROSETTA_MAIN='/Users/brandonfrenz/rosetta-main/'

LS_COLORS=$LS_COLORS:'di=1;34:ow=1;34:'; export LS_COLORS
export PS1="\W \$"
