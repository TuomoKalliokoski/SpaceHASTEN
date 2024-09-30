__conda_setup="$('/data/programs/oce/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/programs/oce/etc/profile.d/conda.sh" ]; then
        . "/data/programs/oce/etc/profile.d/conda.sh"
    else
        export PATH="/data/programs/oce/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

