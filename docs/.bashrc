# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
# Make emacs the default editor
EDITOR=/usr/bin/emacs
export EDITOR
# Set the shell editor to emacs
set -o emacs
# Make control-d not end a shell
set -o ignoreeof

#=======Spack=======
#export SPACK_ROOT=/data/ctm/em440/legcfs01/Libraries/spack
#. $SPACK_ROOT/share/spack/setup-env.sh

#==============================================================================
# %%%%% Settings for colorization  %%%%%
#==============================================================================
export GREP_OPTIONS=--color=auto
export GREP_COLOR=32
export LS_COLORS='no=00:fi=00:di=01;33:ln=01;36:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:su=37;41:sg=30;43:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;37:*.tgz=01;37:*.arj=01;37:*.taz=01;37:*.lzh=01;37:*.zip=01;37:*.z=01;37:*.Z=01;37:*.gz=01;37:*.bz2=01;37:*.deb=01;37:*.rpm=01;37:*.jar=01;37:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.avi=01;35:*.fli=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.flac=01;35:*.mp3=01;35:*.mpc=01;35:*.ogg=01;35:*.wav=01;35:'

#==============================================================================
# %%%%% Environment settings %%%%%%
#==============================================================================
export TERM=xterm                        # Default terminal
export EDITOR=emacs                      # Default editor (emacs or vim)
export VISUAL=emacs                      # Default editor (emacs or vim)
export GIT_EDITOR=emacs                  # Default Git editor (emacs or vim)
export MAIL=/usr/spool/mail/$USER        # Default mail program
mail=$MAIL

#==============================================================================
# %%%%% Shortcut aliases (add your own to .my_personal_settings) %%%%%
#==============================================================================
# Many of these aliases are specific to the Harvard Odyssey servers, so may 
# not work on ALICE/SPECTRE. 

# %%%%% Alias to pull env directory updates %%%%%
alias  getenv="cd ~/env; git pull origin master"

# %%%%% Aliases to load different compiler versions %%%%%
alias  load_15="export LOAD_COMPILER=ifort15; sb; unset LOAD_COMPILER"
alias  load_13="export LOAD_COMPILER=ifort13; sb; unset LOAD_COMPILER"
alias  load_11="export LOAD_COMPILER=ifort11; sb; unset LOAD_COMPILER"
alias  load_pg="export LOAD_COMPILER=pgi; sb; unset LOAD_COMPILER"

# %%%%% Aliases for general Unix commands %%%%
alias  disk="du -h -s -c"
alias  g="grep -in"
alias  gf="gifview"
alias  gt="grep -in --text"
alias  gf="gifview -a"
alias  gvs="gv --orientation=seascape"
alias  m="less"
alias  me="xterm &"
alias  pan="$HOME/bin/panoply.sh &"
alias  proc="ps -ef | grep $USER | sort"
alias  pu="rm *~"
alias  sb=". ~/.bash_profile"
alias  ssh="ssh -X -A"
alias  gsh="ssh -xT"
alias  tf="tail --follow"
alias  zap="kill -9"

# %%%%% Aliases for directory listing %%%%%
alias  ls="ls -CF --time-style=long-iso --color=auto"
alias  l1="ls -1"
alias  ll="ls -l"
alias  llt="ls -lt"
alias  lltm="ls -lt | less"
alias  la="ls -a"
alias  lla="ls -la"
alias  llh="ls -lh"

# %%%%% Aliases for Git and GitHub %%%%%
alias  gui="git gui &"
alias  gk="gitk &"
alias  gka="gitk --all &"
alias  gl="git log"
alias  glo="git log --oneline"
alias  glp="git log --pretty=format:'%h : %s' --topo-order --graph"

# %%%%% Aliases for IDL and GAMAP %%%%%
alias  I="cd $HOME/IDL"
alias  IG="cd $HOME/IDL/gamap2"
alias  IS="cd $HOME/IDL/tests"

# %%%%% Aliases for MET data %%%%%
alias  XC="cd $dataDir/GEOS_0.5x0.666_CH"
alias  XC5="cd $dataDir/GEOS_0.5x0.666_CH.d/GEOS_5"
alias  XS="cd $dataDir/GEOS_0.25x0.3125_NA/"
alias  XNFP="cd $dataDir/GEOS_0.25x0.3125_NA.d/GEOS_FP/"
alias  XN="cd $dataDir/GEOS_0.5x0.666_NA"
alias  XN5="cd $dataDir/GEOS_0.5x0.666_NA.d/GEOS_5"
alias  X1="cd $dataDir/GEOS_1x1"
alias  XNN="cd $dataDir/GEOS_NATIVE"
alias  X2="cd $dataDir/GEOS_2x2.5"
alias  X24="cd $dataDir/GEOS_2x2.5.d/GEOS_4_v4"
alias  X25="cd $dataDir/GEOS_2x2.5.d/GEOS_5"
alias  X2FP="cd $dataDir/GEOS_2x2.5.d/GEOS_FP"
alias  X2M="cd $dataDir/GEOS_2x2.5.d/MERRA"
alias  X4="cd $dataDir/GEOS_4x5"
alias  X44="cd $dataDir/GEOS_4x5/GEOS_4_v4"
alias  X45="cd $dataDir/GEOS_4x5/GEOS_5"
alias  X4FP="cd $dataDir/GEOS_4x5/GEOS_FP"
alias  X4M="cd $dataDir/GEOS_4x5/MERRA"

#==============================================================================
# %%%%% Source a file with your personal aliases and settings %%%%%
#==============================================================================

if [[ -f $HOME/.my_personal_settings ]] ; then
 . $HOME/.my_personal_settings
fi

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/e/em440/miniconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/e/em440/miniconda/etc/profile.d/conda.sh" ]; then
        . "/home/e/em440/miniconda/etc/profile.d/conda.sh"
    else
        export PATH="/home/e/em440/miniconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

export PATH=~/miniconda/bin:$PATH

conda config --set auto_activate_base false
