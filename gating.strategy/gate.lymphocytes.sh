
#strict mode for shell, makes sure all commands exit with 0 and that all
set -e
set -u
# for debugging
#set -x 

#to turn on the extended pattern matching features
shopt -s extglob

GATE_LYMPH=gate.lymphocytes.flowClust.R
GATE_LYMPH=`dirname $0`/$GATE_LYMPH


function usage() {
    echo "syntax: $0"
    echo " -b : <basedir> Basedir under which to write Lymphocytes/Out, Lymphocytes/Err, Lymphocytes/Plot and Lymphocytes/RData"
    echo " -c : <channels> Channels comma separated. Default is FSCA,SSCA,CD4"
    echo " -K : <K> Number of clusters. Mandatory argument must be larger than zero. Default is 5."
    #echo " -d : Downsample."
    #echo "One of these 3 options:"
    #echo "      -p : Make prior from list of files in csv file"
    #echo "      -g : Gate list of files specified in csv file"
    #echo "      -l : Build lymph count file"
    echo " -g : <file of FCS filenames> Looks for prior under <basedir>/RData/prior.RData. If no prior is found, first compute prior from FCS files. Then gate lymphocytes."
    #echo " -p : <CD4 Lymph CD4 MFI> The expected MFI of the CD4+ lymph cluster."
    #echo " -p : <population number> Number of population in prior."
    echo " -h : Prints this message."
    exit 1
}



alias qsub="qsub -cwd -b y"

function make_prior() {
    K=$1
    downsample=$2
    INFILE=$3
    p=$4
    rm --force $basedir/Err/prior.err
    rm --force $basedir/Out/prior.out
    rm --force $basedir/Plot/prior.fcs.png
    chmod ug+rx $GATE_LYMPH
    echo qsub -N make.prior$K -cwd -b y -u $USER -e $basedir/Err/prior.err -o $basedir/Out/prior.out $GATE_LYMPH --in.file $INFILE --out.dir $basedir/RData --plot.dir $basedir/Plot --post.threshold .8 --down.sample $downsample -K $K --channels $channels --cd4.mfi $p
    qsub -N make.prior$K -cwd -b y -u $USER -e $basedir/Err/prior.err -o $basedir/Out/prior.out $GATE_LYMPH --in.file $INFILE --out.dir $basedir/RData --plot.dir $basedir/Plot --post.threshold .8 --down.sample $downsample -K $K --channels $channels --cd4.mfi $p
}


function gate_lymph() {
    K=$1
    downsample=$2
    INFILE=$3
    p=$4
    rm --force $basedir/Out/*
    rm --force $basedir/Err/*
    # wait until prior plot file exists
    prior_err=${basedir}/Err/prior.err
    #prior=${basedir}/RData/prior.RData
    if [ ! -e $prior ]
    then
        echo $prior still does not exist!
        cat $prior_err
        exit 1
    fi
    for x in `cat $INFILE`
        do
        chmod ug+rx $GATE_LYMPH
        y=`basename $x`
        echo qsub -hold_jid make.prior$K -N J$y -cwd -b y -u $USER -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out $GATE_LYMPH --in.file $x --out.dir ${basedir}/RData --plot.dir ${basedir}/Plot --prior $prior --post.threshold .5 --down.sample $downsample -K $K --channels $channels --cd4.mfi $p
        qsub -hold_jid make.prior$K -N J$y -cwd -b y -u $USER -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out $GATE_LYMPH --in.file $x --out.dir ${basedir}/RData --plot.dir ${basedir}/Plot --prior $prior --post.threshold .5 --down.sample $downsample -K $K --channels $channels --cd4.mfi $p
        sleep 1
    done
}



# compares to manually gated phenotype
function compare() {
    K=$1
    outfile=$basedir/Out/lymph.count$K.csv
    echo "fcsFile,lymph.count$K" > $outfile
    cat $basedir/Out/*.out | grep lymph.count | cut -d'/' -f9 | cut -d, -f1,3 | sort >> $outfile
    echo $outfile
    echo $basedir/manual-agreement.pdf
    echo $basedir/manual-agreement.csv
    Rscript $HOME/Projects/IL2RA/bin/compare.lymph.count.R --file1 $HOME/Projects/IL2RA/CellPhenotypes/manual.csv --file2 $outfile --out.file $basedir/manual-agreement.pdf > $basedir/manual-agreement.csv
}


# args are set in README of each project
function main() {
    K=
    channels=FSCA,SSCA,CD4
    basedir=
    gate=
    p=
    downsample=.95
    ## IL2 stim
    #basedir=~/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3
    #gate=${basedir}/PSTAT5-CD25-CD45RA-CD4-FOXP3.csv
    #gate=${basedir}/head20-PSTAT5-CD25-CD45RA-CD4-FOXP3.csv
    #gate=${basedir}/failed.csv
    #p=2.5
    ## IL2RA
    #basedir=~/dunwich/Projects/IL2RA/
    #gate=${basedir}/prior-FCS.csv
    #gate=${basedir}/all-FCS.csv
    #p=3
    ## DILT1D
    #basedir=~/dunwich/Projects/DILT1D/donor7/EFF/
    #gate=${basedir}/EFF-fcs.csv
    #p=3
    #parse the args
    while getopts "c:b:K:g:p:h" optionName
    do
        case "$optionName" in
        c) channels=$OPTARG;;
        b) basedir=$OPTARG;;
        K) K=$OPTARG;;
        g) gate=$OPTARG;;
        p) p=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [[ "$K" -le 0 || "$basedir" == "" ]]
    then
        usage
    fi
    basedir=${basedir}/Lymphocytes$K
    echo basedir: $basedir
    echo infile: $gate
    mkdir -p $basedir
    mkdir -p $basedir/RData
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err
    prior=${basedir}/RData/prior.RData
    #check that all files in $gate are valid
    if [[ "$gate" != "" ]]
    then
        if [[ ! -e $gate ]]
        then
            echo $gate file does not exist!
            exit 1
        fi
        while read -r f
        do
            if [[ ! -e $(eval echo $f) ]]
            then
                echo $f in $gate does not exist!
                exit 1
            fi
        done < $gate
    fi
    if [ ! -e $prior ]
    then
        echo $prior does not exist!
        echo will create prior
        n=`wc -l $gate | cut -f1 -d ' '`
        #n=`cat $gate | wc -l`
        echo $n
        downsample=`echo "scale=0; 10000/$n" | bc -l`
        make_prior $K $downsample $gate $p
        prior_plot=${basedir}/Plot/prior.fcs.png
        while [ ! -f $prior_plot ]
        do
          sleep 60
        done
        gate_lymph $K .1 $gate $p
    elif [ "$gate" != "" ]
    then
        echo $prior exists
        gate_lymph $K $downsample $gate $p
    else
        compare $K
    fi
}

### MAIN
main $*
