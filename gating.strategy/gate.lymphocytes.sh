
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
    echo " -h : Prints this message."
    exit 1
}



alias qsub="qsub -cwd -b y"

function make_prior() {
    K=$1
    downsample=$2
    INFILE=$3
    rm --force $basedir/Err/prior.err
    rm --force $basedir/Out/prior.out
    chmod ug+rx $GATE_LYMPH
    echo qsub -N make.prior$K -cwd -b y -u nikolas -e $basedir/Err/prior.err -o $basedir/Out/prior.out $GATE_LYMPH --in.file $INFILE --out.dir $basedir/RData --plot.dir $basedir/Plot --post.threshold .8 --down.sample $downsample -K $K --channels $channels
    qsub -N make.prior$K -cwd -b y -u nikolas -e $basedir/Err/prior.err -o $basedir/Out/prior.out $GATE_LYMPH --in.file $INFILE --out.dir $basedir/RData --plot.dir $basedir/Plot --post.threshold .8 --down.sample $downsample -K $K --channels $channels
}


function gate_lymph() {
    K=$1
    downsample=$2
    INFILE=$3
    rm --force $basedir/Out/*
    rm --force $basedir/Err/*
    # wait until prior error file exists
    prior_err=${basedir}/Err/prior.err
    while [ ! -f $prior_err ]
    do
      sleep 60
    done
    prior=${basedir}/RData/prior.RData
    if [ ! -e $prior ]
    then
        echo $prior still does not exist!
        echo check $prior_err
    fi
    for x in `cat $INFILE`
        do
        chmod ug+rx 
        y=`basename $x`
        echo qsub -hold_jid make.prior$K -N $y -cwd -b y -u $USER -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out $GATE_LYMPH --in.file $x --out.dir ${basedir}/RData --plot.dir ${basedir}/Plot --prior $prior --post.threshold .5 --down.sample $downsample -K $K --channels $channels
        qsub -hold_jid make.prior$K -N $y -cwd -b y -u $USER -e ${basedir}/Err/${y%.fcs}.err -o ${basedir}/Out/${y%.fcs}.out $GATE_LYMPH --in.file $x --out.dir ${basedir}/RData --plot.dir ${basedir}/Plot --prior $prior --post.threshold .5 --down.sample $downsample -K $K --channels $channels
        sleep 1
    done
}



function build() {
    K=$1
    outfile=$basedir/Out/lymph.count$K.csv
    echo "fcsFile,lymph.count$K" > $outfile
    cat $basedir/Out/*.out | grep lymph.count | cut -d'/' -f9 | cut -d, -f1,3 | sort >> $outfile
    echo $outfile
    Rscript ~nikolas/IL2RA/bin/compare.lymph.count.R --file1 ~/IL2RA/CellPhenotypes/manual.csv --file2 $outfile --out.file $basedir/manual-agreement.pdf
}


function build2() {
    K=$1
    outfile=$basedir/Out/lymph.count.mix$K.csv
    echo "fcsFile,lymph.count$K" > $outfile
    for f in $basedir/RData/*.RData
    do
        echo $f
        x=`echo "load('${f}'); Sys.sleep(1); count <- round(dim(d)[1]*res@w[which.max(res@mu[,3])]); print(count);" | R --no-save --no-restore --silent | grep '[1]'`
        x=`echo $x | awk '{print($(NF))}'`
        echo `basename "${f%.RData}"`,$x >> $outfile
    done
    echo $outfile
    #need to delete last line since it contains some artefact
    head -n -1 $outfile > $outfile
    Rscript ~nikolas/IL2RA/bin/compare.lymph.count.R --file1 ~/IL2RA/CellPhenotypes/manual.csv --file2 $outfile --out.file $basedir/manual-agreement-mix.pdf
}




function main() {
    channels=FSCA,SSCA,CD4
    basedir=
    K=0
    gate=
    build=0
    #parse the args
    while getopts "c:b:K:g:h" optionName
    do
        case "$optionName" in
        c) channels=$OPTARG;;
        b) basedir=$OPTARG;;
        K) K=$OPTARG;;
        g) gate=$OPTARG;;
        ?) usage 0;;
        esac
    done
    if [[ "$K" -le 0 || "$basedir" == "" || "$gate" == "" ]]
    then
        usage
    fi
    basedir=${basedir}/Lymphocytes$K
    echo $basedir
    mkdir -p $basedir
    mkdir -p $basedir/RData
    mkdir -p $basedir/Plot
    mkdir -p $basedir/Out
    mkdir -p $basedir/Err
    prior=${basedir}/RData/prior.RData
    if [ ! -e $prior ]
    then
        echo $prior does not exist!
        echo will create prior
        n=`wc -l $gate | cut -f1 -d ' '`
        downsample=`echo "scale=0; 10000/$n" | bc -l`
        make_prior $K $downsample $gate
        gate_lymph $K .1 $gate
    elif [ "$gate" != 0 ]
    then
        echo $prior exists
        gate_lymph $K .1 $gate
    elif [ "$build" = 1 ]
    then
        build $K
    elif [ "$build" = 2 ]
    then
        build2 $K
    fi
}

### MAIN
main $*
