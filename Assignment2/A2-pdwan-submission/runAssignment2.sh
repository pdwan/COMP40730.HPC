
# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# GIT : https://github.com/pdwan/COMP40730.HPC.git
# DUE DATE : 30-June-2014 (extended to : July-2014)
# ASSIGNMENT : 2
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
Now=$(date +"%Y%m%d.%H%M%S")
logDir="logDir"
txtSuffix=".txt"
datSuffix=".dat"
stdLogFile="runAssignment2-${Now}.log"
compileUsingCblas="false"
compileUsingAtlas="false"
buildManual="false"
buildSolo="false"
initRandom="false"
initIncrement="false"
matrixEnabled="false"
defaultValuesRange="false"
let matrixSize=0
let maxMatrixSize=1000 
let threadSize=0
let maxThreadSize=100

# ##################################################################################

# function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <enter> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0) \ \n\t-1|--manual -2|--solo -d1|--atlas -d2|--cblas -r|--random -i|--increment -m|--matrix <n> -t|--thread <t> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| and then infintiy norm using pthreads \n" 
    $_ECHO -e "LOGS :\tCreated in current dir and moved to [ ${logDir} ] : \n\t<file>.txt : \tmatrix values for matrices |A| |B| & |C| \n\t<file>.dat :\ttiming data for each computation \n\t<file>.log : \tsummary of stdout. \n"
    $_ECHO -e "WHERE :\t-1|--manual \tCompile A2-pthreads-manual.c : straight-forward IJK and DGEMM computations only"
    $_ECHO -e "\t-2|--solo \tCompile A2-pthreads-solo.c   : straight-forward IJK and pThreads computations only, only valid on 'yeats.ucd.ie' "
    $_ECHO -e "\t\t\t'-1|--manual' and '-2|--solo' are mutually exclusive \n"
    $_ECHO -e "\t-d1|--atlas\tCompile  A2-pthreads-manual.c source files using dgemm atlas"
    $_ECHO -e "\t-d2|--cblas\tCompile  A2-pthreads-manual.c source files using dgemm cblas\n"
    $_ECHO -e "\t-r|--random \tInitialize |A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <column> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive \n"
    $_ECHO -e "\t-m|--matrix <n>\tMatrix dimension, if odd number +1 added or if invalid set to [ ${maxMatrixSize} ], thread count set to [ ${maxThreadSize} ] "
    $_ECHO -e "\t-t|--thread <t>\tnumber of threads, if invalid set to  [ ${maxThreadSize} ]  and matrix size set to [ ${maxMatrixSize} ]"
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and for <np> as follows :"
    $_ECHO -e "\t\t\tRange 1\t<matrixArray> : { 50 50 50 100 100 100 500 500 500 1000 1000 1000 }"
    $_ECHO -e "\t\t\t\t<threadArray> : { 10 10 10 10 10 10 20 20 20 20 20 20 20 20 }"
    $_ECHO -e "\t\t\tRange 2\t<matrixArray> : { 50 50 50 50 50 50 100 100 100 100 100 100 }"
    $_ECHO -e "\t\t\t\t<threadArray> : { 2 2 2 5 5 5 10 10 10 20 20 20 }"
    $_ECHO -e "\t\t\tRange 3\t<matrixArray> : { 50 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 }"
    $_ECHO -e "\t\t\t\t<threadArray> : { 10 10 20 20 30 30 40 40 50 50 60 60 70 70 80 80 90 90 100 }"
    $_ECHO -e "\t\t\t'-m|--matrix <n>' & '-t|--thread <t>' are mutually exclusive of  '-v|--values'.\n"
    $_ECHO -e "\t-?|-h|--help \tusage \n"
}

# function : error message and then usage
error() 
{
    err=$1
    case ${err} in
        1)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Unknown parameter '${1}'." >> ${stdLogFile}
            ;;
        2)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating directory '${1}'." >> ${stdLogFile}
            ;;
        3)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating file '${1}'."  >> ${stdLogFile}
            ;;     
        4)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Missing parameter for '${2}'."  >> ${stdLogFile}
            ;;     
        5)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${2}', Values entered are not valid or not a number."  >> ${stdLogFile}
            ;;     
        6)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${2}', Mutually exclusive switches." >> ${stdLogFile}
            ;;     
        7)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Compilation failed." >> ${stdLogFile}
            ;;     
        *)
            $_ECHO -e "ERROR :\tUnknown error."  >> ${stdLogFile}
            $_ECHO ${err}  >> ${stdLogFile}
            $_ECHO $*  >> ${stdLogFile}
            ;;
    esac
    $_ECHO -e
    usage
    exit ${err}
}

# function : build using pthreads only
compile_pthreads()
{
    localProgramToCompile=$1
    $_ECHO -e "pThreads :\t\tCompiling ${localProgramToCompile} \n"  >> ${stdLogFile}  
    cc -Wall -I/home/cs/khasanov/libs/CBLAS/src -o ${localProgramToCompile} ${localProgramToCompile}.c  /home/cs/khasanov/libs/cblas_LINUX.a  /usr/lib/libblas.a -lgfortran -pthread
}

# function : build applying dgemm cblas
compile_dgemm_atlas()
{
    localProgramToCompile=$1
    $_ECHO -e "ATLAS :\t\tCompiling ${localProgramToCompile}-atlas using atlas \n"  >> ${stdLogFile}
    gcc -Wall -o ${localProgramToCompile}-atlas ${localProgramToCompile}.c -I/home/cs/khasanov/libs/ATLAS/include/ -L/home/cs/khasanov/libs/ATLAS/lib/Linux_UNKNOWNSSE2_4/ -lcblas -latlas -lm -O3 
}

# function : build applying dgemm cblas
compile_dgemm_cblas()
{
    localProgramToCompile=$1
    $_ECHO -e "CBLAS :\t\tCompiling ${localProgramToCompile} using cblas \n"  >> ${stdLogFile}
    gcc -Wall -I/home/cs/khasanov/libs/CBLAS/src ${localProgramToCompile}.c -o ${localProgramToCompile}-cblas  /home/cs/khasanov/libs/cblas_LINUX.a  /usr/lib/libblas.a -lgfortran
}

# function : create directory, if it does not exist & validate creation
init_dir() 
{
    creationDir=$1
    if [ ! -d ${creationDir} ] || [ ! -e ${creationDir} ] ; then 
        mkdir ${creationDir}        
        $_ECHO -e "WARNING :\tCreating $creationDir"  >> ${stdLogFile}
        if [[ $? -ne 0 ]] ; then 
            error 2 $creationDir
        fi
    fi  
}

# function : create log files (.txt : matrix values, .dat : timing data for each computation & .log : stderr, stdout) to store values for data for each alogrithim computation
init_log_file() 
{
    localLogFile=$1
    if  [ -e ${localLogFile} ]  ; then
        $_ECHO -e "WARNING :\tFile backup : ${localLogFile} to ${localLogFile}.bup" >> ${localLogFile}
        mv "${localLogFile}" "${localLogFile}.bup"
    fi
}

# function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
    localCmd="$1"
    localOptions="$2"
    localFileMatrix="$3"
    localFileTime="$4"
    $_ECHO -e "RUNNING :\t${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}"  >> ${stdLogFile}
    # $_ECHO -e "DEBUG : RUNNING :\t${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}"
    ./${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}
}

# ##################################################################################

#clear

# Process parameters
if [[ $# -eq 0 ]]; then
    usage 
    exit
else 
    while [ "$1" != "" ]; do
        case ${1} in
            "-1" | "--manual")
                buildManual="true"
                ;;
            "-2" | "--solo")
                buildSolo="true"
                ;;        
		    "-d1" | "--atlas")
				compileUsingAtlas="true"
                ;;
		    "-d2" | "--cblas")
				compileUsingCblas="true"
                ;;
            "-i" | "--increment") 
                initIncrement="true"
                ;;
            "-r" | "--random")
                initRandom="true"
                ;;
            "-i" | "--increment") 
                initIncrement="true"
                ;;
            "-v" | "--values")
                defaultValuesRange="true"
                ;;                    
            "-m" | "--matrix")
                matrixEnabled="true"
                if [[ $2 =~ "^[0-9]+$" ]] ; then 
                    let matrixSize=$2 
                else
                    error 5 "${1} ${2}"
                fi 
                shift 
                ;;
            "-t" | "--thread")
                threadEnabled="true"
                if  [[ $2 =~ "^[ 0-9 ]+$" ]] ; then 
                    let threadSize=$2 
                else
                    error 5 "${1} ${2}"
                fi  
                shift 
                ;;                        
            "-?" | "-h" | "--help")
                usage
                exit
                ;;
            *)
                error 1 $1
                ;;
        esac
        shift
    done
fi

# process and validate parameter values for matrix sizes, if applicable
if [ "${defaultValuesRange}" == "true" ] && ( [ "${matrixEnabled}" == "true" ] || [ "${threadEnabled}" == "true" ]  ) ; then 
    error 6 "<-m> & <-t> and <-v>"
fi
if  [ "${defaultValuesRange}" == "true" ] &&  ( [ "${matrixEnabled}" == "false" ] || [ "${threadEnabled}" == "false" ]  ); then
#   Matrix & pThreads - range 1
    declare -a NXArray=( 50 50 50 50 50 50 100 100 100 100 100 100 )
    declare -a NPArray=( 2 2 2 5 5 5 10 10 10 20 20 20 )
#   Matrix & pThreads - range 2
#    declare -a NXArray=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
#    declare -a NPArray=( 10 10 10 10 10 10 20 20 20 20 20 20 20 20 )
#   Matrix & pThreads - range 3
#     declare -a NXArray=( 50, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 )
#     declare -a NPArray=( 10, 10, 20, 20, 30, 30, 40, 40, 50, 50, 60, 60, 70, 70, 80, 80, 90, 90, 100 )
    matrixEnabled="false"
    threadEnabled="false"
else
    if [ "${matrixEnabled}" == "true" ] || [ "${threadEnabled}" == "true" ]  ; then
#   Validate matrix size initialized    
		if [ "${matrixSize}" == "" ] ; then 
		    error 4 "matrix size"
		fi
#   Validate thread size initialized    
		if [ "${threadSize}" == "" ] ; then 
		    error 4 "thread size"
		fi
#   Validate matrix range		
		if [ ${matrixSize} -le 0 ] || [ ${matrixSize} -gt ${maxMatrixSize} ] ; then
		    $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is invalid, default values used for matrixSize : $maxMatrixSize" >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize        
	        matrixEnabled="true"
	    fi   
#   Validate thread range		
		if [ ${threadSize} -le 0 ] || [ ${threadSize} -gt ${maxthreadSize} ] ; then
		    $_ECHO -e "WARNING :\t$($_BASENAME $0): thread size <np> is invalid, default values used for threadSize : $maxThreadSize" >> ${stdLogFile}
	        let threadSize=$maxThreadSize        
            threadEnabled="true"
	    fi   
#   ensure both are enabled	    	     
    	if  [ "${matrixEnabled}" == "false" ] ; then 
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is now enabled and set to default of : $maxMatrixSize"  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize
	        matrixEnabled="true"
        fi 
        if  [ "${threadEnabled}" == "false" ] ; then 
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): thread size <np> is now enabled and set to default of : $maxThreadSize"  >> ${stdLogFile}
	        let threadSize=$maxThreadSize
            threadEnabled="true"
        fi 
#   ensure matrix is even in size   
        if [[ $(expr ${matrixSize} % 2 ) -eq 0 ]] ; then 
                let matrixSize=$matrixSize 
            else 
                let matrixSize=$(( matrixSize++ )) 
        fi       
		declare -a NXArray=( ${matrixSize} )
		declare -a NPArray=( ${threadSize} )
    fi
fi

# validate manual and solo --> pThreads - mutually exclusive
if  [ "${buildManual}" == "true" ] && [ "${buildSolo}" == "true" ] ; then 
        error 6 "-1|--manual and -2|--solo"
fi 

# validate atlas and cblas - mutually exclusive
if  [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
        error 6 "-d1|--atlas and -d2|--cblas"
fi 

# build up commands to run
algorithmOptions=""
if [ "${initRandom}" == "true" ] && [ "${initIncrement}" == "false" ] ; then 
    algorithmOptions="-r"
fi
if [ "${initRandom}" == "false" ] && [ "${initIncrement}" == "true" ] ; then 
    algorithmOptions="-i"
fi
if [ "${initRandom}" == "true" ] && [ "${initIncrement}" == "true" ] ; then 
    error 6 "<-i> & <-r>" 
fi

# execute algorithms

init_dir ${logDir}
init_log_file ${stdLogFile}

matrixFileRoot="${Now}-values"
dataFileRoot="${Now}-data"

# build Manual : straight-forward IJK and DGEMM
if  [ "${buildManual}" == "true" ]  ; then
    algorithmManual="A2-pthreads-manual" 
    matrixFileManual="${matrixFileRoot}-${algorithmManual}"
    dataFileManual="${dataFileRoot}-${algorithmManual}"
    dataFileManualTiming="${dataFileManual}${datSuffix}" # append to existing file for graphing : retain as separate as may way to create one per iteration
    init_log_file ${dataFileManualTiming}
    if [[ ${#NXArray[*]} -le 0 ]] ; then 
    	error 5 "matrix size"
    else        
    	executeOptions=""
        if  [ "${compileUsingCblas}" == "false" ] ; then 
            if [ "${compileUsingAtlas}" == "false" ] ; then
    	        algorithmManual="${algorithmManual}-cblas" # default
            fi 	        
        fi
        if  [ "${compileUsingCblas}" == "true" ] ; then
    	    compile_dgemm_cblas ${algorithmManual}
    	    algorithmManual="$algorithmManual-cblas"
        fi 
        if  [ "${compileUsingAtlas}" == "true" ] ; then
    	    compile_dgemm_atlas ${algorithmManual}
    	    algorithmManual="${algorithmManual}-atlas"
        fi 
    	for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
    	    matrixFileManualValues="${matrixFileManual}-$i${txtSuffix}" # different file for each run
    	    init_log_file ${matrixFileManualValues} 
    	    executeOptions="${algorithmOptions} ${NXArray[$i]} "
    	    # echo "DEBUG : algorithm_execute  ./${algorithmManual} ${executeOptions} ${matrixFileManualValues} ${dataFileManualTiming}"
    	    algorithm_execute "./${algorithmManual}" "${executeOptions}" "${matrixFileManualValues}" "${dataFileManualTiming}"
    	done
    fi
fi 

# build Solo : straight-forward IJK and pThreads
if  [ "${buildSolo}" == "true" ]  ; then
    algorithmSolo="A2-pthreads-solo" 
    matrixFileSolo="${matrixFileRoot}-${algorithmSolo}"
    dataFileSolo="${dataFileRoot}-${algorithmSolo}"
    dataFileSoloTiming="${dataFileSolo}${datSuffix}" # append to existing file for graphing : retain as separate as may way to create one per iteration
    init_log_file ${dataFileSoloTiming}
    if [[ ${#NXArray[*]} -le 0 ]] ; then 
    	error 5 "matrix size"
    else        
    	executeOptions=""
        if  [ "${compileUsingCblas}" == "false" ] ; then 
            if [ "${compileUsingAtlas}" == "false" ] ; then
#    	        compile_pthreads ${algorithmSolo}
    	        algorithmSolo="${algorithmSolo}" # default
            fi 	        
        fi
        if  [ "${compileUsingCblas}" == "true" ] ||   [ "${compileUsingAtlas}" == "true" ] ; then
  		    $_ECHO -e "ERROR :\t$($_BASENAME $0): cannot compile using DGEMM, existing. \n"    
  		    usage
        fi 
    	for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
    	    matrixFileSoloValues="${matrixFileSolo}-$i${txtSuffix}" # different file for each run
    	    init_log_file ${matrixFileSoloValues} 
    	    executeOptions="${algorithmOptions} ${NXArray[$i]} ${NPArray[$i]} "
    	    # echo "DEBUG : algorithm_execute  ./${algorithmSolo} ${executeOptions} ${matrixFileSoloValues} ${dataFileSoloTiming}"
    	    algorithm_execute "${algorithmSolo}" "${executeOptions}" "${matrixFileSoloValues}" "${dataFileSoloTiming}" "${valueNP}"
    	done
    fi
fi 

# move log files to <logDir>
mv -f *.dat ${logDir} 2> /dev/null
mv -f *.txt ${logDir} 2> /dev/null
mv -f *.log ${logDir} 2> /dev/null
mv -f *.bup ${logDir} 2> /dev/null

pause
exit 0
