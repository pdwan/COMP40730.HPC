#!/bin/bash 

# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# GIT : https://github.com/pdwan/COMP40730.HPC.git
# DUE DATE : 30-June-2014 (extended to : 13-July-2014)
# ASSIGNMENT : 4
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
Now=$(date +"%Y%m%d.%H%M%S")
logDir="logDir"
logPrefix="pdwan-"
txtSuffix=".txt"
DatSuffix=".dat"
pngSuffix=".png"
stdLogFile="runAssignment3-${Now}.log"
compileUsingCblas="false"
initRandom="false"
initIncrement="false"
matrixEnabled="false"
threadEnabled="false"
defaultMatrixRange="false"
let matrixSize=0
let threadSize=0
let maxMatrixSize=1000 
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
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0) \ \n\t -d2|--cblas -r|--random -i|--increment -m|--matrix <n> -t|--thread <t> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| and then infintiy norm using mpi \n" 
    $_ECHO -e "LOGS :\tCreated in current dir and moved to <${logDir}> : \n\t<file>.txt : \tmatrix values for matrices |A| |B| & |C|, \n\t<file>.dat :\ttiming data for each computation \n\t<file>.log : \tsummary of stdout. \n"
    $_ECHO -e "WHERE :\t-d2|--cblas\tCompile .c source files using dgemm cblas & openmp"
    $_ECHO -e "\t-r|--random \tInitialize |A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <column> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive \n"
    $_ECHO -e "\t-m|--matrix <n>\tMatrix dimension, if odd number +1 added or if invalid set to [ ${maxMatrixSize} ], thread count set to [ ${maxThreadSize} ] "
    $_ECHO -e "\t-t|--thread <t>\tnumber of threads, if invalid set to  [ ${maxThreadSize} ]  and matrix size set to [ ${maxMatrixSize} ]"
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t<matrixArray> :\t{ 50, 50, 50, 100, 100, 100, 500, 500, 500, 1000, 1000, 1000 } \n\t\t\t<threadArray> :\t{ 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20 }"
    $_ECHO -e "\t\t\t'-m|--matrix <n>' | -t|--thread <t>' and '-v|--values' mutually exclusive.\n"
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
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Missing parameter for '${1}'."  >> ${stdLogFile}
            ;;     
        5)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Values entered are not valid or not a number."  >> ${stdLogFile}
            ;;     
        6)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Mutually exclusive switches." >> ${stdLogFile}
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

# function : build applying dgemm cblas
compile_dgemm_cblas()
{
    localProgramToCompile=$1
    $_ECHO -e "CBLAS :\t\tCompiling ${localProgramToCompile} using cblas \n"  >> ${stdLogFile}
    gcc -Wall -I/home/cs/khasanov/libs/CBLAS/src ${localProgramToCompile}.c -o ${localProgramToCompile}  /home/cs/khasanov/libs/cblas_LINUX.a  /usr/lib/libblas.a -lgfortran -fopenmp
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
    $_ECHO -e "# LOG :\t${localLogFile} \n\tCreated on ${Now} by ${USER}." >> ${localLogFile}
}

# function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
    localCmd="$1"
    localOptions="$2"
    localFileMatrix="$3"
    localFileTime="$4"
    $_ECHO -e "RUNNING :\t${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}"  >> ${stdLogFile}
    ${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}
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
                defaultMatrixRange="true"
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
                if [[ $2 =~ "^[0-9]+$" ]] ; then 
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
if [ "${defaultMatrixRange}" == "true" ] && [ "${matrixEnabled}" == "true" ] ; then 
    error 6 "<-m> & <-v>"
fi
if [ "${defaultMatrixRange}" == "true" ] && [ "${threadEnabled}" == "true" ] ; then 
    error 6 "<-t> & <-v>"
fi
if  [ "${defaultMatrixRange}" == "true" ] && [ "${matrixEnabled}" == "false" ]  && [ "${threadEnabled}" == "false" ] ; then 
    declare -a NXArray=( 50 50 50 100 100 100 200 200 300 300 300 )
  	declare -a threadArray=( 10 10 10 10 10 10 20 20 20 20 20 20  )
    matrixEnabled="false"
    threadEnabled="false"
else
    if [ "${matrixEnabled}" == "true" ] || [ "${threadEnabled}" == "true" ]  ; then
		if [ "${matrixSize}" == "" ] ; then 
		    error 4 "matrix size"
		fi
		if [ "${threadSize}" == "" ] ; then 
		    error 4 "thread size"
		fi
		if [ ${matrixSize} -le 0 ] || [ ${matrixSize} -gt ${maxMatrixSize} ] ; then
		    $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is invalid, default values used for matrixSize : $maxMatrixSize and threadSize :$maxThreadSize  "  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize        
	        let threadSize=$maxThreadSize        
	        matrixEnabled="true"
            threadEnabled="true"
	    fi       
		if [ ${threadSize} -le 0 ] || [ ${threadSize} -gt ${maxThreadSize} ] ; then
		    $_ECHO -e "WARNING :\t$($_BASENAME $0): thread size <nt> is invalid, default values used for matrixSize : $maxMatrixSize and threadSize :$maxThreadSize  "  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize        
	        let threadSize=$maxThreadSize        
	        matrixEnabled="true"
            threadEnabled="true"
	    fi       
        if [[ $(expr ${matrixSize} % 2 ) -eq 0 ]] ; then 
                let matrixSize=$matrixSize 
            else 
                let matrixSize=$(( matrixSize++ )) 
        fi  
        if [[ ! ( $(expr ${matrixSize}  % ${threadSize} ) -eq 0 ) ]] ; then 
		    $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size is not multiple of thread size, default values used for matrixSize : $maxMatrixSize and threadSize :$maxThreadSize  "  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize        
	        let threadSize=$maxThreadSize        
	        matrixEnabled="true"
            threadEnabled="true"
        fi
		declare -a NXArray=( ${matrixSize} )
		declare -a threadArray=( ${threadSize} )
    fi
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

matrixFileRoot="${logPrefix}${Now}-values"
dataFileRoot="${logPrefix}${Now}-data"
algorithmName="A4-mpi-1D" 
matrixFileName="${matrixFileRoot}-${algorithmName}"
dataFileName="${dataFileRoot}-${algorithmName}"
dataFileNameTiming="${dataFileName}${DatSuffix}" # append to existing file for graphing : retain as spearate as may way to create one per iteration
if [[ ${#NXArray[*]} -le 0 ]] ; then 
	error 5 "matrix size"
else        
	executeOptions=""
	if [ "${compileUsingCblas}" == "true" ] ; then
	    compile_dgemm_cblas ${algorithmName}
	fi
	for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	    matrixFileNameValues="${matrixFileName}-$i${txtSuffix}" # different file for each run
	    init_log_file ${matrixFileNameValues} 
	    init_log_file ${dataFileNameTiming}
	    executeOptions="${algorithmOptions} ${NXArray[$i]}  ${threadArray[$i]}"
	    // echo "DEBUG : algorithm_execute  ./${algorithmName} ${executeOptions} ${matrixFileNameValues} ${dataFileNameTiming}"
	    algorithm_execute "./${algorithmName}" "${executeOptions}" "${matrixFileNameValues}" "${dataFileNameTiming}"
	done
fi

# move log files to <logDir>
mv -f *.dat ${logDir} 2> /dev/null
mv -f *.txt ${logDir} 2> /dev/null
mv -f *.log ${logDir} 2> /dev/null
mv -f *.bup ${logDir} 2> /dev/null

pause
exit 0
