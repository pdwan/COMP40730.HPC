#!/bin/bash 

# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# GIT : https://github.com/pdwan/COMP40730.HPC.git
# DUE DATE : 30-June-2014 (extended to : 13-July-2014)
# ASSIGNMENT : 1
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
stdLogFile="runAssignment1-${Now}.log"
compileUsingCblas="false"
compileUsingAtlas="false"
buildSimple="false"
buildBlockedIJK="false"
buildBlockedKIJ="false"
initRandom="false"
initIncrement="false"
matrixEnabled="false"
blockEnabled="false"
defaultMatrixRange="false"
let matrixSize=0
let maxMatrixSize=1000 
let blockSize=0
let maxBlockSize=100

# ##################################################################################

# function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <enter> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0) \ \n\t\t\t-a|--all -1|--simple -2|--ijk -3|--kij -d1|--atlas -d2|--cblas -r|--random -i|--increment \ \n\t\t\t-m|--matrix<n> -b|--block <b> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| using 1 -> 3 algoritms : Straight-forward IJK, Blocked IJK and Blocked KIJ." 
    $_ECHO -e "LOGS :\tCreated in current dir and moved to <${logDir}> : \n\t<file>.txt : \tmatrix values for matrices |A| |B| & |C|, \n\t<file>.dat :\ttiming data for each computation \n\t<file>.log : \tsummary of stdout. \n"
    $_ECHO -e "WHERE :\t-a|--all \tCalculate data for all algorithms via separate .c programs to multiply |A|x|B| -> |C| "
    $_ECHO -e "\t\t\tStraightforward IJK algorithm :\tA1-Sijk-1D.c \n\t\t\tBlocked IJK algorithm using square bxb blocks :\tA1-Bijk-1D.c \n\t\t\tBlocked KIJ algorithm using square bxb blocks :\tA1-Bkij-1D.c "
    $_ECHO -e "\t-1|--simple \tCalculate data for only the algorithm \n\t\t\tStraightforward IJK algorithm :\tA1-Sijk-1D.c "
    $_ECHO -e "\t-2|--bijk \tCalculate data for only the algorithm \n\t\t\tBlocked IJK algorithm using square bxb blocks :\tA1-Bijk-1D.c " 
    $_ECHO -e "\t-3|--bkij \tCalculate data for only the algorithm \n\t\t\tBlocked KIJ algorithm using square bxb blocks :\tA1-Bkij-1D.c  \n"
    $_ECHO -e "\t-d1|--atlas\tCompile .c source files using dgemm atlas"
    $_ECHO -e "\t-d2|--cblas\tCompile .c source files using dgemm cblas"
    $_ECHO -e "\t-r|--random \tInitialize |A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <column> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive \n"
    $_ECHO -e "\t-m|--matrix <N>\tMatrix size, if  invalid, (matrix size > max) set to [${maxMatrixSize}]"
    $_ECHO -e "\t-b|--block <B>\tBlock matrix size, if invalid (matrix % block != 0 or block size > max),  set to [${maxBlockSize}] and matrix size set to [${maxMatrixSize}]."
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t<NXArray> \t( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )"
    $_ECHO -e "\t\t\t<NBArray> \t( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 )"
    $_ECHO -e "\t\t\t'-m|--matrix <n>' & '-b|--block <b>' are mutually exclusive of '-v|--values'.\n"
    $_ECHO -e "\t-?|-h|--help \tusage \n"
}

# function : error message and then usage
error() 
{
    echo "error $1"
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

# function : build applying dgemm atlas
compile_dgemm_atlas()
{
    localProgramToCompile=$1
    $_ECHO -e "CBLAS :\t\tCompiling ${localProgramToCompile} using atlast \n"  >> ${stdLogFile}
    gcc -o ${localProgramToCompile} ${localProgramToCompile}.c -I/home/cs/khasanov/libs/ATLAS/include/ -L/home/cs/khasanov/libs/ATLAS/lib/Linux_UNKNOWNSSE2_4/ -lcblas -latlas -lm -O3
}

# function : build applying dgemm cblas
compile_dgemm_cblas()
{
    localProgramToCompile=$1
    $_ECHO -e "CBLAS :\t\tCompiling ${localProgramToCompile} using cblas \n"  >> ${stdLogFile}
    gcc -Wall -I/home/cs/khasanov/libs/CBLAS/src ${localProgramToCompile}.c -o ${localProgramToCompile}  /home/cs/khasanov/libs/cblas_LINUX.a  /usr/lib/libblas.a -lgfortran -fopenmp
}

# function : create directory, if  it does not exist & validate creation
init_dir() 
{
    creationDir=$1
    if  [ ! -d ${creationDir} ] || [ ! -e ${creationDir} ] ; then 
        mkdir ${creationDir}        
        $_ECHO -e "WARNING :\tCreating $creationDir"  >> ${stdLogFile}
        if  [[ $? -ne 0 ]] ; then 
            error 2 $creationDir
        fi 
    fi   
}

# function : create log files (.txt : matrix values, .dat : timing data for each computation & .log : stderr, stdout) to store values for data for each alogrithim computation
init_log_file() 
{
    localLogFile=$1
    if   [ -e ${localLogFile} ]  ; then
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
if  [[ $# -eq 0 ]] ; then
    usage 
    exit
else 
    while [ "$1" != "" ]; do
        case ${1} in
            "-a" | "--all")
                buildSimple="true"
                buildBlockedIJK="true"
                buildBlockedKIJ="true" 
                ;;
            "-1" | "--simple")
                buildSimple="true"
                ;;
            "-2" | "--ijk")
                buildBlockedIJK="true"
                ;;
            "-3" | "--kij")
                buildBlockedKIJ="true"
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
                defaultMatrixRange="true"
                ;;                    
            "-m" | "--matrix")
                matrixEnabled="true"
                if  [[ $2 =~ "^[ 0-9 ]+$" ]] ; then 
                    let matrixSize=$2 
                else
                    error 5 "${1} ${2}"
                fi  
                shift 
                ;;
            "-b" | "--block")
                blockEnabled="true"
                if  [[ $2 =~ "^[ 0-9 ]+$" ]] ; then 
                    let blockSize=$2 
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

#   process and validate parameter values for matrix and block sizes, if  applicable
if  [ "${defaultMatrixRange}" == "true" ] && ( [ "${matrixEnabled}" == "true" ] || [ "${blockEnabled}" == "true" ] ) ; then 
    error 6 "<-b>, <-m> & <-v>"
fi 
if   [ "${defaultMatrixRange}" == "true" ] && ( [ "${matrixEnabled}" == "false" ] || [ "${blockEnabled}" == "false" ] ) ; then 
    declare -a NXArray=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
    declare -a NBArray=( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 )
    matrixEnabled="false"
    blockEnabled="false"
elif  [ "${matrixEnabled}" == "true" ] || [ "${blockEnabled}" == "true" ] ; then
#   Validate matrix size initialized
	    if  [ "${matrixSize}" == "" ] ; then 
	        error 4 "matrix size"
    	fi 
#   Validate block size initialized if  applicable
	    if  [ "${blockSize}" == "" ]  && [ buildBlockedIJK="true" ]  ; then 
	        error 4 "block size"
	    fi 
	    if  [ "${blockSize}" == "" ]  && [ buildBlockedJIK="true" ]  ; then 
	        error 4 "block size"
	    fi 
#   Validate matrix range
	    if  [ ${matrixSize} -le 0 ] || [ ${matrixSize} -gt ${maxMatrixSize} ] ; then
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): Invalid matrix size <nx>, now set to default of : $maxMatrixSize"  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize
	        matrixEnabled="true"
	    fi         
#   Validate block range
        if  [ ${blockSize} -le 0 ] || [ ${blockSize} -gt ${maxBlockSize} ] ; then
    	    $_ECHO -e "WARNING :\t$($_BASENAME $0): Invalid block size <nb>, now set to default of : $maxBlockSize"  >> ${stdLogFile}
    	    let blockSize=$maxBlockSize
    	    blockEnabled="true"
    	fi 
#   ensure both are enabled if  one is unless Simple only
	    if  [ "${matrixEnabled}" == "false" ] ; then 
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is now enabled and set to default of : $maxMatrixSize"  >> ${stdLogFile}
	        let matrixSize=$maxMatrixSize
	        matrixEnabled="true"
    	fi 
	    if  [ "${blockEnabled}" == "false" ] && ( [ "${buildBlockedIJK}"="true" ] || [ "${buildBlockedKIJ}"="true" ] ) ; then 
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): block size <nb> is now enabled and set to default of : $maxBlockSize"  >> ${stdLogFile}
	        let blockSize=$maxBlockSize
	        blockEnabled="true"
	    fi 
#   Validate matrix / block remainder
        if  [ $(( ${matrixSize} % ${blockSize} )) -ne 0 ]; then
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): block size needs to be an even multiple of matrix size. Using defaults."  >> ${stdLogFile}
	        let blockSize=$maxBlockSize
	        let matrixSize=$maxMatrixSize
    	fi 
#   Validate Simple IJK & set matrix and block arrays
        if  [ buildSimple="true" ] && [ buildBlockedIJK="false" ]  && [ buildBlockedKIJ="false" ] ; then 
            if  [ "${matrixEnabled}" == "true" ] && [ "${blockEnabled}" == "true" ] ; then
                error 6 "Simple IJK only : <-m> & <-b>"
            fi 
            declare -a NXArray=( $matrixSize )
        else 
            declare -a NXArray=( $matrixSize )
	        declare -a NBArray=( $blockSize )
        fi 
fi 

# validate atlas and cblas - mutually exclusive
if  [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
        error 7 "Atlas & cBlas"
fi 

# build up commands to run
algorithmOptions=""
if  [ "${initRandom}" == "true" ] && [ "${initIncrement}" == "false" ] ; then 
    algorithmOptions="-r"
fi 
if  [ "${initRandom}" == "false" ] && [ "${initIncrement}" == "true" ] ; then 
    algorithmOptions="-i"
fi 
if  [ "${initRandom}" == "true" ] && [ "${initIncrement}" == "true" ] ; then 
    error 6 "<-i> & <-r>" 
fi 

# execute algorithms

init_dir ${logDir}
init_log_file ${stdLogFile}

matrixFileRoot="${logPrefix}${Now}-values"
dataFileRoot="${logPrefix}${Now}-data"

# build Simple IJK
if  [ "${buildSimple}" == "true" ]  || [ "${buildAll}" == "true" ] ; then
    algorithmSimple="A1-Sijk-1D" 
    matrixFileSimple="${matrixFileRoot}-${algorithmSimple}"
    dataFileSimple="${dataFileRoot}-${algorithmSimple}"
    dataFileSimpleTiming="${dataFileSimple}${DatSuffix}" # append to existing file for graphing
    if  [[ ${#NXArray[*]} -le 0 ]] ; then 
	    error 5 "matrix size"
    else        
	    executeOptions=""
        if  [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmSimple}
	    fi 
	    if  [ "${compileUsingAtlas}" == "true" ] ; then
	        compile_dgemm_atlas ${algorithmSimple}
	    fi 
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	        matrixFileSimpleValues="${matrixFileSimple}-$i${txtSuffix}" # different file for each iteration
	        init_log_file ${matrixFileSimpleValues} 
	        init_log_file ${dataFileSimpleTiming}
	        executeOptions="${algorithmOptions} ${NXArray[$i]}  ${threadArray[$i]}"
            echo " DEBUG : algorithm_execute ./${algorithmSimple} ${executeOptions} ${matrixFileSimpleValues} ${dataFileSimpleTiming}"
	        algorithm_execute "./${algorithmSimple}" "${executeOptions}" "${matrixFileSimpleValues}" "${dataFileSimpleTiming}"
	        matrixFileSimpleValues="${matrixFileSimple}"
	    done
    fi 
fi 

# build Blocked IJK
if  [ "${buildBlockedIJK}" == "true" ]  || [ "${buildAll}" == "true" ] ; then
    algorithmBlockIJK="A1-Bijk-1D" 
    matrixFileBlockIJK="${matrixFileRoot}-${algorithmBlockIJK}"
    dataFileBlockIJK="${dataFileRoot}-${algorithmBlockIJK}"
    dataFileBlockIJKTiming="${dataFileBlockIJK}${DatSuffix}" # append to existing file for graphing
    if  [[ ${#NXArray[*]} -le 0 ]] ; then 
	    error 5 "matrix size"
    else        
	    executeOptions=""
        if  [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmBlockIJK}
	    fi 
	    if  [ "${compileUsingAtlas}" == "true" ] ; then
	        compile_dgemm_atlas ${algorithmBlockIJK}
	    fi 
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	        matrixFileBlockIJKValues="${matrixFileBlockIJK}-$i${txtSuffix}" # different file for each iteration
	        init_log_file ${matrixFileBlockIJKValues} 
	        init_log_file ${dataFileBlockIJKTiming}
	        executeOptions="${algorithmOptions} ${NXArray[$i]}  ${threadArray[$i]}"
	        echo "DEBUG : algorithm_execute ./${algorithmBlockIJK} ${executeOptions} ${matrixFileBlockIJKValues} ${dataFileBlockIJKTiming}"
	        algorithm_execute "./${algorithmBlockIJK}" "${executeOptions}" "${matrixFileBlockIJKValues}" "${dataFileBlockIJKTiming}"
	        matrixFileBlockIJKValues="${matrixFileBlockIJK}"
	    done
    fi 
fi 

# build Blocked KIJ
if  [ "${buildBlockedKIJ}" == "true" ]  || [ "${buildAll}" == "true" ] ; then
    algorithmBlockKIJ="A1-Bkij-1D" 
    matrixFileBlockKIJ="${matrixFileRoot}-${algorithmBlockKIJ}"
    dataFileBlockKIJ="${dataFileRoot}-${algorithmBlockKIJ}"
    dataFileBlockKIJTiming="${dataFileBlockKIJ}${DatSuffix}" # append to existing file for graphing
    if  [[ ${#NXArray[*]} -le 0 ]] ; then 
	    error 5 "matrix size"
    else        
	    executeOptions=""
        if  [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmBlockKIJ}
	    fi 
	    if  [ "${compileUsingAtlas}" == "true" ] ; then
	        compile_dgemm_atlas ${algorithmBlockKIJ}
	    fi 
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )) ; do
	        matrixFileBlockKIJValues="${matrixFileBlockKIJ}-$i${txtSuffix}" # different file for each iteration
	        init_log_file ${matrixFileBlockKIJValues} 
	        init_log_file ${dataFileBlockKIJTiming}
	        executeOptions="${algorithmOptions} ${NXArray[$i]} ${threadArray[$i]}"
	        echo "DEBUG : algorithm_execute ./${algorithmBlockKIJ} ${executeOptions} ${matrixFileBlockKIJValues} ${dataFileBlockKIJTiming}"
	        algorithm_execute "./${algorithmBlockKIJ}" "${executeOptions}" "${matrixFileBlockKIJValues}" "${dataFileBlockKIJTiming}"
	        matrixFileBlockKIJValues="${matrixFileBlockKIJ}"
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
