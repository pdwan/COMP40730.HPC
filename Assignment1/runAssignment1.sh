#!/bin/bash 

# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# GIT : https://github.com/pdwan/COMP40730.HPC.git
# DUE DATE : 30-June-2014 (extended to : 09-July-2014)
# ASSIGNMENT : 1
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
_TEE="tee"
Now=$(date +"%Y%m%d.%H%M%S")
logDir="logDir"
graphDir="graphDir"
logPrefix="pdwan-"
txtSuffix=".txt"
DatSuffix=".dat"
pngSuffix=".png"
stdLogFile="runAssignment1-${Now}.log"
buildAll="false"
buildSimple="false" 
buildBlockedIJK="false"
buildBlockedKIJ="false" 
compileUsingAtlas="false"
compileUsingCblas="false"
plotGraphUsingGnuPlot="false"
initRandom="false"
initIncrement="false"
matrixEnabled="false"
blockEnabled="false"
defaultMatrixRange="false"
let matrixSize=0
let blockSize=0 
let maxMatrixSize=1000 
let maxBlockSize=50

# ##################################################################################

# function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <enter> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0) -a|--all -1|--simple -2|--ijk -3|--kij -d1|--atlas -p|--plot -r|--random -i|--increment \ \n\t\t\t-m|--matrix<n> -b|--block <b> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| using 1 -> 3 algoritms : Straight-forward IJK, Blocked IJK and Blocked KIJ." 
    $_ECHO -e "LOGS :\tCreated in <${logDir}> : <file>.txt matrix values for matrices |A| |B| & |C|, \n\t<file>.dat : timing of each computation & <${logDir}/runAssignment1-timestamp.log> summary of stdout. \n"
    $_ECHO -e "WHERE :\t-a|--all \tCalculate data for all algorithms via separate .c programs to multiply |A|x|B| -> |C| "
    $_ECHO -e "\t\t\tStraightforward IJK algorithm :\tA1-Sijk-1D.c \n\t\t\tBlocked IJK algorithm using square bxb blocks :\tA1-Bijk-1D.c \n\t\t\tBlocked KIJ algorithm using square bxb blocks :\tA1-Bkij-1D.c "
    $_ECHO -e "\t-1|--simple \tCalculate data for only the algorithm \n\t\t\tStraightforward IJK algorithm :\tA1-Sijk-1D.c "
    $_ECHO -e "\t-2|--bijk \tCalculate data for only the algorithm \n\t\t\tBlocked IJK algorithm using square bxb blocks :\tA1-Bijk-1D.c " 
    $_ECHO -e "\t-3|--bkij \tCalculate data for only the algorithm \n\t\t\tBlocked KIJ algorithm using square bxb blocks :\tA1-Bkij-1D.c  \n"
    $_ECHO -e "\t-d1|--atlas\tCompile .c source files using dgemm ATLAS "
    $_ECHO -e "\t-d2|--cblas\tCompile .c source files using dgemm cBLAS"
    $_ECHO -e "\t\t\tEach is mutually exclusive of the other" 
    $_ECHO -e "\t-p|--plot\tPlot graphs using GnuPlot creating .png for each algorithm and store in <${logDir}> for \n\t\t\t(i) matrix size -v- time taken  &  (ii) block size -v- time taken \n"
    $_ECHO -e "\t-r|--random \tInitialize |A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <row> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive \n"
    $_ECHO -e "\t-m|--matrix <n>\tMatrix dimension, if invalid set to '1,000' "
    $_ECHO -e "\t-b|--block <b> \tBlock size, with '<nb> < <nx>' & '<nx> % <nb> = 0', if invalid set to '50'. <nx> set to '1,000' " 
    $_ECHO -e "\t\t\tMutually exclusive with '-v|--values'" 
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t<nx> = { 50, 50, 50, 100, 100, 100, 500, 500, 500, 1000, 1000, 1000 } "
    $_ECHO -e "\t\t\t<nb> = { 2, 5, 10, 5, 10, 20, 10, 20, 50, 10, 50, 100 } "
    $_ECHO -e "\t\t\tMutually exclusive with '-m|--matrix <n>' and '-b|--block <b>' \n"
    $_ECHO -e "\t-?|-h|--help \tusage \n"
}

# function : error message and then usage
error() 
{
    err=$1
    mess=$2
    case ${err} in
        1)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Unknown parameter '${mess}'." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;
        2)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating directory '${mess}'." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;
        3)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating file '${mess}'." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        4)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Missing parameter for '${mess}'." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        5)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${mess}', Values entered are not valid or not a number." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        6)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${mess}', Mutually exclusive switches." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        7)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${mess}', Compilation failed." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        8)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${mess}', GnuPlot - graph export error." |& $_TEE -a ${logDir}/${stdLogFile}
            ;;     
        *)
            $_ECHO -e "ERROR :\tUnknown error." |& $_TEE -a ${logDir}/${stdLogFile}
            $_ECHO ${err} |& $_TEE -a ${logDir}/${stdLogFile}
            $_ECHO $* |& $_TEE -a ${logDir}/${stdLogFile}
            ;;
    esac
    $_ECHO -e
    usage
    exit ${err}
}

#function : build applying dgemm atlas
compile_dgemm_atlas()
{
    localProgramToCompile=$1
    $_ECHO -e "ATLAS :\t\tCompiling ${localProgramToCompile} using atlas \n" |& $_TEE -a ${logDir}/${stdLogFile}
    gcc -o ATLAS ${localProgramToCompile}.c -I/home/cs/khasanov/libs/ATLAS/include/ -L/home/cs/khasanov/libs/ATLAS/lib/Linux_UNKNOWNSSE2_4/ -lcblas -latlas -lm -O3
}

# function : build applying dgemm cblas
compile_dgemm_cblas()
{
    localProgramToCompile=$1
    $_ECHO -e "CBLAS :\t\tCompiling ${localProgramToCompile} using cblas \n" |& $_TEE -a ${logDir}/${stdLogFile}
    gcc -Wall -I/home/cs/khasanov/libs/CBLAS/src ${localProgramToCompile}.c -o ${localProgramToCompile}  /home/cs/khasanov/libs/cblas_LINUX.a  /usr/lib/libblas.a -lgfortran
}

# function : plot graph using GnuPlot
plot_graph()
{
    localDatFileToGraph=$1
    localPngGraph=$2
    $_ECHO -e "GNUplot :\tCreating ${localPngGraph} from ${localDatFileToGraph}" |& $_TEE -a ${logDir}/${stdLogFile}
}

# function : create directory, if it does not exist & validate creation
init_dir() 
{
    creationDir=$1
    if [ ! -d ${creationDir} ] || [ ! -e ${creationDir} ] ; then 
        mkdir ${creationDir}        
        $_ECHO -e "WARNING :\tCreating $creationDir" |& $_TEE -a ${logDir}/${stdLogFile}
        if [[ $? -ne 0 ]] ; then 
            error 2 $creationDir
        fi
    fi  
}

# function : create log files (.txt : matrix values, .dat : timing of each computation & .log : stderr, stdout) to store values for data for each alogrithim computation
init_log_file() 
{
    localLogFile=$1
    if  [ -e ${localLogFile} ]  ; then
        $_ECHO -e "WARNING :\tFile backup : ${localLogFile} to ${localLogFile}.bup" 
        mv "${localLogFile}" "${localLogFile}.bup"
    fi
    $_ECHO -e "# LOG FILE :\t${localLogFile} \tcreated on ${Now} by ${USER}." |& $_TEE ${localLogFile}
}

# function : add initial comments to matrix .txt and to timing .dat file
add_comments_to_log_file() 
{
    localLogFile=$1
    localProgramName=$2
    FileTypeDat='dat'
    FileTypeSimple='Sijk'

    $_ECHO -e "# -----------------------------------------------------------------------------------------------------------------------------------------  \n# " >> ${localLogFile}
    $_ECHO -e "# Program : ${localProgramName} \n# \n# Log file : ${localLogFile} \n# where :\t.dat contains timing data & .txt contains matrix values \n#" >> ${localLogFile}
    $_ECHO -e "# Computation of Straight-forward IJK, Blocked IJK or Blocked KIJ using square NxN \n# & square bxb block (as applicable) for matrices |C| += |A| * |B| \n# " >> ${localLogFile}
    $_ECHO -e "# -----------------------------------------------------------------------------------------------------------------------------------------  \n# " >> ${localLogFile}
    if [[ $localLogFile == *"$FileTypeDat"* ]] ; then
	    $_ECHO -e "# Time taken to compute \n#" >> ${localLogFile} # dat
	    if [[ "${localLogFile}" == *"$FileTypeSimple"* ]] ; then
	        $_ECHO -e "# Matrix Size \tTime/Simple \tTime/Complex \tTime/dgenn \n# " |& $_TEE -a ${localLogFile}
	    else 
	        $_ECHO -e "# Matrix Size \tBlock Size \tTime/Simple \tTime/Complex \tTime/dgenn \n# " |& $_TEE -a ${localLogFile}
        fi
    else 
        $_ECHO -e "# Summary of values added to each matrix - retained for later reference and validation \n#" >> ${localLogFile} # txt 
    fi
}

# function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
    localCmd="$1"
    localOptions="$2"
    localFileMatrix="$3"
    localFileTime="$4"
    $_ECHO -e "RUNNING :\t${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}" |& $_TEE -a ${logDir}/${stdLogFile}
    # ${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}
}

# ##################################################################################

clear

# Process parameters
if [[ $# -eq 0 ]]; then
    usage 
    exit
else 
    while [ "$1" != "" ]; do
        case ${1,,} in
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
		    "-p" | "--plot")
				plotGraphUsingGnuPlot="true"
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
                if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    let matrixSize=$2 
                else 
                    error 5 "${1} ${2}" 
                fi                    
                shift 
                ;;
            "-b" | "--block")
                blockEnabled="true"
                if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    let matrixSize=$2 
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

# process and validate parameter values for matrix and block sizes, if applicable
if [ "${defaultMatrixRange}" == "true" ] && ( [ "${matrixEnabled}" == "true" ] || [ "${blockEnabled}" == "true" ] ) ; then 
    error 6 "<-b>, <-m> & <-v>"
fi
if  [ "${defaultMatrixRange}" == "true" ] && ( [ "${matrixEnabled}" == "false" ] || [ "${blockEnabled}" == "false" ] ) ; then 
    declare -a NXArray=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
    declare -a NBArray=( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 )
    matrixEnabled="false"
    blockEnabled="false"
else
    if [ "${matrixEnabled}" == "true" ] || [ "${blockEnabled}" == "true" ] ; then
	if [ "${matrixSize}" == "" ] ; then 
	    error 4 "matrix size"
	fi
	if [ "${blockSize}" == "" ] ; then 
	    error 4 "block size"
	fi
	# Validate matrix range
	if [ ${matrixSize} -le 0 ] || [ ${matrixSize} -gt ${maxMatrixSize} ] ; then
	    $_ECHO -e "WARNING :\t$($_BASENAME $0): Invalid matrix size <nx>, now set to default of : $maxMatrixSize" |& $_TEE -a ${logDir}/${stdLogFile}
	    let matrixSize=$maxMatrixSize
	    matrixEnabled="true"
	fi        
	# Validate block range
	if [ ${blockSize} -le 0 ] || [ ${blockSize} -gt ${maxBlockSize} ] ; then
	    $_ECHO -e "WARNING :\t$($_BASENAME $0): Invalid block size <nb>, now set to default of : $maxBlockSize" |& $_TEE -a ${logDir}/${stdLogFile}
	    let blockSize=$maxBlockSize
	    blockEnabled="true"
	fi
	# ensure both are enabled if one is unless Simple only
	if [ "${matrixEnabled}" == "false" ] ; then 
	    $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is now enabled and set to default of : $maxMatrixSize" |& $_TEE -a ${logDir}/${stdLogFile}
	    let matrixSize=$maxMatrixSize
	    matrixEnabled="true"
	fi
	if [ "${blockEnabled}" == "false" ] && ( ["${buildBlockedIJK}"="true"] || ["${buildBlockedKIJ}"="true" ] ) ; then 
	    $_ECHO -e "WARNING :\t$($_BASENAME $0): block size <nb> is now enabled and set to default of : $maxBlockSize" |& $_TEE -a ${logDir}/${stdLogFile}
	    let blockSize=$maxBlockSize
	    blockEnabled="true"
	fi
	# Validate matrix / block remainder
    	if [ $(( ${matrixSize} % ${blockSize} )) -ne 0 ]; then
	    $_ECHO -e "WARNING :\t$($_BASENAME $0): block size needs to be an even multiple of matrix size. Using defaults." |& $_TEE -a ${logDir}/${stdLogFile}
	    let blockSize=$maxBlockSize
	    let matrixSize=$maxMatrixSize
    	fi
	declare -a NXArray=( $matrixSize )
	declare -a NBArray=( $blockSize )
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
    error 6 "<-i> and <-r>" 
fi

# execute algorithms

init_dir ${logDir}
init_dir ${graphDir}
init_log_file ${logDir}/${stdLogFile}

matrixFileRoot="${logDir}/${logPrefix}${Now}-values"
dataFileRoot="${logDir}/${logPrefix}${Now}-timing"
graphFileRoot="${graphDir}/${logPrefix}${Now}-graph"


if [ "${buildSimple}" == "true" ]  || [ "${buildAll}" == "true" ] ; then
    algorithmSimple="A1-Sijk-1D" 
    matrixFileSimple="${matrixFileRoot}-${algorithmSimple}"
    dataFileSimple="${dataFileRoot}-${algorithmSimple}"
    graphFileSimple="${graphFileRoot}-${algorithmSimple}"
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
	    error 5 "matrix size and block size"
    else        
	    executeOptions=""
	    if [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
            error 7 "Atlas & cBlas"
        fi
        if [ "${compileUsingAtlas}" == "true" ] ; then
	        compile_dgemm_atlas ${algorithmSimple}
	    elif [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmSimple}
	    fi
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	        matrixFileSimpleValues="${matrixFileSimple}-$i${txtSuffix}"
	        dataFileSimpleTiming="${dataFileSimple}-$i${DatSuffix}"
	        graphFileSimplePlotted="${graphFileSimple}-$i${pngSuffix}"
	        init_log_file ${matrixFileSimpleValues} 
	        init_log_file ${dataFileSimpleTiming}
	        add_comments_to_log_file "${matrixFileSimpleValues}" "${algorithmSimple}"
	        add_comments_to_log_file "${dataFileSimpleTiming}" "${algorithmSimple}"
	        executeOptions="${algorithmOptions} ${NXArray[$i]}"
	        algorithm_execute "${algorithmSimple}" "${executeOptions}" "${matrixFileSimpleValues}" "${dataFileSimpleTiming}"
 	        if [ "${plotGraphUsingGnuPlot}" == "true" ] ; then 
	    	    plot_graph ${dataFileSimpleTiming} ${graphFileSimplePlotted}
	        fi  
	        matrixFileSimpleValues="${matrixFileSimple}"
	        dataFileSimpleTiming="${dataFileSimple}"
	        graphFileSimplePlotted="${graphFileSimple}"
	    done
    fi
fi

if [ "${buildBlockedIJK}" == "true" ]  || [ "${buildAll}" == "true" ] ; then 
    algorithmBlockedIJK="A1-Bijk-1D"
    matrixFileBlockedIJK="${matrixFileRoot}-${algorithmBlockedIJK}"
    dataFileBlockedIJK="${dataFileRoot}-${algorithmBlockedIJK}"
    graphFileBlockedIJK="${graphFileRoot}-${algorithmBlockedIJK}"
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
	    error 5 "matrix size and block size"
    else 
	    executeOptions=""
	    if [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
            error 7 "Atlas & cBlas"
        fi
	    if [ "${compileUsingAtlas}" == "true" ] ; then 
	        compile_dgemm_atlas ${algorithmBlockedIJK}
	    elif [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmBlockedIJK}
	    fi
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	        matrixFileBlockedIJKValues="${matrixFileBlockedIJK}-$i${txtSuffix}"
	        dataFileBlockedIJKTiming="${dataFileBlockedIJK}-$i${DatSuffix}"
	        graphFileBlockedIJKPlotted="${graphFileBlockedIJK}-$i${pngSuffix}"
	        init_log_file $matrixFileBlockedIJKValues 
	        init_log_file $dataFileBlockedIJKTiming
	        add_comments_to_log_file "${matrixFileBlockedIJKValues}" "${algorithmBlockedIJK}"
	        add_comments_to_log_file "${dataFileBlockedIJKTiming}" "${algorithmBlockedIJK}"
	        executeOptions="${algorithmOptions} ${NXArray[$i]} ${NBArray[$i]}"
	        algorithm_execute "${algorithmBlockedIJK}"       "${executeOptions}"             \
                                  "${matrixFileBlockedIJKValues}"  "${dataFileBlockedIJKTiming}"
    	    if [ "${plotGraphUsingGnuPlot}" == "true" ] ; then 
	    	    plot_graph ${dataFileBlockedIJKTiming} ${graphFileBlockedIJKPlotted}
	        fi  
	        matrixFileBlockedIJKValues="${matrixFileBlockedIJK}"
	        dataFileBlockedIJKTiming="${dataFileBlockedIJK}"
	        graphFileBlockedIJKPlotted="${graphFileBlockedIJK}"
	    done
    fi
fi

if [ "${buildBlockedKIJ}" == "true" ]  || [ "${buildAll}" == "true" ] ; then 
    algorithmBlockedKIJ="A1-Bkij-1D"
    matrixFileBlockedKIJ="${matrixFileRoot}-${algorithmBlockedKIJ}"
    dataFileBlockedKIJ="${dataFileRoot}-${algorithmBlockedKIJ}"
    graphFileBlockedKIJ="${graphFileRoot}-${algorithmBlockedKIJ}"
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
	    error 5 "matrix size and block size"
    else 
	    executeOptions=""
	    if [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
            error 7 "Atlas & cBlas"
        fi
	    if [ "${compileUsingAtlas}" == "true" ] ; then 
	        compile_dgemm_atlas ${algorithmBlockedKIJ}
	    elif [ "${compileUsingCblas}" == "true" ] ; then
	        compile_dgemm_cblas ${algorithmBlockedKIJ}
	    fi
	    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
	        matrixFileBlockedKIJValues="${matrixFileBlockedKIJ}-$i${txtSuffix}"
    	    dataFileBlockedKIJTiming="${dataFileBlockedKIJ}-$i${DatSuffix}"
	        graphFileBlockedKIJPlotted="${graphFileBlockedKIJ}-$i${pngSuffix}"
	        init_log_file ${matrixFileBlockedKIJValues} 
	        init_log_file ${dataFileBlockedKIJTiming}
	        add_comments_to_log_file "${matrixFileBlockedKIJValues}" "${algorithmBlockedKIJ}" 
	        add_comments_to_log_file "${dataFileBlockedKIJTiming}" "${algorithmBlockedKIJ}" 
	        executeOptions="${algorithmOptions} ${NXArray[$i]} ${NBArray[$i]}"
	        algorithm_execute "${algorithmBlockedKIJ}" "${executeOptions}" \
                            "${matrixFileBlockedKIJValues}" "${dataFileBlockedKIJTiming}"
	        if [ "${plotGraphUsingGnuPlot}" == "true" ] ; then 
	    	    plot_graph ${dataFileBlockedKIJTiming} ${graphFileBlockedKIJPlotted}
	        fi
	        matrixFileBlockedKIJValues="${matrixFileBlockedKIJ}"
	        dataFileBlockedKIJTiming="${dataFileBlockedKIJ}"
	        graphFileBlockedKIJPlotted="${graphFileBlockedKIJ}"
	    done
    fi
fi

pause
exit 0
