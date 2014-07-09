#!/bin/bash 

# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# GIT : https://github.com/pdwan/COMP40730.HPC.git
# DUE DATE : 30-June-2014 (extended to : 09-July-2014)
# ASSIGNMENT : 2
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
stdLogFile="runAssignment2-${Now}.log"
compileUsingAtlas="false"
compileUsingCblas="false"
plotGraphUsingGnuPlot="false"
initRandom="false"
initIncrement="false"
matrixEnabled="false"
defaultMatrixRange="false"
let matrixSize=0
let maxMatrixSize=1000 
let maxThreads=124

# ##################################################################################

# function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <enter> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0)  -d1|--atlas -p|--plot -r|--random -i|--increment -m|--matrix<n> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| and then infinity norm using pthreads. \n" 
    $_ECHO -e "LOGS :\tCreated in <${logDir}> : <file>.txt matrix values for matrices |A| |B| & |C|, \n\t<file>.dat : timing of each computation & <${logDir}/runAssignment1-timestamp.log> summary of stdout. \n"
    $_ECHO -e "WHERE :\t-d1|--atlas\tCompile .c source files using dgemm ATLAS "
    $_ECHO -e "\t-d2|--cblas\tCompile .c source files using dgemm cBLAS"
    $_ECHO -e "\t\t\tEach is mutually exclusive of the other" 
    $_ECHO -e "\t-p|--plot\tPlot graphs using GnuPlot creating .png for each algorithm and store in <${logDir}> for matrix size -v- time taken \n"
    $_ECHO -e "\t-r|--random \tInitialize |A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <row> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive \n"
    $_ECHO -e "\t-m|--matrix <n>\tMatrix dimension, if odd number +1 added or if invalid set to '1,000', thread count set to { 2 }  "
    $_ECHO -e "\t\t\tMutually exclusive with '-v|--values'" 
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t<nx> \t\t{ 50, 50, 50, 100, 100, 100, 500, 500, 500, 1000, 1000, 1000 } and \n\t\t\t<threadArray> \t{  10 10 10 10 10 10 20 20 20 20 20 20 20 20  } "
    $_ECHO -e "\t\t\t'-m|--matrix <n>' & '-v|--values' mutually exclusive. \n"
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

    $_ECHO -e "# -----------------------------------------------------------------------------------------------------------------------------------------  \n# " >> ${localLogFile}
    $_ECHO -e "# Program : ${localProgramName} \n# \n# Log file : ${localLogFile} \n# where :\t.dat contains timing data & .txt contains matrix values \n#" >> ${localLogFile}
    $_ECHO -e "# Calculate |C| = |A| x |B| and then infinity norm using pthreads.| \n# " >> ${localLogFile}
    $_ECHO -e "# -----------------------------------------------------------------------------------------------------------------------------------------  \n# " >> ${localLogFile}
    if [[ $localLogFile == *"$FileTypeDat"* ]] ; then
	    $_ECHO -e "# Time taken to compute \n#" >> ${localLogFile} # dat
	    $_ECHO -e "# Matrix Size \tNo Threads \tTime/manual \tTime/dgenn \n# " |& $_TEE -a ${localLogFile}
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

# process and validate parameter values for matrix, if applicable
if [ "${defaultMatrixRange}" == "true" ] && [ "${matrixEnabled}" == "true" ] ; then 
    error 6 "<-m> & <-v>"
fi
if  [ "${defaultMatrixRange}" == "true" ] && [ "${matrixEnabled}" == "false" ]  ; then 
    declare -a NXArray=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
	declare -a threadArray=( 10 10 10 20 20 20 50 50 50 50 50 50 50 50 )
	matrixEnabled="false"
else
    if [ "${matrixEnabled}" == "true" ] ; then
	    if [ "${matrixSize}" == "" ] ; then 
	        error 4 "matrix size"
    	fi
	    if [ ${matrixSize} -le 0 ] || [ ${matrixSize} -gt ${maxMatrixSize} ] ; then
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): Invalid matrix size <nx>, now set to default of : $maxMatrixSize" |& $_TEE -a ${logDir}/${stdLogFile}
	        let matrixSize=$maxMatrixSize
	        matrixEnabled="true"
	    fi        
        if [ ( ${matrixSize} % 2 ) == 0 ] ; then 
                let matrixSize=$matrixSize 
            else 
                let matrixSize=$(( matrixSize++ )) 
        fi  
		declare -a NXArray=( $matrixSize )
		declare -a threadArray=( 2 )
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

algorithmName="A2-pthreads-1D" 
matrixFileName="${matrixFileRoot}-${algorithmName}"
dataFileName="${dataFileRoot}-${algorithmName}"
graphFileName="${graphFileRoot}-${algorithmName}"

if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
    error 5 "matrix size"
else        
    executeOptions=""
    if [ "${compileUsingAtlas}" == "true" ] && [ "${compileUsingCblas}" == "true" ] ; then 
        error 7 "Atlas & cBlas"
    fi
    if [ "${compileUsingAtlas}" == "true" ] ; then
	    compile_dgemm_atlas ${algorithmName}
    elif [ "${compileUsingCblas}" == "true" ] ; then
	    compile_dgemm_cblas ${algorithmName}
	fi
    for (( i = 0 ; i < ${#NXArray[@]} ; i++ )); do
        matrixFileNameValues="${matrixFileName}-$i${txtSuffix}"
	    dataFileNameTiming="${dataFileName}-$i${DatSuffix}"
	    graphFileNamePlotted="${graphFileName}-$i${pngSuffix}"
	    init_log_file ${matrixFileNameValues} 
	    init_log_file ${dataFileNameTiming}
	    add_comments_to_log_file "${matrixFileNameValues}" "${algorithmName}"
	    add_comments_to_log_file "${dataFileNameTiming}" "${algorithmName}"
	    executeOptions="${algorithmOptions} ${NXArray[$i]} ${threadArray[$i]}"
	    algorithm_execute "${algorithmName}" "${executeOptions}" "${matrixFileNameValues}" "${dataFileNameTiming}"
 	    if [ "${plotGraphUsingGnuPlot}" == "true" ] ; then 
	        plot_graph ${dataFileNameTiming} ${graphFileNamePlotted}
	    fi  
	    matrixFileNameValues="${matrixFileName}"
	    dataFileNameTiming="${dataFileName}"
	    graphFileNamePlotted="${graphFileName}"
	done
fi

pause
exit 0
