
# ##################################################################################
# 
# DESC : Script to multiply two nxn matrices using three algorithms.
# AUTHOR : Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# DATE : 30-June-2014
# ASSIGNMENT : 1
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
_TEE="tee"
Now=$(date +"%Y%m%d-%H%M%S")
LogDir="logDir"
LogPrefix="pdwan-"
TxtSuffix=".txt"
DatSuffix=".dat"
StdLogFile=$(_BASENAME $0 +"-${Now}.log")
BuildSimple="false" 
BuildBlockedIJK="false"
BuildBlockedKIJ="false" 
InitRandom="false"
InitIncrement="false"
MatrixEnabled="false"
BlockEnabled="false"
DefaultMatrixRange="false"
let MatrixSize=0
let BlockSize=0 
let MaxMatrixSize=1000 
let MaxBlockSize=50

# ##################################################################################

# function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <any key> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0) \ \n\t\t\t-a|--all -1|--simple -2|--ijk -3|--kij -r|--random -i|--increment \ \n\t\t\t-m|--matrix<n> -b|--block <b> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO :\tCalculate |C| = |A| x |B| using one or three algoritms : Straight-forward IJK, Blocked IJK and Blocked KIJ. \n"   
    $_ECHO -e "WHERE:\t-a | --all \tCalculate data for all algorithms via separate .c programs to multiply |A|x|B| -> |C|"
    $_ECHO -e "\t\t\tStraightforward non-blocked ijk algorithm :\tA1-Sijk-1D.c \n\t\t\tBlocked ijk algorithm using square bxb blocks :\tA1-Bijk-1D.c \n\t\t\tBlocked kij algorithm using square bxb blocks :\tA1-Bkij-1D.c"
    $_ECHO -e "\t-1|--simple \tCalculate data for only the algorithm \n\t\t\tStraightforward non-blocked ijk algorithm :\tA1-Sijk-1D.c "
    $_ECHO -e "\t-2|--bijk \tCalculate data for only the algorithm \n\t\t\tBlocked ijk algorithm using square bxb blocks :\tA1-Bijk-1D.c " 
    $_ECHO -e "\t-3|--bkij \tCalculate data for only the algorithm \n\t\t\tBlocked kij algorithm using square bxb blocks :\tA1-Bkij-1D.c \n"
    $_ECHO -e "\t-r|--random \tInitialization A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment \tInitialize |A| & |B| incrementally with <row> value and |C| with '0' "
    $_ECHO -e "\t\t\t'-i|--increment' & '-r|--random' are mutually exclusive. \n"
    $_ECHO -e "\t-m|--matrix <n>\tMatrix dimension, set to maximum of '1,000' if invalid or not provided"
    $_ECHO -e "\t-b|--block <b> \tBlock size, with '<nb> < <nx>' and '<nx> % <nb> = 0', \n\t\t\tdefaults to '50' if invalid or not provided and <nx> will be reset to 1,000." 
    $_ECHO -e "\t\t\tMutually exclusive with '-v|--values' \n"    
    $_ECHO -e "\t-v|--values \tUse predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t<nx> = { 50, 50, 50, 100, 100, 100, 500, 500, 500, 1000, 1000, 1000 } "
    $_ECHO -e "\t\t\t<nb> = { 2, 5, 10, 5, 10, 20, 10, 20, 50, 10, 50, 100 } "
    $_ECHO -e "\t\t\tMutually exclusive with '-m|--matrix <n>' and '-b|--block <b>' \n"
    $_ECHO -e "\t-?|-h|--help \tusage"
    $_ECHO ""
}

# function : error message and then usage
error() 
{
    err=$1
    
    shift
    case ${err} in
        1)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Unknown parameter '${1}'."  2>&1 |& $_TEE -a ${StdLogFile}
            ;;
        2)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating directory '${1}'."  2>&1 |& $_TEE -a ${StdLogFile}
            ;;
        3)
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Error creating file '${1}'."  2>&1 |& $_TEE -a ${StdLogFile}
            ;;     
        4)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): Missing parameter for '${1}'."  2>&1 |& $_TEE -a ${StdLogFile}
            ;;     
        5)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Values entered are not valid or not a number."  2>&1 |& $_TEE -a ${StdLogFile}
            ;;     
        6)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Mutually exclusive switches."  2>&1 |& $_TEE  -a${StdLogFile}
            ;;     
        7)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', Compilation failed."  2>&1 |& $_TEE  -a${StdLogFile}
            ;;     
        8)     
            $_ECHO -e "ERROR ${err} :\t$($_BASENAME $0): '${1}', GnuPlot - graph export error."  2>&1 |& $_TEE  -a${StdLogFile}
            ;;     
        *)
            $_ECHO -e "ERROR :\tUnknown error."  2>&1 |& $_TEE -a ${StdLogFile}
            $_ECHO ${err}  2>&1 |& $_TEE -a ${StdLogFile}
            $_ECHO $*  2>&1 |& $_TEE -a ${StdLogFile}
            ;;
    esac
    $_ECHO -e
    usage
    exit ${err}
}

# function : create directory to store data if it does not exist & validate creation
init_log_dir() 
{
    if [ ! -d ${LogDir} ] || [ ! -e ${LogDir} ] ; then 
        mkdir $LogDir        
         $_ECHO -e "WARNING :\tCreating $LogDir"  2>&1 |& $_TEE -a ${StdLogFile}
        if [[ $? -ne 0 ]] ; then 
            error 2 $LogDir
        fi
    fi  
}

# function : create log files (.txt : matrix values, .dat : timing of each computation & .log : stderr, stdout) to store values for data for each alogrithim computation
init_log_file() 
{
    localLogFile=$1
   echo -e "DEBUG : initializing ${localLogFile}"
    if  [ -e $LogDir/${localLogFile} ]  ; then
    {
        $_ECHO -e "WARNING :\tFile backup : ${localLogFile} to ${localLogFile}.bup"  |& ${localLogFile}
        mv "$LogDir/${localLogFile}" "$LogDir/${localLogFile}.bup"
    } 
    fi
    $_ECHO -e "LOG FILE :\tCreated on: ${Now} \n \tby :\t\t${USER}\n ------------------------------------------------------------------------------ \n" |& $_TEE ${localLogFile}
}

# function : add initial comments to matrix .txt and to timing .dat file
add_comments_to_log_file() {

    localLogFile=$1
    localProgramName=$2
    FileTypeDat='dat'
    FileTypeSimple='Sijk'

    $_ECHO -e "# ------------------------------------------------------------------------------------- \n# " 2>&1 |& $_TEE -a $LogDir/${localLogFile}
    $_ECHO -e "# Program : ${localProgramName} \n# \n# Log file : ${localLogFile} \n# where :\t.dat contains timing data & .txt contains matrix values \n#"
    $_ECHO -e "# Computation of Straight-forward IJK, Blocked IJK or Blocked KIJ using square NxN \n# & square bxb block (as applicable) for matrices |C| += |A| * |B| \n# "
    if [[ $localLogFile == *"$FileTypeDat"* ]] ; then
    {
    $_ECHO -e "# Time taken to compute \n#" 2>&1 |& $_TEE -a $LogDir/${localLogFile} # dat
    $_ECHO -e "# ------------------------------------------------------------------------------------- \n# " 2>&1 |& $_TEE -a $LogDir/${localLogFile}
    if [[ $localLogFile == *"$FileTypeSimple"* ]] ; then
        {
            $_ECHO "# Matrix Size \tTime/manual \tTime/manual \tTime/dgenn \n# \tSimple \tComplex \n# " 2>&1 |& $_TEE -a $LogDir/${localLogFile}
        } else {
            $_ECHO -e "# Matrix Size \tBlock Size \tTime/manual \tTime/manual \tTime/dgenn \n# \t \tSimple \t\tComplex \n# " 2>&1 |& $_TEE -a $LogDir/${localLogFile}
        }
        fi
    } else 
    {
        $_ECHO -e "# Summary of values added to each matrix - retained for later reference and validation \n#" 2>&1 |& $_TEE -a $LogDir/${localLogFile} # txt 
    }
    fi
    $_ECHO -e "# ------------------------------------------------------------------------------------- \n# " 2>&1 |& $_TEE -a $LogDir/${localLogFile}
}

# function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
     localCmd="$1"
     localOptions="$2"
     localFileMatrix="$3"
     localFileTime="$4"
     echo -e "DEBUG :\tlocalCmd\t${localCmd} \n\tl_OPTIONS\t${localOptions} \n\tlocalFileMatrix \t${localFileMatrix} \n\tlocalFileTime \t${localFileTime}"
     
     $_ECHO -e "RUNNING :\t${localCmd} ${localOptions} ${localFileMatrix} ${localFileTime}"  2>&1 |& $_TEE -a  ${StdLogFile}
     # ${localCmd}  ${localOptions} ${localFileMatrix} ${localFileTime}
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
                    BuildSimple="true"
                    BuildBlockedIJK="true"
                    BuildBlockedKIJ="true" 
                    ;;
                "-1" | "--simple")
                    BuildSimple="true"
                    ;;
                "-2" | "--ijk")
                    BuildBlockedIJK="true"
                    ;;
                "-3" | "--kij")
                    BuildBlockedKIJ="true"
                    ;;                     
                "-r" | "--random")
                    InitRandom="true"
                    ;;
                "-i" | "--increment") 
                    InitIncrement="true"
                    ;;
                "-v" | "--values")
                  	DefaultMatrixRange="true"
                      ;;                    
                "-m" | "--matrix")
                    MatrixEnabled="true"
                    if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    { 
                        let MatrixSize=$2 
                    } else 
                    { 
                    	echo "DEBUG : ${1} ${2}" 
                        error 5 "${1} ${2}" 
                    } 
                    fi                    
                    shift 
                    ;;
                "-b" | "--block")
                    BlockEnabled="true"
                    if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    { 
                        let MatrixSize=$2 
                    } else 
                    { 
                    	echo "DEBUG : ${1} ${2}" 
                        error 5 "${1} ${2}" 
                    } 
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
# TODO add check for BuildBlockedIJK="false" and BuildBlockedKIJ="false" and no block size added --> correct
# TODO check multiple -m -b -v : error not working
if [ "${DefaultMatrixRange}" == "true" ] && ( [ "${MatrixEnabled}" == "true" ] || [ "${BlockEnabled}" == "true" ] ) ; then 
{
    echo -e "DEBUG : all true \n\tDefaultMatrixRange = ${DefaultMatrixRange}\n\tMatrixEnabled=${MatrixEnabled}\n\tBlockEnabled=${BlockEnabled} "
    error 6 "<-b>, <-m> & <-v>"
}   elif  [ "${DefaultMatrixRange}" == "true" ] && ( [ "${MatrixEnabled}" == "false" ] || [ "${BlockEnabled}" == "false" ] ) ; then 
{
    echo -e "DEBUG : first true \n\tDefaultMatrixRange = ${DefaultMatrixRange}\n\tMatrixEnabled=${MatrixEnabled}\n\tBlockEnabled=${BlockEnabled}"
     declare -a NXArray=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
    declare -a NBArray=( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 )
    MatrixEnabled="false"
    BlockEnabled="false"
    echo ${#NXArray[@]} 
    echo ${#NBArray[@]} 
}   else 
    {
    echo -e "DEBUG : matrix & block only \n\tDefaultMatrixRange = ${DefaultMatrixRange}\n\tMatrixEnabled=${MatrixEnabled}\n\tBlockEnabled=${BlockEnabled} "
    if [ "${MatrixEnabled}" == "true" ] || [ "${BlockEnabled}" == "true" ] ; then
        if [ "${MatrixSize}" == "" ] ; then 
            error 4 "matrix size"
        fi
	   # Validate matrix range
        if [ ${MatrixSize} -le 0 ] || [ ${MatrixSize} -gt ${MaxMatrixSize} ] ; then
            $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is invalid, now set to default of : $MaxMatrixSize"  2>&1 |& $_TEE -a  ${StdLogFile}
            let MatrixSize=$MaxMatrixSize
            MatrixEnabled="true"
        fi        
        if [ "${BlockSize}" == "" ] ; then 
            error 4 "block size"
        fi
	   # Validate block range
        if [ ${BlockSize} -le 0 ] || [ ${BlockSize} -gt ${MaxBlockSize} ] ; then
            $_ECHO -e "WARNING :\t$($_BASENAME $0): block size <nb> is invalid, now set to default of : $MaxBlockSize"  2>&1 |& $_TEE -a  ${StdLogFile}
            let BlockSize=$MaxBlockSize
            BlockEnabled="true"
        fi
        # ensure both are enabled if one is
        if [ "${MatrixEnabled}" == "false" ] ; then 
            $_ECHO -e "WARNING :\t$($_BASENAME $0): matrix size <nx> is not enabled, now enabled and set to default of : $MaxMatrixSize"  2>&1 |& $_TEE -a  ${StdLogFile}
            let MatrixSize=$MaxMatrixSize
            MatrixEnabled="true"
        fi
        if [ "${BlockEnabled}" == "false" ] ; then 
            $_ECHO -e "WARNING :\t$($_BASENAME $0): block size <nb> is not enabled, now enabled and set to default of : $MaxBlockSize"  2>&1 |& $_TEE -a  ${StdLogFile}
            let BlockSize=$MaxBlockSize
            BlockEnabled="true"
        fi
        # Validate matrix / block remainder
    	if [ $(( ${MatrixSize} % ${BlockSize} )) -ne 0 ]; then
	        $_ECHO -e "WARNING :\t$($_BASENAME $0): block size needs to be an even multiple of matrix size. Using defaults."  2>&1 |& $_TEE -a  ${StdLogFile}
	        let BlockSize=$MaxBlockSize
	        let MatrixSize=$MaxMatrixSize
    	fi
        declare -a NXArray=( $MatrixSize )
        declare -a NBArray=( $BlockSize )
    fi
}
fi

# build up commands to run
AlgorithmOptions=""
if [ "$InitRandom" == "true" ] && [ "$InitIncrement" == "false" ] ; then 
    AlgorithmOptions="-r"
fi
if [ "$InitRandom" == "false" ] && [ "$InitIncrement" == "true" ] ; then 
    AlgorithmOptions="-i"
fi
if [ "$InitRandom" == "true" ] && [ "$InitIncrement" == "true" ] ; then 
        error 6 "<-i> and <-r>"    
fi

# execute algorithms

init_log_dir
init_log_file  $StdLogFile

if [ "$BuildSimple" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then
    Algorithm_SIMPLE="./Programs/A1-Sijk-1D"    
    LogType_SIMPLE="A1.Sijk-"
    MAT_FILE_SIMPLE="${LogPrefix}${LogType_SIMPLE}${Now}"
    DAT_FILE_SIMPLE="${MAT_FILE_SIMPLE}"
    STD_FILE_SIMPLE="${MAT_FILE_SIMPLE}"
    STD_FILE_SIMPLE="${STD_FILE_SIMPLE}${STD_SUFFIX}"
    
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
        error 5 "matrix size and block size"
    else        
        ExecuteOptions=""
        for (( i = 0 ; i < ${#NXArray[@]} ; i++ )) do
            MAT_FILE_SIMPLE_VALUES="${MAT_FILE_SIMPLE}-values-$i${TxtSuffix}"
            DAT_FILE_SIMPLE_TIMING="${DAT_FILE_SIMPLE}-timing-$i${DatSuffix}"
            init_log_file $MAT_FILE_SIMPLE_VALUES 
            init_log_file $DAT_FILE_SIMPLE_TIMING
            init_log_file $DAT_FILE_SIMPLE_TIMING

            add_comments_to_log_file "${MAT_FILE_SIMPLE_VALUES}" "${Algorithm_SIMPLE}"
            add_comments_to_log_file "${DAT_FILE_SIMPLE_TIMING}" "${Algorithm_SIMPLE}"
            ExecuteOptions="${AlgorithmOptions} ${NXArray[$i]}"
            algorithm_execute "${Algorithm_SIMPLE}" "${ExecuteOptions}" "${MAT_FILE_SIMPLE_VALUES}" "${DAT_FILE_SIMPLE_TIMING}"
            MAT_FILE_SIMPLE_VALUES="${MAT_FILE_SIMPLE}"
            DAT_FILE_SIMPLE_TIMING="${DAT_FILE_SIMPLE}"
        done
    fi
fi

if [ "$BuildBlockedIJK" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then 
    Algorithm_BLOCKED_IJK="./Programs/A1-Bijk-1D"
    LogType_BIJK="A1-Bijk-"
    MAT_FILE_BIJK="${LogPrefix}${LogType_BIJK}${Now}"
    DAT_FILE_BIJK="${MAT_FILE_BIJK}"
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
        error 5 "matrix size and block size"
    else 
        ExecuteOptions=""
        for (( i = 0 ; i < ${#NXArray[@]} ; i++ )) do
            MAT_FILE_BIJK_VALUES="${MAT_FILE_BIJK}-values-$i${TxtSuffix}"
            DAT_FILE_BIJK_TIMING="${DAT_FILE_BIJK}-timing-$i${DatSuffix}"
            init_log_file $MAT_FILE_BIJK_VALUES 
            init_log_file $DAT_FILE_BIJK_TIMING
            add_comments_to_log_file "${MAT_FILE_BIJK_VALUES}" "${Algorithm_BIJK}"
            add_comments_to_log_file "${DAT_FILE_BIJK_TIMING}" "${Algorithm_BIJK}"
            ExecuteOptions="${AlgorithmOptions} ${NXArray[$i]} ${NBArray[$i]}"
            algorithm_execute "${Algorithm_BLOCKED_IJK}" "${ExecuteOptions}" "${MAT_FILE_BIJK_VALUES}" "${DAT_FILE_BIJK_TIMING}"
            MAT_FILE_BIJK_VALUES="${MAT_FILE_BIJK}"
            DAT_FILE_BIJK_TIMING="${DAT_FILE_BIJK}"
        done
    fi
fi

if [ "$BuildBlockedKIJ" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then 
    Algorithm_BLOCKED_KIJ="./Programs/A1-Bkij-1D"
    LogType_BKIJ="A1-Bkij-"
    MAT_FILE_BKIJ="${LogPrefix}${LogType_BKIJ}${Now}"
    DAT_FILE_BKIJ="${MAT_FILE_BKIJ}"
    if [[ ${#NXArray[*]} -ne ${#NBArray[*]} ]] ; then 
        error 5 "matrix size and block size"
    else 
        ExecuteOptions=""
        for (( i = 0 ; i < ${#NXArray[@]} ; i++ )) do
            MAT_FILE_BKIJ_VALUES="${MAT_FILE_BKIJ}-values-$i${TxtSuffix}"
            DAT_FILE_BKIJ_TIMING="${DAT_FILE_BKIJ}-timing-$i${DatSuffix}"
            init_log_file $MAT_FILE_BKIJ_VALUES 
            init_log_file $DAT_FILE_BKIJ_TIMING
            add_comments_to_log_file "${MAT_FILE_BKIJ_VALUES}" "${Algorithm_BKIJ}" 
            add_comments_to_log_file "${DAT_FILE_BKIJ_TIMING}" "${Algorithm_BKIJ}" 
            ExecuteOptions="${AlgorithmOptions} ${NXArray[$i]} ${NBArray[$i]}"
            algorithm_execute "${Algorithm_BLOCKED_KIJ}" "${ExecuteOptions}" "${MAT_FILE_BKIJ_VALUES}" "${DAT_FILE_BKIJ_TIMING}"
            MAT_FILE_BKIJ_VALUES="${MAT_FILE_BKIJ}"
            DAT_FILE_BKIJ_TIMING="${DAT_FILE_BKIJ}"
        done
    fi
fi

pause
exit 0
