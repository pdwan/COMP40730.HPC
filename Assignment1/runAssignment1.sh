
# ##################################################################################
# 
# DESC:         Script to multiply two nxn matrices using three algorithms.
# AUTHOR :   Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# DATE :     30-June-2014
# ASSIGNMENT : 1
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
_TEE="tee"
NOW=$(date +"%Y%m%d-%H%M%S")
LOG_DIR="logDir"
LOG_PREFIX="pdwan-"
TXT_SUFFIX=".txt"
DAT_SUFFIX=".dat"
STD_LOG_FILE=$(_BASENAME $0 +"-${NOW}.log")
BUILD_SIMPLE="false" 
BUILD_BLOCKED_IJK="false"
BUILD_BLOCKED_KIJ="false" 
INIT_RANDOM="false"
INIT_INCREMENT="false"
MATRIX_ENABLED="false"
BLOCK_ENABLED="false"
DEFAULT_MATRIX_RANGE="false"
let MATRIX_SIZE=0
let BLOCK_SIZE=0 
let MAX_MATRIX_SIZE=1000 
let MAX_BLOCK_SIZE=50

# ##################################################################################

#   function : pause
pause () {
    $_ECHO -e "\n"
    read -p "Press <any key> to continue ... " nothing
}

#   function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE :\t./$($_BASENAME $0)  \ \n\t\t\t-a|--all -1|--simple -2|--ijk -3|--kij -r|--random -i|--increment \ \n\t\t\t-m|--matrix<n> -b|--block <b> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO : \tCalculate |C| = |A| x |B| using one or three algoritms : Straight-forward IJK, Blocked IJK and Blocked KIJ. \n"   
    $_ECHO -e "where:\t-a | --all \tCalculate data for all algorithms via separate .c programs to multiply |A|x|B| -> |C|"
    $_ECHO -e "\t\t\tStraightforward non-blocked ijk algorithm :\tA1-Sijk-1D.c \n\t\t\tBlocked ijk algorithm using square bxb blocks :\tA1-Bijk-1D.c \n\t\t\tBlocked kij algorithm using square bxb blocks :\tA1-Bkij-1D.c"
    $_ECHO -e "\t-1|--simple \tCalculate data for only the algorithm \n\t\t\tStraightforward non-blocked ijk algorithm : \tA1-Sijk-1D.c "
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

#   function : error message and then usage
error() 
{
    ERR=$1
    
    shift
    case ${ERR} in
        1)
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): Unknown parameter <${1}>."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;
        2)
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): Error creating directory <${1}>."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;
        3)
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): Error creating file <${1}>."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;     
        4)     
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): Missing parameter for <${1}>."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;     
        5)     
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): <${1}>, Values entered are not valid or not a number."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;     
        6)     
            $_ECHO -e "ERROR ${ERR} : \t$($_BASENAME $0): <${1}>, Mutually exclusive switches."  2>&1 |& $_TEE  -a${STD_LOG_FILE}
            ;;     
        *)
            $_ECHO -e "ERROR : \tUnknown error."  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            $_ECHO ${ERR}  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            $_ECHO $*  2>&1 |& $_TEE -a ${STD_LOG_FILE}
            ;;
    esac
    $_ECHO -e
    usage
    exit ${ERR}
}

#   function : create directory to store data if it does not exist & validate creation
init_log_dir() 
{
    if [ ! -d ${LOG_DIR} ] || [ ! -e ${LOG_DIR} ] ; then 
        mkdir $LOG_DIR        
         $_ECHO -e "WARNING : \tCreating $LOG_DIR"  2>&1 |& $_TEE -a ${STD_LOG_FILE}
        if [[ $? -ne 0 ]] ; then 
            error 2 $LOG_DIR
        fi
    fi  
}

#   function : create log files (.txt : matrix values, .dat : timing of each computation & .log : stderr, stdout) to store values for data for each alogrithim computation
init_log_file() 
{
    l_LOG_FILE=$1
   echo -e "DEBUG : initializing ${l_LOG_FILE}"
    if  [ -e $LOG_DIR/${l_LOG_FILE} ]  ; then
    {
        $_ECHO -e "WARNING : \tFile backup : ${l_LOG_FILE} to ${l_LOG_FILE}.bup"  |& ${l_LOG_FILE}
        mv "$LOG_DIR/${l_LOG_FILE}" "$LOG_DIR/${l_LOG_FILE}.bup"
    } 
    fi
    $_ECHO -e "LOG FILE : \tCreated on: ${NOW} \n \tby : \t\t${USER}\n ------------------------------------------------------------------------------ \n" |& $_TEE ${l_LOG_FILE}
}

#   function : add initial comments to matrix .txt and to timing .dat file
add_comments_to_log_file() {

    l_LOG_FILE=$1
    l_PROG_NAME=$2
    file_type_dat='dat'
    file_type_simple='Sijk'

    $_ECHO -e "# Program : ${l_PROG_NAME} \n# \n# Log file : ${l_LOG_FILE} \n# where : \t.dat contains timing data & .txt contains matrix values \n#"
    $_ECHO -e "# Computation of Straight-forward IJK, Blocked IJK or Blocked KIJ using square NxN \n# & square bxb block (as applicable) for matrices |C| += |A| * |B| \n# "
    if [[ $l_LOG_FILE == *"$file_type_dat"* ]] ; then
    {
    $_ECHO -e "# Time taken to compute \n#" 2>&1 |& $_TEE -a $LOG_DIR/${l_LOG_FILE} # dat
    if [[ $l_LOG_FILE == *"$file_type_simple"* ]] ; then
        {
            $_ECHO "# Matrix Size \tTime/manual \tTime/manual \tTime/dgenn \n# \tSimple \tComplex \n# " 2>&1 |& $_TEE -a $LOG_DIR/${l_LOG_FILE}
        } else {
            $_ECHO -e "# Matrix Size \tBlock Size \tTime/manual \tTime/manual \tTime/dgenn \n# \t \tSimple \t\tComplex \n# " 2>&1 |& $_TEE -a $LOG_DIR/${l_LOG_FILE}
        }
        fi
    } else 
    {
        $_ECHO -e "# Summary of values added to each matrix - retained for later reference and validation \n#" 2>&1 |& $_TEE -a $LOG_DIR/${l_LOG_FILE} # txt 
    }
    fi
    $_ECHO -e "# ------------------------------------------------------------------------------------- \n# " 2>&1 |& $_TEE -a $LOG_DIR/${l_LOG_FILE}
}

#   function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
     l_CMD="$1"
     l_OPTIONS="$2"
     l_FILE_MATRIX="$3"
     l_FILE_TIME="$4"
     echo -e "DEBUG : \tl_CMD\t${l_CMD} \n\tl_OPTIONS\t${l_OPTIONS} \n\tl_FILE_MATRIX \t${l_FILE_MATRIX} \n\tl_FILE_TIME \t${l_FILE_TIME}"
     
     $_ECHO -e "RUNNING :\t${l_CMD} ${l_OPTIONS} ${l_FILE_MATRIX} ${l_FILE_TIME}"  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
     # ${l_CMD}  ${l_OPTIONS} ${l_FILE_MATRIX} ${l_FILE_TIME}
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
                    BUILD_SIMPLE="true"
                    BUILD_BLOCKED_IJK="true"
                    BUILD_BLOCKED_KIJ="true" 
                    ;;
                "-1" | "--simple")
                    BUILD_SIMPLE="true"
                    ;;
                "-2" | "--ijk")
                    BUILD_BLOCKED_IJK="true"
                    ;;
                "-3" | "--kij")
                    BUILD_BLOCKED_KIJ="true"
                    ;;                     
                "-r" | "--random")
                    INIT_RANDOM="true"
                    ;;
                "-i" | "--increment") 
                    INIT_INCREMENT="true"
                    ;;
                "-v" | "--values")
                     DEFAULT_MATRIX_RANGE="true"
                      ;;                    
                "-m" | "--matrix")
                    MATRIX_ENABLED="true"
                    if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    { 
                        let MATRIX_SIZE=$2 
                    } else 
                    { 
                    echo "DEBUG : ${1} ${2}" 
                        error 5 "${1} ${2}" 
                    } 
                    fi                    
                    shift 
                    ;;
                "-b" | "--block")
                    BLOCK_ENABLED="true"
                    if [ "${2}" -eq "${2}" ] 2>/dev/null ; then
                    { 
                        let MATRIX_SIZE=$2 
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
# TODO add check for BUILD_BLOCKED_IJK="true" and BUILD_BLOCKED_KIJ="true" and no block size added --> correct
# TODO check multiple -m -b -v : error not working
if [ "${DEFAULT_MATRIX_RANGE}" == "true" ] && ( [ "${MATRIX_ENABLED}" == "true" ] || [ "${BLOCK_ENABLED}" == "true" ] ) ; then 
{
    echo -e "DEBUG : all true \n\tDEFAULT_MATRIX_RANGE = ${DEFAULT_MATRIX_RANGE}\n\tMATRIX_ENABLED=${MATRIX_ENABLED}\n\tBLOCK_ENABLED=${BLOCK_ENABLED} "
    error 6 "<-b>, <-m> & <-v>"
}   elif  [ "${DEFAULT_MATRIX_RANGE}" == "true" ] && ( [ "${MATRIX_ENABLED}" == "false" ] || [ "${BLOCK_ENABLED}" == "false" ] ) ; then 
{
    echo -e "DEBUG : first true \n\tDEFAULT_MATRIX_RANGE = ${DEFAULT_MATRIX_RANGE}\n\tMATRIX_ENABLED=${MATRIX_ENABLED}\n\tBLOCK_ENABLED=${BLOCK_ENABLED}"
     declare -a NX_ARRAY=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 )
    declare -a NB_ARRAY=( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 )
    MATRIX_ENABLED="false"
    BLOCK_ENABLED="false"
    echo ${#NX_ARRAY[@]} 
    echo ${#NB_ARRAY[@]} 
}   else 
    {
    echo -e "DEBUG : matrix & block only  \n\tDEFAULT_MATRIX_RANGE = ${DEFAULT_MATRIX_RANGE}\n\tMATRIX_ENABLED=${MATRIX_ENABLED}\n\tBLOCK_ENABLED=${BLOCK_ENABLED} "
    if [ "${MATRIX_ENABLED}" == "true" ] || [ "${BLOCK_ENABLED}" == "true" ] ; then
        if [ "${MATRIX_SIZE}" == "" ] ; then 
            error 4 "matrix size"
        fi
	   # Validate matrix range
        if [ ${MATRIX_SIZE} -le 0 ] || [ ${MATRIX_SIZE} -gt ${MAX_MATRIX_SIZE} ] ; then
            $_ECHO -e "WARNING : \t$($_BASENAME $0): matrix size <nx> is invalid, now set to default of : $MAX_MATRIX_SIZE"  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
            let MATRIX_SIZE=$MAX_MATRIX_SIZE
            MATRIX_ENABLED="true"
        fi        
        if [ "${BLOCK_SIZE}" == "" ] ; then 
            error 4 "block size"
        fi
	   # Validate block range
        if [ ${BLOCK_SIZE} -le 0 ] || [ ${BLOCK_SIZE} -gt ${MAX_BLOCK_SIZE} ] ; then
            $_ECHO -e "WARNING : \t$($_BASENAME $0): block size <nb> is invalid, now set to default of : $MAX_BLOCK_SIZE"  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
            let BLOCK_SIZE=$MAX_BLOCK_SIZE
            BLOCK_ENABLED="true"
        fi
        # ensure both are enabled if one is
        if [ "${MATRIX_ENABLED}" == "false" ] ; then 
            $_ECHO -e "WARNING : \t$($_BASENAME $0): matrix size <nx> is not enabled, now enabled and set to default of : $MAX_MATRIX_SIZE"  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
            let MATRIX_SIZE=$MAX_MATRIX_SIZE
            MATRIX_ENABLED="true"
        fi
        if [ "${BLOCK_ENABLED}" == "false" ] ; then 
            $_ECHO -e "WARNING : \t$($_BASENAME $0): block size <nb> is not enabled, now enabled and set to default of : $MAX_BLOCK_SIZE"  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
            let BLOCK_SIZE=$MAX_BLOCK_SIZE
            BLOCK_ENABLED="true"
        fi
        # Validate matrix / block remainder
    	if [ $(( ${MATRIX_SIZE} % ${BLOCK_SIZE} )) -ne 0 ]; then
	        $_ECHO -e "WARNING : \t$($_BASENAME $0): block size needs to be an even multiple of matrix size. Using defaults."  2>&1 |& $_TEE -a  ${STD_LOG_FILE}
	        let BLOCK_SIZE=$MAX_BLOCK_SIZE
	        let MATRIX_SIZE=$MAX_MATRIX_SIZE
    	fi
        declare -a NX_ARRAY=( $MATRIX_SIZE )
        declare -a NB_ARRAY=( $BLOCK_SIZE )
    fi
}
fi

# build up commands to run
ALGORITHIM_OPTIONS=""
if [ "$INIT_RANDOM" == "true" ] && [ "$INIT_INCREMENT" == "false" ] ; then 
    ALGORITHIM_OPTIONS="-r"
fi
if [ "$INIT_RANDOM" == "false" ] && [ "$INIT_INCREMENT" == "true" ] ; then 
    ALGORITHIM_OPTIONS="-i"
fi
if [ "$INIT_RANDOM" == "true" ] && [ "$INIT_INCREMENT" == "true" ] ; then 
        error 6 "<-i> and <-r>"    
fi

# execute algorithms

init_log_dir
init_log_file  $STD_LOG_FILE

if [ "$BUILD_SIMPLE" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then
    ALGORITHIM_SIMPLE="./Programs/A1-Sijk-1D"    
    LOG_TYPE_SIMPLE="A1.Sijk-"
    MAT_FILE_SIMPLE="${LOG_PREFIX}${LOG_TYPE}${NOW}"
    DAT_FILE_SIMPLE="${MAT_FILE_SIMPLE}"
    STD_FILE_SIMPLE="${MAT_FILE_SIMPLE}"
    STD_FILE_SIMPLE="${STD_FILE_SIMPLE}${STD_SUFFIX}"
    
    if [[ ${#NX_ARRAY[*]} -ne ${#NB_ARRAY[*]} ]] ; then 
        error 5 "matrix size and block size"
    else        
        EXECUTE_OPTIONS=""
        for (( i = 0 ; i < ${#NX_ARRAY[@]} ; i++ )) do
            MAT_FILE_SIMPLE_VALUES="${MAT_FILE_SIMPLE}-values-$i${TXT_SUFFIX}"
            DAT_FILE_SIMPLE_TIMING="${DAT_FILE_SIMPLE}-timing-$i${DAT_SUFFIX}"
            init_log_file $MAT_FILE_SIMPLE_VALUES 
            init_log_file $DAT_FILE_SIMPLE_TIMING
            init_log_file $DAT_FILE_SIMPLE_TIMING

            add_comments_to_log_file "${MAT_FILE_SIMPLE_VALUES}" "${ALGORITHIM_SIMPLE}"
            add_comments_to_log_file "${DAT_FILE_SIMPLE_TIMING}" "${ALGORITHIM_SIMPLE}"
            EXECUTE_OPTIONS="${ALGORITHIM_OPTIONS} ${NX_ARRAY[$i]}"
            algorithm_execute "${ALGORITHIM_SIMPLE}" "${EXECUTE_OPTIONS}" "${MAT_FILE_SIMPLE_VALUES}" "${DAT_FILE_SIMPLE_TIMING}"
            MAT_FILE_SIMPLE_VALUES="${MAT_FILE_SIMPLE}"
            DAT_FILE_SIMPLE_TIMING="${DAT_FILE_SIMPLE}"
        done
    fi
fi

if [ "$BUILD_BLOCKED_IJK" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then 
    ALGORITHIM_BLOCKED_IJK="./Programs/A1-Bijk-1D"
    LOG_TYPE="A1-Bijk-"
    MAT_FILE_BIJK="${LOG_PREFIX}${LOG_TYPE}${NOW}"
    DAT_FILE_BIJK="${MAT_FILE_BIJK}"
    if [[ ${#NX_ARRAY[*]} -ne ${#NB_ARRAY[*]} ]] ; then 
        error 5 "matrix size and block size"
    else 
        EXECUTE_OPTIONS=""
        for (( i = 0 ; i < ${#NX_ARRAY[@]} ; i++ )) do
            MAT_FILE_BIJK_VALUES="${MAT_FILE_BIJK}-values-$i${TXT_SUFFIX}"
            DAT_FILE_BIJK_TIMING="${DAT_FILE_BIJK}-timing-$i${DAT_SUFFIX}"
            init_log_file $MAT_FILE_BIJK_VALUES 
            init_log_file $DAT_FILE_BIJK_TIMING
            add_comments_to_log_file "${MAT_FILE_BIJK_VALUES}" "${ALGORITHIM_BIJK}"
            add_comments_to_log_file "${DAT_FILE_BIJK_TIMING}" "${ALGORITHIM_BIJK}"
            EXECUTE_OPTIONS="${ALGORITHIM_OPTIONS} ${NX_ARRAY[$i]} ${NB_ARRAY[$i]}"
            algorithm_execute "${ALGORITHIM_BLOCKED_IJK}" "${EXECUTE_OPTIONS}" "${MAT_FILE_BIJK_VALUES}" "${DAT_FILE_BIJK_TIMING}"
            MAT_FILE_BIJK_VALUES="${MAT_FILE_BIJK}"
            DAT_FILE_BIJK_TIMING="${DAT_FILE_BIJK}"
        done
    fi
fi

if [ "$BUILD_BLOCKED_KIJ" == "true" ]  || [ "$BUILD_ALL" == "true" ] ; then 
    ALGORITHIM_BLOCKED_KIJ="./Programs/A1-Bkij-1D"
    LOG_TYPE="A1-Bkij-"
    MAT_FILE_BKIJ="${LOG_PREFIX}${LOG_TYPE}${NOW}"
    DAT_FILE_BKIJ="${MAT_FILE_BKIJ}"
    if [[ ${#NX_ARRAY[*]} -ne ${#NB_ARRAY[*]} ]] ; then 
        error 5 "matrix size and block size"
    else 
        EXECUTE_OPTIONS=""
        for (( i = 0 ; i < ${#NX_ARRAY[@]} ; i++ )) do
            MAT_FILE_BKIJ_VALUES="${MAT_FILE_BKIJ}-values-$i${TXT_SUFFIX}"
            DAT_FILE_BKIJ_TIMING="${DAT_FILE_BKIJ}-timing-$i${DAT_SUFFIX}"
            init_log_file $MAT_FILE_BKIJ_VALUES 
            init_log_file $DAT_FILE_BKIJ_TIMING
            add_comments_to_log_file "${MAT_FILE_BKIJ_VALUES}" "${ALGORITHIM_BKIJ}" 
            add_comments_to_log_file "${DAT_FILE_BKIJ_TIMING}" "${ALGORITHIM_BKIJ}" 
            EXECUTE_OPTIONS="${ALGORITHIM_OPTIONS} ${NX_ARRAY[$i]} ${NB_ARRAY[$i]}"
            algorithm_execute "${ALGORITHIM_BLOCKED_KIJ}" "${EXECUTE_OPTIONS}" "${MAT_FILE_BKIJ_VALUES}" "${DAT_FILE_BKIJ_TIMING}"
            MAT_FILE_BKIJ_VALUES="${MAT_FILE_BKIJ}"
            DAT_FILE_BKIJ_TIMING="${DAT_FILE_BKIJ}"
        done
    fi
fi

pause
exit 0
