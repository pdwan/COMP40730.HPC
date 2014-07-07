
# ##################################################################################
# 
# DESC:         Script to calculate infinity norm of |C| where |C| =+ |A| * |B|
# AUTHOR :      Paula Dwan (paula.dwan@ericsson.com | paula.dwan@gmail.com)
# DATE :        30-June-2014
# ASSIGNMENT :  2
#
# ##################################################################################

# constants and variables and executables
_ECHO="echo"
_BASENAME="basename"
NOW=$(date +"%Y%m%d-%H%M%S")
LOG_DIR="logDir"
LOG_PREFIX="pdwan-"
TXT_SUFFIX=".txt"
DAT_SUFFIX=".dat"
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

#  function pause
pause () {
    read -p "Press <any key> to continue ... " nothing
}

# function : usage instructions
usage() 
{
    $_ECHO -e "\nUSAGE : ./$($_BASENAME $0) -a|--all -1|--simple -2|--ijk -3|--kij -r|--random -i|--increment \ \n\t\t\t  -m|--matrix<n> -b|--block <b> -v|--values -?|-h|--help \n"
    $_ECHO -e "TO : \tCalculate infinity norm using pthreads and single calculation using dgemm ATLAS / BLAS. \n"
    $_ECHO -e "\twhere:"
    $_ECHO -e "\t-r|--random :     initialization A| & |B| with random numbers and |C| with '0' "
    $_ECHO -e "\t-i|--increment :  initialize |A| & |B| incrementally with <row> value and |C| with '0' "
    $_ECHO -e "\t\t\t  '-i|--increment' & '-r|--random' are mutually exclusive. \n"
    $_ECHO -e "\t-m|--matrix <n> : matrix dimension, set to maximum of '1,000' if invalid or not provided"
    $_ECHO -e "\t-b|--block <b> :  block size, with '<nb> < <nx>' and '<nx> % <nb> = 0', \n\t\t\t  defaults to '50' if invalid or not provided and <nx> will be reset to 1,000." 
    $_ECHO -e "\t\t\t  Mutually exclusive with '-v|--values' \n"    
    $_ECHO -e "\t-v|--values :     use predefined range of valid values for <nx> and <nb> as follows :"
    $_ECHO -e "\t\t\t  <nx> = { 50, 50, 50, 100, 100, 100, 500, 500, 500, 1000, 1000, 1000 } "
    $_ECHO -e "\t\t\t  <nb> = { 2, 5, 10, 5, 10, 20, 10, 20, 50, 10, 50, 100 } "
    $_ECHO -e "\t\t\t  Mutually exclusive with '-m|--matrix <n>' and '-b|--block <b>' \n"
    $_ECHO -e "\t-?|-h|--help :    usage"
    $_ECHO ""
}

# function : error message and then usage
error() 
{
    ERR=$1
    shift
    case ${ERR} in
        1)
            $_ECHO -e "ERROR : \t$($_BASENAME $0): Unknown parameter $1."
            ;;
        2)
            $_ECHO -e "ERROR : \t$($_BASENAME $0): Error creating directory $1."
            ;;
        3)
            $_ECHO -e "ERROR : \t$($_BASENAME $0): Error creating file $1."
            ;;     
        4)     
            $_ECHO -e "ERROR : \t$($_BASENAME $0): Missing parameter for $1."
            ;;     
        5)     
            $_ECHO -e "ERROR : \t$($_BASENAME $0): $1, Values entered are not valid."
            ;;     
        6)     
            $_ECHO -e "ERROR : \t$($_BASENAME $0): $1, Mutually exclusive switches."
            ;;     
        7)     
            $_ECHO -e "ERROR : \t$($_BASENAME $0): $1, Not a number."
            ;;     
        *)
            $_ECHO -e "ERROR : \tUnknown error."
            $_ECHO ${ERR}
            $_ECHO $*
            ;;
    esac
    $_ECHO -e
    usage
    exit ${ERR}
}

# function : create directory to store data if it does not exist & validate creation
init_log_dir() 
{
    if [ ! -d ${LOG_DIR} ] || [ ! -e ${LOG_DIR} ] ; then 
        mkdir $LOG_DIR        
         $_ECHO -e "WARNING : \tcreating $LOG_DIR" 
        if [[ $? -ne 0 ]] ; then 
            error 2 $LOG_DIR
        fi
    fi  
}

# function : create log file to store values for data for each alogrithim computation
init_log_file() 
{
    l_LOG_FILE=$1
    if  [ -e $LOG_DIR/${l_LOG_FILE} ]  ;
    then
        $_ECHO -e "WARNING : \tFile backup  :  ${l_LOG_FILE} to ${l_LOG_FILE}.bup"
        mv "$LOG_DIR/${l_LOG_FILE}" "$LOG_DIR/${l_LOG_FILE}.bup"
    fi
}

# function : add initial comments to matrix .txt and to timing .dat file
add_comments_to_log_file() {

    l_LOG_FILE=$1
    l_PROG_NAME=$2
    file_type_dat='dat'
    
    $_ECHO -e "# Program : ${l_PROG_NAME} \n# Log file  : ${l_LOG_FILE} \n# where :  \n#     .dat contains timing data & .txt contains matrix values \n#" > $LOG_DIR/${l_LOG_FILE}
    $_ECHO -e "# Computation of Straight-forward IJK, Blocked IJK or Blocked KIJ using square NxN \n# & square bxb block (as applicable) for matrices |C| += |A| * |B| \n# " >> $LOG_DIR/${l_LOG_FILE}
    if [[ $l_LOG_FILE == *"$file_type_dat"* ]] ; then
    {
        $_ECHO -e "#Matrix Size \\tBlock Size \\tTime/manual \\tTime/manual \\tTime/dgenn \n# \\t \\tSimple \tComplex \n#" >> $LOG_DIR/${l_LOG_FILE}
    } else 
    {
        $_ECHO -e "# Summary of values added to each matrix - retained for later reference and validation \n# \n#" >> $LOG_DIR/${l_LOG_FILE} # txt 
    }
    fi
    $_ECHO -e "# ------------------------------------------------------------------------------------- \n# " >> $LOG_DIR/${l_LOG_FILE}
}

# function : execute each algorithm in turn wih specified parameters / options
algorithm_execute() 
{
     l_CMD="$1"
     l_OPTIONS="$2"
     l_FILE_MATRIX="$3"
     l_FILE_TIME="$4"
     
     $_ECHO -e "RUNNING : \t${l_CMD} ${l_OPTIONS} ${l_FILE_MATRIX} ${l_FILE_TIME}"
     # ${l_CMD}  ${l_OPTIONS} ${l_FILE_MATRIX} ${l_FILE_TIME}
}

# ##################################################################################

# clear

# Process parameters
if [[ $# -eq 0 ]]; then
    usage
    exit
else 
    while [ "$1" != "" ]; do
        case ${1,,} in
                "-r" | "--random")
                    INIT_RANDOM="true"
                    ;;
                "-i" | "--increment") 
                    INIT_INCREMENT="true"
                    ;;
                "-m" | "--matrix")
                    MATRIX_ENABLED="true"
                    # if ! [ $2 =~ $'^[0-9]+$' ] ; then
                    if [ "$2" -eq "$2" ] 2>/dev/null ; then
                    { 
                        let MATRIX_SIZE=$2 
                    } else 
                    { 
                        error 7 "<${2}>" 
                    } 
                    fi                    
                    shift 
                    ;;
                "-b" | "--block")
                    BLOCK_ENABLED="true"
                    if ! [[ $2 =~ $'^[0-9]+$' ]] ; then
                    { 
                        let BLOCK_SIZE=$2
                    } else 
                    { 
                        error 7 "<${2}>" 
                    } 
                    fi                    
                    shift 
                    ;;                     
                "-v" | "--values")
                     DEFAULT_MATRIX_RANGE="true"
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
if [ "${DEFAULT_MATRIX_RANGE}" == "true" ] ; then 
    declare -a NX_ARRAY=( 50 50 50 100 100 100 500 500 500 500 1000 1000 1000 1000 ) ;
    declare -a NB_ARRAY=( 2 5 10 5 10 20 5 10 20 50 5 10 50 100 ) ;
    MATRIX_ENABLED="false"
    BLOCK_ENABLED="false"
else
    if [ "${MATRIX_ENABLED}" == "true" ] || [ "${BLOCK_ENABLED}" == "true" ] ; then
        if [ "${MATRIX_SIZE}" == "" ] ; then 
            error 4 "matrix size"
        fi
	   # Validate matrix range
        if [ ${MATRIX_SIZE} -le 0 ] || [ ${MATRIX_SIZE} -gt ${MAX_MATRIX_SIZE} ] ; then
            $_ECHO -e "WARNING : \t$($_BASENAME $0): matrix size <nx> is invalid, now set to default of : $MAX_MATRIX_SIZE"
            let MATRIX_SIZE=$MAX_MATRIX_SIZE
            MATRIX_ENABLED="true"
        fi        
        if [ "${BLOCK_SIZE}" == "" ] ; then 
            error 4 "block size"
        fi
	   # Validate block range
        if [ ${BLOCK_SIZE} -le 0 ] || [ ${BLOCK_SIZE} -gt ${MAX_BLOCK_SIZE} ] ; then
            $_ECHO -e "WARNING : \t$($_BASENAME $0): block size <nb> is invalid, now set to default of : $MAX_BLOCK_SIZE"
            let BLOCK_SIZE=$MAX_BLOCK_SIZE
            BLOCK_ENABLED="true"
        fi
        # ensure both are enabled if one is
        if [ "${MATRIX_ENABLED}" == "false" ] ; then 
            $_ECHO -e "WARNING : \t$($_BASENAME $0): matrix size <nx> is not enabled, now enabled and set to default of : $MAX_MATRIX_SIZE"
            let MATRIX_SIZE=$MAX_MATRIX_SIZE
            MATRIX_ENABLED="true"
        fi
        if [ "${BLOCK_ENABLED}" == "false" ] ; then 
            $_ECHO -e "WARNING : \t$($_BASENAME $0): block size <nb> is not enabled, now enabled and set to default of : $MAX_BLOCK_SIZE"
            let BLOCK_SIZE=$MAX_BLOCK_SIZE
            BLOCK_ENABLED="true"
        fi
        # Validate matrix / block remainder
    	if [ $(( ${MATRIX_SIZE} % ${BLOCK_SIZE} )) -ne 0 ]; then
	        $_ECHO -e "WARNING : \t$($_BASENAME $0): block size needs to be an even multiple of matrix size. Using defaults."
	        let BLOCK_SIZE=$MAX_BLOCK_SIZE
	        let MATRIX_SIZE=$MAX_MATRIX_SIZE
    	fi
        declare -a NX_ARRAY=( $MATRIX_SIZE )
        declare -a NB_ARRAY=( $BLOCK_SIZE )
    fi
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

PROG_NAME="./Programs/A1-pthreads-1D"    
LOG_TYPE="A1-pthreads-1D-"
LOG_FILE_PROG="${LOG_PREFIX}${LOG_TYPE}${NOW}"
DAT_FILE_PROG="${LOG_FILE_PROG}"
if [[ ${#NX_ARRAY[*]} -ne ${#NB_ARRAY[*]} ]] ; then 
    error 5 "matrix size and block size"
else        
    EXECUTE_OPTIONS=""
    for (( i = 0 ; i < ${#NX_ARRAY[@]} ; i++ )) do
        LOG_FILE_PROG_VALUES="${LOG_FILE_PROG}-values-$i${TXT_SUFFIX}"
        DAT_FILE_PROG_TIMING="${DAT_FILE_PROG}-timing-$i${DAT_SUFFIX}"
        init_log_file $LOG_FILE_PROG_VALUES 
        init_log_file $DAT_FILE_PROG_TIMING
        add_comments_to_log_file "${LOG_FILE_PROG_VALUES}" "${PROG_NAME}"
        add_comments_to_log_file "${DAT_FILE_PROG_TIMING}" "${PROG_NAME}"
        EXECUTE_OPTIONS="${ALGORITHIM_OPTIONS} ${NX_ARRAY[$i]}"
        algorithm_execute "${PROG_NAME}" "${EXECUTE_OPTIONS}" "${LOG_FILE_PROG_VALUES}" "${DAT_FILE_PROG_TIMING}"
        LOG_FILE_PROG_VALUES="${LOG_FILE_PROG}"
        DAT_FILE_PROG_TIMING="${DAT_FILE_PROG}"
    done
fi

pause
exit 0
