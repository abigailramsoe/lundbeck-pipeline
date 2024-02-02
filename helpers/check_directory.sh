# Function to check if directory exists and is writable
if [[ ! -d $1 ]]; then
    echo "Directory $1 does not exist. Attempting to create it..."
    mkdir -p "$1"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create directory $1."
        exit 1
    fi
elif [[ ! -w $1 ]]; then
    echo "Error: Directory $1 is not writable."
    exit 1
fi
