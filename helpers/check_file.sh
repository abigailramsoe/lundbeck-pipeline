# Function to check if file exists
if [[ ! -f $1 ]]; then
    echo "Error: File $1 does not exist."
    exit 1
fi
