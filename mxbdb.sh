#!/bin/bash
run_layout=false
run_chen=false
run_all=false

# Defining arguments
if [[ "$#" -eq 1 && "$1" == "--help" ]]; then
    cat << 'EOF'
In development
EOF
    exit 0
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --layout)
            run_layout=true
            shift
            ;;
        --chen)
            run_chen=true
            shift
            ;;
        --all)
            run_all=true
            shift
            ;;
        *)
            echo "Invalid argument: $1" >&2
            exit 1
            ;;
    esac
done

if $run_all; then
    if $run_chen; then
        echo "Error: Argument --all was used alongside to specific subprojects, this argument enables all."
        exit 1
    fi
    run_chen=true
fi

if ! $run_chen; then
    echo "Error: The whole database or specific subprojects are needed to run MXD"
    exit 1
fi

echo "Starting"
echo "Creating Data Directories"

data_dir="Data"
tmp_dir="tmp"
log_dir="log"

sub_data_dirs=("$data_dir/Entities" "$data_dir/Raw" "$data_dir/Raw/Bins" "$data_dir/Raw/Processed" "$data_dir/Reports") 

if [ -d "$data_dir" ]; then
    echo "Directory '$data_dir' already exists."
    read -p "Replace? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Removing existing directory..."
            if rm -rf "$data_dir"; then
                echo "Data directories removed successfully."
            
            elif rm -rf "$tmp_dir"; then
                echo "Tmp directory removed sucessfully."

            elif rm -rf "$log_dir"; then
                echo "Log directory removed sucessfully."

            else
                echo "Error removing the directories."
                exit 1
            fi
            ;;
        *)
            echo "Keeping existing directory. Skipping creation."
            skip_data=true
            ;;
    esac
fi

if [ -z "$skip_data" ]; then
    for dir in "${sub_data_dirs[@]}"; do
        if mkdir -p "$dir"; then
            echo "Created: $dir"
        else
            echo "Failed to create: $dir"
            exit 1
        fi
    done

    if mkdir -p "$tmp_dir"; then
        echo "Tmp directory created"

        cat > "$tmp_dir"/metadata.csv <<EOF
        bin,sample,project,id_study,coord,date
EOF
    else
        echo "Failed to create: $tmp_dir"
    fi
fi

if mkdir -p "$log_dir"; then
    echo "Log directory created"

    cat > "$log_dir"/run.log <<EOF
##########     LOG FILE     ##########

This file records the duration of each checkpoint, memory and space usage, it is intended for developers,
feel free to delete it when the database is built on your machine.

Run: $(date)

Process: $0


EOF

else
    echo "Failed to create: $log_dir"
fi

if $run_chen; then
    python3 "source/chen_data/run_chen.py"
fi













