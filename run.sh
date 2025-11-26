#!/bin/bash

echo "Starting..."

data_dir="Data"

if [ -d "$data_dir" ]; then
    echo "Directory '$data_dir' already exists."
    read -p "Replace? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Removing existing directory..."
            if rm -rf "$data_dir"; then
                echo "Directory removed succesfully."
            else
                echo "Error in removing the directory."
                exit 1
            fi
            ;;
        *)
    esac
fi

if mkdir -p "$data_dir"; then
    echo "Directory '$data_dir' created sucessfully."
else
    echo "Error creating directory '$data_dir'."
    exit 1
fi

bash Scripts/Chen_Analysis/run_chen.sh