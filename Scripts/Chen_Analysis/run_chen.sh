#!/bin/bash

echo "RUNNING CHEN ET AL. DATA ANALYSIS"

scripts_dir="Scripts/Chen_Analysis"
chen_data_dir="Data/Chen_Analysis"

if [ -d "$chen_data_dir" ]; then
    echo "Directory '$chen_data_dir' already exists."
    read -p "Replace? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Removing existing directory..."
            if rm -rf "$chen_data_dir"; then
                echo "Directory removed succesfully."
            else
                echo "Error in removing the directory."
                exit 1
            fi
            ;;
        *)
    esac
fi

if mkdir -p "$chen_data_dir"; then
    echo "Directory '$chen_data_dir' created sucessfully."
else
    echo "Error creating directory '$chen_data_dir'."
    exit 1
fi

csv_file="$chen_data_dir/chen_data.csv"
json_file="$chen_data_dir/chen_data.json"
if [ ! -f "$csv_file" ]; then
    wget "https://db.cngb.org/maya/api/get_export_data?dataset_id=MDB0000002&table_name=sample" -O "$json_file"
    python3 "$scripts_dir/convert_json_csv.py"
fi

echo "RUNNING DATA FILTERING"

python3 "$scripts_dir/filtering.py"

echo "DOWNLOADING SELECTED FASTAS"

bins_dir="Data/Chen_Analysis/Bins"
if [ ! -d "$bins_dir" ]; then
    bash "$scripts_dir/links.sh"
    echo "Done!"
else
    echo "Done!"
fi

echo "GENERATING SUMMARY STATISTICS"

bash "$scripts_dir/summary_gen.sh"

echo "CREATING SUMMARY REPORTS NOTEBOOK"

nb_dir="Data/Chen_Analysis/Reports"
mkdir -p "$nb_dir"

jupyter nbconvert --to notebook --execute \
    "$chen_data_dir/Reports/summary.ipynb" \
    --output-dir "$(realpath "$nb_dir")" \
    --output "summary"

echo "REMOVING THE OUTLIERS"

python3 "$scripts_dir/iqr_apply.py"

echo "CHECKM2 QUALITY CONTROL ANALYSIS"

if [ -d "$chen_data_dir/CheckM2" ]; then
    echo "Directory '$chen_data_dir/CheckM2' already exists."
    read -p "Overwrite? (y/N): " RESP

    case $RESP in
        [yY])
            echo "Overwriting CheckM2 analysis..."
            rm -rf "$chen_data_dir/CheckM2"
            bash "$scripts_dir/checkm2_analysis.sh"
            ;;
        *)
            echo "Skipping CheckM2 analysis."
            ;;
    esac
else
    echo "Running CheckM2 analysis..."
    bash "$scripts_dir/checkm2_analysis.sh"
fi

jupyter nbconvert --to notebook --execute \
    "$chen_data_dir/Reports/qc_report.ipynb" \
    --output-dir "$(realpath "$nb_dir")" \
    --output "qc_report"

echo "REMOVING LOW QUALITY BINS"

python3 "$scripts_dir/qc_apply.py"

