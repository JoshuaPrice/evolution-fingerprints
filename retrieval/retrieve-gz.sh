#!/bin/bash
# Bash script that downloads data from the link provided, stores it in the location designated

echo -e "Please enter the file address you would like to download."
read file_url
echo -e "Please designate the directory in which to store this file, e.g. /shared/data"
read directory_loc
echo -e "What would you like to call this file? Should end in .gz if the original file does."
read file_name
echo "Downloading $file_url..."

wget -O "$directory_loc/$file_name" "$file_url"

echo -e "Finished downloading file. Now unzipping file."

gunzip "$directory_loc/$file_name"

echo -e "Done unzipping file. Enjoy your newfound data!"



