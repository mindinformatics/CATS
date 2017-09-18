# remove all of the GSE and GPL files from the directory
rm GSE*
rm GPL*

# set Datasets.csv to the column names
head -1 Datasets.csv > DatasetsCopy.csv
rm Datasets.csv
mv DatasetsCopy.csv Datasets.csv

