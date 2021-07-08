#!/usr/bin/env bash

# Check that only one argument was received
if [ $# -ne 5 ]; then
    echo "Five arguments (DIR, USER, PASS, HOST, PORT) expected, $# received"
    exit 1
fi

# Parsing arguments
USER=$2
PASSWORD=$3
HOST=$4
PORT=$5

for f in $(ls $1/*.db)
do
echo "Processing: $f"
filename=$(basename -- "$f")
DB="${filename%.*}"
createdb -h $HOST -p $PORT -O $USER -w -U $USER $DB
echo "pgloader $f postgresql://$USER:$PASSWORD@$HOST:$PORT/$DB"
sqlite3 $f .dump >> $f.dump
#psql -d $DB -U $USER -w < $f.dump
#pgloader $f postgresql://$USER:$PASSWORD@$HOST:$PORT/$DB
done