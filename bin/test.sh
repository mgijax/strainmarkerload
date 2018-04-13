#!/bin/sh

count=10
min_records=11
echo "count: $count min_records: $min_records"
if [ $count -lt $min_records ]
then
    echo '<'
else
    echo '>='
fi
