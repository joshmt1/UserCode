#!/bin/sh

#to use: just set old and new below, and do ./switchCastorDir.sh
old='V00-01-01'
new='V00-01-02'

mkdir cfg_files_$old

thelist=`ls crab.*.cfg`

for filespec in $thelist; do

    echo $filespec
    #note --double quotes tells the shell to expand the variables....
    thedir=`awk -F= '/user_remote_dir/ {print $2}' $filespec | sed 's@/user/j/joshmt/@/castor/cern.ch/user/j/joshmt/@' | sed "s/$old/$new/"`
    echo $thedir
    nsmkdir -p $thedir
    sed "s/$old/$new/" $filespec > $filespec.new
    mv $filespec cfg_files_$old
    mv $filespec.new $filespec
done
