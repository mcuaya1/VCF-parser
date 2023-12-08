#!/usr/bin/bash

GENOME=S288C_reference_sequence_R64-4-1_20230823.fsa
awk '/>/{print}' S288C_reference_sequence_R64-4-1_20230823.fsa | awk '{gsub(/\[|\]|=/, " "); print}' | awk '{print $9,$10}' | awk '{gsub(/ /,"_"); print'} > tmp_new.txt
awk '/>/{print}' S288C_reference_sequence_R64-4-1_20230823.fsa | awk '{sub(">","");print}' | awk '{print $1}' > tmp_old.txt

old_chrom=tmp_old.txt
new_chrom=tmp_new.txt

backup_file=ORIGINAL_S288C_reference_sequence_R64-4-1_20230823.fsa

if ! [ -f $backup_file ]; then
	echo "Backing up original file..."
	cat $GENOME > $backup_file
fi

echo "Changing chromosome location header using the slowest possible method..."
while read -r -u 3 new_name && read -r -u 4 old_name; do
	sed -i "s/${old_name}/${new_name}/" "${GENOME}"
done  3<$new_chrom 4<$old_chrom
echo "Removing tmp files..."
rm tmp_old.txt
rm tmp_new.txt
echo "Finished.."
