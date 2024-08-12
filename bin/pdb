#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: $0 <PDB ID>"
  exit 1
fi

pdb_id=$1

curl -s https://files.rcsb.org/view/${pdb_id}.pdb > ${pdb_id}.pdb

curl -s https://files.rcsb.org/view/${pdb_id}.pdb | grep -E '^HET\s' | awk '{print $2}' | uniq | while read lig; do
    result=$(curl -s https://files.rcsb.org/ligands/view/${lig}_model.sdf | grep -A 1 'ISO_SMILES' | sed 's/^> <//;s/>//' | tail -n 1)
    echo "$lig: $result"
done

echo "PDB 파일이 다운로드되었습니다."
