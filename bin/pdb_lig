#!/bin/bash

if [ $# -ne 2 ]; then
  echo "Usage: $0 <PDB ID> <ligand name>"
  exit 1
fi

pdb_id=$1
ligand_name=$2

# 리간드 SDF 파일 다운로드 (ideal 및 model)
curl -s https://files.rcsb.org/ligands/view/${ligand_name}_ideal.sdf > ${pdb_id}_${ligand_name}_ideal.sdf
curl -s https://files.rcsb.org/ligands/view/${ligand_name}_model.sdf > ${pdb_id}_${ligand_name}_model.sdf

echo "${pdb_id}_${ligand_name}_ideal.sdf 파일이 다운로드되었습니다."
echo "${pdb_id}_${ligand_name}_model.sdf 파일이 다운로드되었습니다."

# PDB 파일 다운로드
curl -s https://files.rcsb.org/view/${pdb_id}.pdb > ${pdb_id}.pdb

# pdb_delresname 명령어 실행 (pdb-tools 필요)
pdb_delresname -${ligand_name} ${pdb_id}.pdb > ${pdb_id}_no_${ligand_name}.pdb
echo "${pdb_id}_no_${ligand_name}.pdb 파일이 생성되었습니다."

