import requests
import json
from io import StringIO
import pandas as pd
import numpy as np
import os
import requests
import pandas as pd
import xml.etree.ElementTree as ET
import math

ignore_list = open('ignore_ligands.txt', 'r').read().split('\n')

def is_json_key_present(json, key):
    try:
        buf = json[key]
    except KeyError:
        return False

    return True

def pdb2ligands(pdb):
    res = requests.get(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')

    if 'rcsb_entry_container_identifiers' not in res.json():
        pass  # pdb가 아직 릴리즈 전인 경우
    elif 'non_polymer_entity_ids' not in res.json()['rcsb_entry_container_identifiers']:
        pass  # ligand가 없는 경우
    else:
        comp_list = dict()
        het_entities = res.json()['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
        for het in het_entities:
            res = requests.get(f'https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb}/{het}')
            comp_id = res.json()['pdbx_entity_nonpoly']['comp_id']
            if comp_id not in ignore_list:
                res = requests.get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}')
                smiles = res.json()['rcsb_chem_comp_descriptor']['smiles']
                if len(smiles) > 7: # smiles 글자 수가 7개 미만이면 빼는 걸로
                    comp_list[comp_id] = smiles
                else:
                    pass
            ## 여기에 binding affinity 찾아서 집어넣기
            else:
                pass
        
        return comp_list

def uniprot2pdbs(uniprot_id):
    res = requests.get(f'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}')
    df = pd.read_json(StringIO(json.dumps(res.json()[f"{uniprot_id}"])))
    df = df.groupby('pdb_id', as_index=False).first()
    df = df[['pdb_id', 'start', 'end', 'coverage', 'resolution', 'experimental_method']]

    return df

def uniprot2pdbligands(uniprot_id, count_cutoff=100):
    from check_pdb import uniprot2pdbs
    df = uniprot2pdbs(uniprot_id)
    df = df.sort_values('coverage', ascending=False)    # coverage 순서로 정렬
    df['numligs'] = None    # pdb2ligands 검색 결과가 없는 경우 생기는 오류 방지
    print(uniprot_id, 'Total PDBs:', df.shape[0])
    print('----------------------')

    from check_pdb import pdb2ligands
    lig_count = 0
    for idx in df.index:
        pdb = df.pdb_id[idx]
        coverage = df.coverage[idx]

        if pdb2ligands(pdb) is not None:
            lig_dict = dict()    # 이전 loop에서 온 변수를 그대로 가져가지 않도록 초기화
            lig_dict = pdb2ligands(pdb)
            count = len(lig_dict.values())
            print(pdb, coverage, count, lig_dict)
            
            # 데이터 프레임에 집어넣기 
            df.loc[df['pdb_id'] == pdb, 'numligs'] = int(count)    
            for num, key in zip(range(0, count), lig_dict.keys()):
                df.loc[df['pdb_id'] == pdb, 'het'+str(num)] = lig_dict[key]
        else:
            lig_dict = dict()
        
        usable = df.loc[(df['coverage'] > 0.5) & (df['numligs'] > 0)]
        if usable.shape[0] >= count_cutoff:
            break
        else:
            pass
    
    return df

def uniprot2domaininfo(uniprot_id):
    # 불러오기
    URL = f'https://www.uniprot.org/uniprot/{uniprot_id}.xml' 
    response = requests.get(URL) 
    status = response.status_code 
    text = response.text
    root = ET.fromstring(response.text)

    # Protein 이름 프린트
    protein_name = root.findall("./{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}protein/{http://uniprot.org/uniprot}recommendedName/{http://uniprot.org/uniprot}fullName")[0].text
    print("Target name: ", protein_name)

    # PDB 이름들을 읽어서 리스트로 만들고 빈 데이터 프레임에 집어넣기 
    df = pd.DataFrame()
    pdb_list = []
    for n in root.iter('{http://uniprot.org/uniprot}dbReference'):
        if n.attrib['type'] == 'PDB':
            pdb_id = n.attrib['id']
            pdb_list.append(pdb_id)
    df['PDB'] = pdb_list

    # 첫번째, Uniprot domain 정보를 이름:서열 딕셔너리로 저장하기 
    domains_dict = dict()
    for n in root.iter('{http://uniprot.org/uniprot}feature'):
        if n.attrib['type'] == 'domain':    # type == domain인 값을 불러옴
            # begin 숫자, end 숫자로부터 리스트를 만들어서 딕셔너리 값으로 저장
            name = n.attrib['description']
            begin = int(n[0][0].attrib['position'])
            end = int(n[0][1].attrib['position'])
            seq = list(range(begin, end))    
            domains_dict[name] = seq    # 딕셔너리에 저장
            domain_info = name+' '+str(begin)+'-'+str(end)   # CXC 10-30 과 같은 포맷으로 데이터 프레임에 기록
            df[name] = domain_info
            print(domain_info)

    # 두번째, uniprot에 등록된 PDB 정보를 읽어서 그 안에 chain 정보를 가져온다. 
    for n in root.iter('{http://uniprot.org/uniprot}dbReference'):
        if n.attrib['type'] == 'PDB':
            for i in n.findall('{http://uniprot.org/uniprot}property'):
                if i.attrib['type'] == 'chains':
                    pdb_id = n.attrib['id']
                    chain_info = i.attrib['value']
                    df.loc[df['PDB'] == pdb_id, 'PDB_CHAIN'] = chain_info  
                    #print(pdb_id, chain_info)
    
    df['PDB_CHAIN'] = df['PDB_CHAIN'].fillna('NoInfo=0-0')   # chain 정보가 없는 PDB가 있어서 임의의 값을 넣어줌 

    # 세번째, 읽어온 chain 정보를 split하여 begin-end로부터 리스트 형태의 서열을 만든다. 모든 chain 서열을 합친다. 
    for idx in df.index:
        value = df.loc[idx, 'PDB_CHAIN']    # value 값을 뽑아 서열 리스트를 만든다. 
        seq_sum = []    # 리스트가 여러 개일 경우 합치기 위해 빈 리스트를 만든다. 
        seq = value.split(',')    # {A=10-20, B=30-40} 과 같은 형태를 1차 분리
        for s in seq:
            seq = s.split('=')[1]    # 2차 분리
            begin = int(seq.split('-')[0])
            end = int(seq.split('-')[1])
            seq_chain = list(range(begin, end+1))    # range 함수로 나열하여 리스트로 만듬
            seq_sum = seq_sum+seq_chain
        
        # 네번째, uniprot domain 딕셔너리에 있는 모든 서열 리스트와 PDB 서열 리스트(세번째)를 비교
        for name, domain_seq in zip(domains_dict.keys(), domains_dict.values()):
            failed_seq = []
            for seq in domain_seq:
                if seq not in seq_sum:
                    failed_seq.append(seq)
            loss_ratio = len(failed_seq)/len(domain_seq)
            if loss_ratio < 0.9:
                df.loc[idx, name] = 'Y'
                # print(idx, pdb_id, name, 'Success')
            else:
                pass
                # print(idx, pdb_id, name, 'Fail')   

    # 중복제거
    df = df.drop_duplicates()
    
    return df


