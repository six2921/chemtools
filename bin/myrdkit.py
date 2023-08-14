import os
import pandas as pd
import numpy as np
import glob
import itertools
import sys
import time
import math
import rdkit
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit import Chem
from rdkit.Chem import rdFMCS as MCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolsToGridImage
from pathos.multiprocessing import ProcessingPool
num_cores = os.cpu_count()
pool = ProcessingPool(num_cores)

def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def molidx(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def molsgrid(smiles_list, legend=None, max_num=10): 
    mols = [Chem.MolFromSmiles(x) for x in smiles_list][:max_num]
    if legend is not None: legend = [str(i) for i in legend][:max_num]
        
    image = Draw.MolsToGridImage(mols, molsPerRow=5, legends=legend)
    
    return image

def split_list(mylist):
    """You can write docstring here"""
    cores = os.cpu_count()
    list_size = len(mylist)
    n = int(list_size/cores) if list_size > cores*5 else 4
    spl = [mylist[i * n:(i + 1) * n] for i in range((len(mylist) + n - 1) // n )]
    return spl

def sim_search(query_smiles_list, db_smiles_list, cutoff=0.6, query_id_list=None):
    start = time.time()  # 시작 시간 저장

    # 시리즈 데이터가 아닌 경우 시리즈 데이터로 변환
    if isinstance(query_smiles_list, pd.Series) is False: query_smiles_list = pd.Series(query_smiles_list) 
    if isinstance(db_smiles_list, pd.Series) is False: db_smiles_list = pd.Series(db_smiles_list) 

    # query_id가 없는 경우 인덱스 생성
    if query_id_list is None: query_id_list = list(range(0, len(query_smiles_list))) 

    # 병렬 처리할 함수 정의
    def generate_fps(smiles_list):
        fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in smiles_list]
        return fps
    
    # (pathos parallel) generate fps for query
    print('# Generting fingerprints...')    # 계산에 대한 문구 프린트
    
    mylist = split_list(query_smiles_list)    # 사용할 데이터 정의 및 분할 
    myfunc = generate_fps    # 사용할 함수 정의
    calculated = list(pool.map(myfunc, mylist))
    calculated = list(itertools.chain.from_iterable(calculated))
    qfps = calculated    # 계산 결과를 사용하고자 하는 변수에 저장 
    
    # id 혹은 smiles를 잚못 입력한 경우 terminate 조건
    if len(query_smiles_list) != len(qfps): return print('##### [TERMINATED] Num Smiles != Num Genereated Fingerprints #####')
    if len(query_id_list) != len(qfps): return print('##### [TERMINATED] Num IDs != Num Genereated Fingerprints #####')

    # (pathos parallel) generate fps for db
    # print('# Generting fingerprints for db smiles...')    # 계산에 대한 문구 프린트
    mylist = split_list(db_smiles_list)    # 사용할 데이터 정의 및 분할 
    myfunc = generate_fps    # 사용할 함수 정의
    calculated = list(pool.map(generate_fps, mylist))
    calculated = list(itertools.chain.from_iterable(calculated))
    dbfps = calculated    # 계산 결과를 사용하고자 하는 변수에 저장 

    # genreate empty dataframe
    df = pd.DataFrame(columns = ['id', 'smiles', 'matched_freq', 'matched_smiles'])
    df['id'] = query_id_list
    df['smiles'] = query_smiles_list
    
    # Searching
    print('# Searching...: ', end='')
    matched_freq = []
    matched_smiles = []
    for fp, idx in zip(qfps, range(0, len(qfps))):
        print(idx,' ', end='')
        
        # rdkit Similarity
        sim_value = [DataStructs.FingerprintSimilarity(fp, fp2) for fp2 in dbfps]
        
        # sim_value를 구하고 sort, size 체크
        sim_value = pd.Series(sim_value)
        matched = sim_value.loc[sim_value > cutoff]
        matched = matched.sort_values(ascending=False)
        matched_freq.append(len(matched)) # append
        
        # index를 추출하여 db_smiles_series에서 value 가져오기
        smiles = list(db_smiles_list.loc[matched[:5].index].values)
        matched_smiles.append(smiles) # append
        
    # 값 기록하기
    df['matched_freq'] = matched_freq
    df['matched_smiles'] = matched_smiles

    end = time.time()  # 끝나는 시간 저장
    calc_time = round(end-start, 2) # 계산 시간
    print("\n# The time spent:", calc_time)
    
    return df

def ss_search(query_smiles_list, db_smiles_list, query_id_list=None):
    start = time.time()  # 시작 시간 저장

    # 시리즈 데이터가 아닌 경우 시리즈 데이터로 변환
    if isinstance(query_smiles_list, pd.Series) is False: query_smiles_list = pd.Series(query_smiles_list) 
    if isinstance(db_smiles_list, pd.Series) is False: db_smiles_list = pd.Series(db_smiles_list) 

    # query_id가 없는 경우 인덱스 생성
    if query_id_list is None: query_id_list = list(range(0, len(query_smiles_list))) 
    
    # 병렬 처리할 함수 정의
    def generate_mols(smiles_list):
        mols = [Chem.MolFromSmiles(x) for x in smiles_list]
        return mols
    
    # (pathos parallel) generate fps for query
    print('# Generting mols...')    # 계산에 대한 문구 프린트
    mylist = split_list(query_smiles_list)    # 사용할 데이터 정의 및 분할 
    myfunc = generate_mols    # 사용할 함수 정의
    calculated = list(pool.map(myfunc, mylist))
    calculated = list(itertools.chain.from_iterable(calculated))
    qmols = calculated    # 계산 결과를 사용하고자 하는 변수에 저장 
    
    # id 혹은 smiles를 잚못 입력한 경우 terminate 조건
    if len(query_smiles_list) != len(qmols): return print('##### [TERMINATED] Num Smiles != Num Genereated Molecules #####')
    if len(query_id_list) != len(qmols): return print('##### [TERMINATED] Num IDs != Num Genereated Molecules #####')

    # (pathos parallel) generate fps for db
    # print('# Generting mols for db smiles...')    # 계산에 대한 문구 프린트
    mylist = split_list(db_smiles_list)    # 사용할 데이터 정의 및 분할 
    myfunc = generate_mols    # 사용할 함수 정의
    calculated = list(pool.map(myfunc, mylist))
    calculated = list(itertools.chain.from_iterable(calculated))
    dbmols = calculated    # 계산 결과를 사용하고자 하는 변수에 저장 

    # genreate empty dataframe
    df = pd.DataFrame(columns = ['id', 'smiles', 'matched_freq', 'matched_smiles'])
    df['id'] = query_id_list
    df['smiles'] = query_smiles_list
    
    # Searching 
    print('# Searching...: ', end='')
    matched_freq = []
    matched_smiles = []
    for mol, idx in zip(qmols, range(0, len(qmols))):
        print(idx,' ', end='')
        
        # rdkit substructure
        ss_bool = [x.HasSubstructMatch(mol) for x in dbmols]
        
        # matched를 boolean으로 size 체크
        ss_bool = pd.Series(ss_bool)
        matched = ss_bool.loc[ss_bool==True]
        matched_freq.append(len(matched))

        # index를 추출하여 db_smiles_series에서 value 가져오기
        index = ss_bool.loc[ss_bool==True].index
        smiles = list(db_smiles_list.loc[index].values)
        matched_smiles.append(smiles)

    # 값 기록하기
    df['matched_freq'] = matched_freq
    df['matched_smiles'] = matched_smiles

    end = time.time()  # 끝나는 시간 저장
    calc_time = round(end-start, 2) # 계산 시간
    print("\n# The time spent:", calc_time)
    
    return df

def butina_cluster(smiles_list, coefficient_cut=0.6, id_list=None):
    start = time.time()  # 시작 시간 저장

    # 시리즈 데이터가 아닌 경우 시리즈 데이터로 변환
    if isinstance(smiles_list, pd.Series) is False: smiles_list = pd.Series(smiles_list) 

    # query_id가 없는 경우 인덱스 생성
    if id_list is None: id_list = list(range(0, len(smiles_list))) 
        
    # 병렬 처리할 함수 정의
    def generate_fps(smiles_list):
        fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in smiles_list]
        return fps

    # (pathos parallel) generate fps for query
    print('# Generting fingerprints...')    # 계산에 대한 문구 프린트
    mylist = split_list(smiles_list)    # 사용할 데이터 정의 및 분할 
    myfunc = generate_fps    # 사용할 함수 정의
    calculated = list(pool.map(myfunc, mylist))
    calculated = list(itertools.chain.from_iterable(calculated))
    fps = calculated    # 계산 결과를 사용하고자 하는 변수에 저장 

    # id 혹은 smiles를 잚못 입력한 경우 terminate 조건
    if len(smiles_list) != len(fps): return print('##### [TERMINATED] Num Smiles != Num Generated Fingerprints #####')
    if len(id_list) != len(fps): return print('##### [TERMINATED] Num IDs != Num Genereated Molecules #####')

    # generate the distance matrix:
    print('# Generating distance matrix & clustering')
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    
    # now cluster the data:
    cutoff = 1-coefficient_cut
    nfps = len(fps)
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)

    # tuple data to list and sort list by each length of component
    clusters = list(cs)
    clusters.sort(key = lambda x: len(x), reverse=True)    # how to sort tuple data
    clusters_size_list = [len(x) for x in clusters if len(x) != 1]
    print('# Size of Clusters: ', clusters_size_list)

    # genreate empty dataframe
    df = pd.DataFrame(columns = ['id', 'smiles', 'cl_butina', 'cl_size'])
    df['id'] = id_list
    df['smiles'] = smiles_list

    clusters_list = range(0, len(clusters))
    for i in clusters_list:
        size = len(clusters[i])
        if size == 1:
            for idx in clusters[i]:
                df.loc[idx, 'cl_butina'] = 0  # 클러스터 크기 1인 것은 0
                df.loc[idx, 'cl_size'] = size
        else:
            for idx in clusters[i]:
                df.loc[idx, 'cl_butina'] = i+1  # 0이 없도록
                df.loc[idx, 'cl_size'] = size 

    end = time.time()  # 끝나는 시간 저장
    calc_time = round(end-start, 2) # 계산 시간
    print("# The spent time:", calc_time)

    return df 
