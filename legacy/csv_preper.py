#!/Users/siu/miniconda3/bin/python

import pandas as pd
from pandas.api.types import is_numeric_dtype
import os
import numpy as np
from sklearn import metrics
import argparse
import math
from tabulate import tabulate

parser = argparse.ArgumentParser(description = 'Organize your csv file')
parser.add_argument('csv', help='your csv file') # 필요한 인수를 추가
args = parser.parse_args()

csv = args.csv

# ------------------------------------------------------------------------------------------------ #

def stat(df):
    stat = pd.DataFrame(columns = ["Column", "Type", "Null", "Uniq", "Value"])
    
    num = 0
    for col in df.columns:
        uniq = len(df[col].unique())
        
        # stat df에 값 채워넣기
        stat.loc[num]=[col, df[col].dtype, df[col].isna().sum(), uniq, df[col][df.index.min()]] 
        num +=1
    
    for i in stat.index:
        # str(object type)이 아닌 컬럼에서 Uniq 체크되는 것 삭제
        if stat.loc[i,'Type'] != 'object':
            stat.loc[i, 'Uniq'] = ''
        # uniq에서 전부 같을 때는 빼기
        if stat.loc[i, 'Uniq'] ==  df.shape[0]:
            stat.loc[i, 'Uniq'] = ''
        # null이 0인 것 보여주지 않기
        if stat.loc[i, 'Null'] == 0:
            stat.loc[i, 'Null'] = ''

    print(tabulate(stat, stat.columns, tablefmt="simple_outline", numalign="right"))
    return stat

def string_lister(string):
    slist = string.split(' ')
    slist = [x.strip() for x in slist]
    if '' in slist:
        slist.remove('') # 'a,' 이렇게 입력했을 때 오류 수정
    return slist

def int_checker(value):
    if type(value) == list:
        try:
            value = [int(x) for x in value]
        except:
            print("# [WARNING] 정수(인덱스)만 입력해주세요.")
            return
    else:
        try:
            value = int(value)
        except:
            print("# [WARNING] 정수(인덱스)만 입력해주세요.")
            return     
    return value

def find_dups(df, your_column):
    dpf = df.loc[df[your_column].duplicated(keep=False) == True]
    duplicated_items = list(set(dpf[your_column]))

    for i in duplicated_items:
        tf = dpf.loc[df[your_column] == i]
        tf = tf.copy()
        for col in tf.columns:
            if (len(set(tf[col])) == 1) & (col != your_column):
                tf.loc[tf[your_column]==i, col] = None
                tf = tf.dropna(axis=1)
        print(tf)
        print('')

# ------------------------------------------------------------------------------------------------ #

df = pd.read_csv(csv)  
fn = os.path.splitext(csv)[0]

check_point = {}  # 체크포인트 저장을 위한 딕셔너리
check_point['raw'] = df  # reload에 사용됨

# ------------------------------------------------------------------------------------------------ #
print('')
print("  This csv table has the number of rows and columns:", df.shape[0], df.shape[1])
stat(df)

mode = None
while mode != 'exit':
    string = input('### [menu] cols, dup, null, math, stat (exit) : ' )
    temp = string.split(',') # ,는 분리하기
    temp = [x.strip() for x in temp] # 공백제거
    mode = temp[0]

    if mode == '~':
        print('')
        print("This csv table has the number of rows and columns:", df.shape[0], df.shape[1])
        stat(df)
    
    if mode == 'df':
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 200)
        pd.set_option('display.max_colwidth', 20)
        pd.set_option('display.max_rows', 6)
        print('')
        print(df)

        pd.reset_option("display")

    if mode == 'reload':
        df = check_point['raw']
        print('')
        print("This csv table has the number of rows and columns:", df.shape[0], df.shape[1])
        stat(df)

    if mode == 'save':  # df를 원하는 이름의 csv로 export
        df.to_csv(f"{fn}_prep.csv", index=False)
        print(f"{fn}_prep.csv is saved.")

### [cols] ----------------------------------------------------------------------------------------- #

    if mode == 'cols':
        act = None
        while act != '~':
            string = input('### [menu > cols] keep, del, rename, reorder (~: menu): ' )
            temp = string.split(',') # ,는 분리하기
            temp = [x.strip() for x in temp] # 공백제거
            act = temp[0]

            if act == 'df':
                pd.set_option('display.max_columns', None)
                pd.set_option('display.width', 200)
                pd.set_option('display.max_colwidth', 20)
                pd.set_option('display.max_rows', 6)
                print('')
                print(df)

                pd.reset_option("display")
                
            if act == 'keep':
                string = input('# 컬럼을 입력하세요. 나머지는 컬럼은 삭제됩니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                df = df[names]
                

            if act == 'del':
                string = input('# 컬럼을 입력하세요. 삭제됩니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                df = df.drop(names, axis=1)
                stat(df)

            if act == 'rename':
                string = input('# 컬럼을 입력하세요. 이름을 변경할 수 있습니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                string = input(f'# 변경할 이름을 입력하세요 : ')
                if string == 'c':
                    continue

                df = df.rename(columns = {name : string})
                stat(df)

            if act == 'reorder':
                string = input('# 컬럼을 입력하세요. 복수 가능. 맨 앞으로 가져옵니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                
                cols_diff = [x for x in list(df.columns) if x not in names]
                cols_all_new = names + cols_diff

                df = df[cols_all_new]
                stat(df)

### [math] ----------------------------------------------------------------------------------------- #

    if mode == 'math':
        act = None
        while act != '~':
            string = input('### [menu > math] round, -log, reverse, asc, dsc (~: menu): ' )
            temp = string.split(',') # ,는 분리하기
            temp = [x.strip() for x in temp] # 공백제거
            act = temp[0]

            if act == 'df':
                pd.set_option('display.max_columns', None)
                pd.set_option('display.width', 200)
                pd.set_option('display.max_colwidth', 20)
                pd.set_option('display.max_rows', 6)
                print('')
                print(df)

                pd.reset_option("display")

            if act == 'round':
                string = input('# 컬럼을 입력하세요. 소수점 아래 2자리 반올림됩니다. (여러개 입력 가능) (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]

                for i in names:
                    df = df.round({i : 2})

                print(df[names])

            if act == '-log':
                string = input('# 컬럼을 입력하세요. -log(x)-9 계산됩니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue

                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                df.insert(idxs[0]+1, 'log_'+name, np.log(df[name])+9)

                print(df[[name, 'log_'+name]])

            if act == 'reverse':
                string = input('# 컬럼을 입력하세요. -(x) 계산됩니다. (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue

                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                df.insert(idxs[0]+1, 'rev_'+name, -(df[name]))

                print(df[[name, 'rev_'+name]])

            if act == 'asc':
                string = input('# 컬럼을 입력하세요. 오름차순 (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue

                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                df = df.sort_values(name, ascending=True)

                print(df[name])

            if act == 'dsc':
                string = input('# 컬럼을 입력하세요. 내림차순 (c: 취소): ')
                if string == 'c':
                    continue

                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드
                if len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))):  # 교집합 길이로 확인
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue

                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                df = df.sort_values(name, ascending=False)

                print(df[name])

### [dup] ----------------------------------------------------------------------------------------- #

    if mode == 'dup':
        #dups = df.loc[df.duplicated(keep=False) == True]
        #print(dups)
        #print(f"Above dulicated {len(dups)} rows are deleted.")
        #df = df.drop_duplicates(keep='first')
        
        act = None
        while act != '~':
            string = input('### [menu > dup] see, del, avg, subdiv (~: menu): ' )
            temp = string.split(',') # ,는 분리하기
            temp = [x.strip() for x in temp] # 공백제거
            act = temp[0]

            if act == 'df':
                pd.set_option('display.max_columns', None)
                pd.set_option('display.width', 200)
                pd.set_option('display.max_colwidth', 20)
                pd.set_option('display.max_rows', 6)
                print('')
                print(df)

                pd.reset_option("display")

            if act == 'see':
                string = input('# 컬럼을 입력하세요. 중복 확인 기준 (-1: 전체) (c: 취소): ')
                if string == 'c':
                    continue
                
                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드 
                if ( len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))) ) & (idxs != [-1]) : # 교집합 길이로 확인 # -1이 아닌 경우에만 작동하도록 설정
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue
                
                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                if idxs != [-1]:
                    find_dups(df, name)
                elif idxs == [-1]:
                    tf = df.loc[df.duplicated(keep=False) == True] # sort되지 않은 것
                    print(tf)

            if act == 'del':
                string = input('# 컬럼을 입력하세요. 중복 확인 기준 (-1: 전체) (c: 취소): ')
                if string == 'c':
                    continue
                
                #### 입력값을 정수인지 체크하는 코드
                idxs = string_lister(string)
                if int_checker(idxs) is None:
                    print("정수를 입력하세요.")
                    continue
                idxs = int_checker(idxs)

                # 인덱스에 포함되었는지 확인하는 코드 
                if ( len(idxs) != len(list(set(idxs) & set(list(range(0,len(df.columns)))))) ) & (idxs != [-1]) : # 교집합 길이로 확인 # -1이 아닌 경우에만 작동하도록 설정
                    print('컬럼 인덱스 값에 포함되지 않습니다.')
                    continue

                # 1개만 입력 확인하는 코드 
                if len(idxs) > 1:
                    print('컬럼을 하나만 정해주세요.')
                    continue
                #### --------------------------

                names = [df.columns[n] for n in idxs]
                name = names[0]

                if idxs != [-1]:
                    find_dups(df, name)
                    df = df.drop_duplicates(name, keep='first')
                    print(f"{name} 컬럼을 기준으로 중복된 중 첫번째 열만 남고 제거되었습니다.")
                elif idxs == [-1]:
                    df = df.drop_duplicates(keep='first')
                    print(f"전체 컬럼을 기준으로 중복된 항목이 제거되었습니다.")

            if act == 'avg':
                string = input('# 컬럼을 입력하세요. 중복 확인 기준 (-1: 전체) (c: 취소): ')
                if string == 'c':
                    continue
                
            

