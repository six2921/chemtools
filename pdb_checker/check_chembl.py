
def assay(target_chembl_id, pchembl_value__isnull=False):
    import pandas as pd
    from chembl_webresource_client.new_client import new_client
    
    # CHEMBL API로 assay를 json으로 받아오기 
    activities = new_client.activity
    res = activities.filter(target_chembl_id=target_chembl_id, pchembl_value__isnull=pchembl_value__isnull)
    print(f"# Molecules in {target_chembl_id}: ", len(res))

    # json 형태를 합쳐야 하는데, value가 [] 형태이지 않으면 DataFrame.from_dict를 할 수 없다. 그래서 각 프로퍼티마다 값을 모아서 리스트로 만들고, 모든 리스트를 key:value 형태로 딕셔너리로 저장. 그리고 DataFrame.from_dict 실행

    # 프로퍼티를 리스트 형태로 만들어두고
    props = list(res[0].keys())

    # 각 프로퍼티마다 빈 리스트를 생성하고, 빈 리스트를 딕셔너리에 모음
    props_dict = {}
    for i in props:
        props_dict[i] = []

    # res[0], res[1]...이 물질들인데 여기서 프로퍼티마다 값들을 빈 리스트에 저장 
    print('# Loading(%): ', end='')
    n=0
    while res[n] != None:
        for i in props:
            value = res[n][i]
            props_dict[i].append(value)
        n+=1
        if n%50 == 0 : print(str(int((n/len(res))*100))+'%,', end=' ')

    # 각 프로퍼티에 대한 리스트가 저장된 딕셔너리(props_dict)를 데이터프레임으로 저장하기 
    df = pd.DataFrame.from_dict(props_dict)
    
    return df