import argparse
import pandas as pd

def read_and_display_columns(file_path):
    df = pd.read_csv(file_path)
    print("CSV 파일의 컬럼들:")
    for idx, col in enumerate(df.columns):
        print(f"{idx}: {col}")
    return df

def select_columns_and_save(df, columns, output_file):
    selected_columns = [df.columns[int(col)] for col in columns]
    df_selected = df[selected_columns]
    df_selected.to_csv(output_file, index=False)
    print(f"선택한 컬럼들로 새로운 CSV 파일 '{output_file}'이(가) 생성되었습니다.")

def main():
    parser = argparse.ArgumentParser(description="CSV 파일에서 컬럼들을 선택하여 새로운 CSV 파일을 생성합니다.")
    parser.add_argument('input_file', type=str, help="입력 CSV 파일 경로")
    args = parser.parse_args()

    df = read_and_display_columns(args.input_file)
    
    columns = input("저장할 컬럼들을 띄어쓰기로 입력 (예: 1 3 10 14): ").split()
    output_file_name = args.input_file.replace('.csv', '') + '-subset.csv'
    select_columns_and_save(df, columns, output_file_name)

if __name__ == "__main__":
    main()
