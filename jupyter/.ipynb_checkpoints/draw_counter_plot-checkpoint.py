# plot_histogram.py

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ipywidgets import interactive, widgets, HBox, VBox, Layout
from IPython.display import display, clear_output

# 히스토그램 그리는 함수 정의
def plot_histogram(df, column_name, class_column, min_val, max_val, bin_width, aspect_ratio, save_path=None):
    # 데이터 정렬
    if class_column != 'None':
        data = df[(df[column_name] >= min_val) & (df[column_name] <= max_val)].sort_values(by=class_column)
    else:
        data = df[(df[column_name] >= min_val) & (df[column_name] <= max_val)]
    
    plt.figure(figsize=(6 * aspect_ratio, 6))
    ax = plt.gca()
    
    # width 매개변수를 사용하여 막대의 너비를 조정하고 간격을 둠
    bar_width = bin_width * 0.8  # 막대의 너비를 줄여서 간격 생성
    
    if class_column != 'None':
        sns.histplot(data, x=column_name, hue=class_column, bins=range(min_val, max_val + bin_width, bin_width), 
                     multiple='dodge', kde=False, ax=ax, palette='tab10', shrink=0.8)
    else:
        sns.histplot(data[column_name], bins=range(min_val, max_val + bin_width, bin_width), kde=False, color='blue', ax=ax, shrink=0.8)
    
    # 막대 위에 빈도 수 표시
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.text(p.get_x() + p.get_width() / 2., height, int(height), ha='center', va='bottom')
    
    # x축 눈금 간격 설정
    plt.xticks(range(min_val, max_val + bin_width, bin_width))
    
    # Y축 그리드 설정
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5)
    
    plt.xlabel(column_name)
    plt.ylabel('Count')
    
    plt.show()

def create_widgets_and_display(df):
    # int 또는 float 타입의 컬럼만 선택 가능하도록 필터링
    numeric_columns = df.select_dtypes(include=['int', 'float']).columns.tolist()

    # 유니크 값이 10개 이하인 컬럼 필터링
    class_columns = ['None'] + [col for col in df.columns if df[col].nunique() <= 10]

    # 위젯 정의
    value_widget = widgets.Dropdown(
        options=numeric_columns,
        value=numeric_columns[0],
        description='Value:',
        style={'description_width': 'initial'},
        layout=Layout(width='70%')
    )

    type_widget = widgets.Dropdown(
        options=class_columns,
        value='None',
        description='Type:',
        style={'description_width': 'initial'},
        layout=Layout(width='70%')
    )

    min_val_widget = widgets.IntText(description='Min Value:', style={'description_width': 'initial'}, layout=Layout(width='70%'))
    max_val_widget = widgets.IntText(description='Max Value:', style={'description_width': 'initial'}, layout=Layout(width='70%'))
    bin_width_widget = widgets.IntText(value=1, description='Bin Width:', style={'description_width': 'initial'}, layout=Layout(width='70%'))
    image_ratio_widget = widgets.FloatText(value=1.0, description='Image Ratio:', style={'description_width': 'initial'}, layout=Layout(width='70%'))

    # 그래프를 저장하는 함수
    def save_graph(*args):
        save_path = f"{value_widget.value}.png"
        plot_histogram(df, value_widget.value, type_widget.value, min_val_widget.value, max_val_widget.value, bin_width_widget.value, image_ratio_widget.value, save_path)

    # 컬럼 선택 시 최소값과 최대값을 업데이트하는 함수
    def update_min_max(*args):
        column_name = value_widget.value
        min_val_widget.value = df[column_name].min()
        max_val_widget.value = df[column_name].max()

    # 컬럼 이름 선택 시 트리거
    value_widget.observe(update_min_max, names='value')

    # 초기 값 설정
    update_min_max()

    # 입력란과 그래프를 각각 세로로 가운데 정렬
    input_widgets = VBox(
        [value_widget, type_widget, min_val_widget, max_val_widget, bin_width_widget, image_ratio_widget],
        layout=Layout(align_items='center')
    )

    output = widgets.Output()

    def update_graph(*args):
        with output:
            clear_output(wait=True)
            plot_histogram(df, value_widget.value, type_widget.value, min_val_widget.value, max_val_widget.value, bin_width_widget.value, image_ratio_widget.value)

    # 각 위젯의 값이 변경될 때마다 그래프 업데이트
    for widget in [value_widget, type_widget, min_val_widget, max_val_widget, bin_width_widget, image_ratio_widget]:
        widget.observe(update_graph, 'value')

    # 저장 버튼 생성 및 클릭 이벤트 연결
    save_button = widgets.Button(description="Export PNG", layout=Layout(width='70%'))
    save_button.on_click(save_graph)

    # 초기 그래프 생성
    update_graph()

    # 입력란과 그래프를 나란히 배치
    ui = HBox([VBox([input_widgets, save_button]), output], layout=Layout(align_items='center'))
    display(ui)
