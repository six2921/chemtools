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
    
    fig, ax = plt.subplots(figsize=(6 * aspect_ratio, 6))
    
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
    
    # x축 레이블과 눈금 폰트 사이즈 설정
    ax.set_xlabel(column_name, fontsize=14, fontweight='heavy')  # x축 레이블 폰트 사이즈 설정
    
    # Y축 그리드 설정
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5)
    
    plt.xlabel(column_name)
    plt.ylabel('Count')

    if save_path:
        plt.savefig(save_path)
    
    plt.show()

def chart_num(df):
    # int 또는 float 타입의 컬럼만 선택 가능하도록 필터링
    numeric_columns = df.select_dtypes(include=['int', 'float']).columns.tolist()

    # 유니크 값이 10개 이하인 컬럼 필터링
    class_columns = ['None'] + [col for col in df.columns if df[col].nunique() <= 10]

    # 위젯 정의
    value_widget = widgets.Dropdown(
        options=numeric_columns,
        value=numeric_columns[0],
        description='Column:',
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
        layout=Layout(align_items='flex-start')
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

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ipywidgets import interactive, widgets, HBox, VBox, Layout
from IPython.display import display, clear_output

# 누적 바 차트를 그리는 함수 정의
def plot_stacked_bar_chart(df, column):
    data = df[column].astype(str).value_counts(normalize=True).reset_index()
    data.columns = [column, 'Proportion']
    data['Count'] = df[column].astype(str).value_counts().values

    fig, ax = plt.subplots(figsize=(3, 3))
    colors = plt.cm.tab10(range(len(data)))
    ax.bar(data[column], data['Proportion'], color=colors)
    
    # Bar labels
    for i in range(len(data)):
        ax.text(i, data['Proportion'][i] / 2, f"{data['Proportion'][i] * 100:.1f}%\n({data['Count'][i]})", 
                ha='center', va='center', color='black')
                
    ax.set_ylabel('Proportion')
    ax.set_ylim(0, 1)
    
    # x축 레이블과 눈금 폰트 사이즈 및 굵기 설정
    ax.set_xlabel(column, fontsize=14, fontweight='heavy')  # x축 레이블 폰트 사이즈 및 굵기 설정

    plt.xlabel(column)
    plt.ylabel('Proportion')
    
    plt.show()
    return fig

def chart_category(df):
    # 유니크 값이 10개 이하인 컬럼 필터링
    unique_value_columns = [col for col in df.columns if df[col].nunique() <= 10]

    if not unique_value_columns:
        print("No columns with 10 or fewer unique values.")
        return

    # 위젯 정의
    column_widget = widgets.Dropdown(
        options=unique_value_columns,
        value=unique_value_columns[0],
        description='Column:',
        style={'description_width': 'initial'},
        layout=Layout(width='70%')
    )

    output = widgets.Output()
    save_button = widgets.Button(description="Export PNG", layout=Layout(width='70%'))

    def update_graph(*args):
        with output:
            clear_output(wait=True)
            fig = plot_stacked_bar_chart(df, column_widget.value)
            save_button.on_click(lambda b: fig.savefig(f"{column_widget.value}.png"))

    def on_button_click(b):
        fig = plot_stacked_bar_chart(df, column_widget.value)
        fig.savefig(f"{column_widget.value}.png")

    # 컬럼 선택 시 그래프 업데이트
    column_widget.observe(update_graph, names='value')
    save_button.on_click(on_button_click)

    # 초기 그래프 생성
    update_graph()

    # 입력란과 그래프를 각각 세로로 가운데 정렬
    input_widgets = VBox(
        [column_widget, save_button],
        layout=Layout(align_items='flex-start', width='25%')
    )

    # UI 구성
    ui = HBox(
        [input_widgets, output], 
        layout=Layout(align_items='center')
    )
    display(ui)