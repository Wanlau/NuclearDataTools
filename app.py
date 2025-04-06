import gradio as gr
import pandas as pd
import json
from tempfile import NamedTemporaryFile

import NDfilter
import NDplot

## 核素筛选
def process_filters(Z_min, Z_max, Z_oe_idx, N_min, N_max, N_oe_idx, A_min, A_max, A_oe_idx, hl_enable_idx, hl_min, hl_min_unit, hl_max, hl_max_unit, dm_enable_idx, decay_modes):

    filtered_data = NDfilter.nuclidesFilterZNA(nuclides_data, Z_min, Z_max, Z_oe_idx, N_min, N_max, N_oe_idx, A_min, A_max, A_oe_idx)

    # 根据母核半衰期进行筛选
    if hl_enable_idx == 1:
        hl_min_sec = HLunit_convert(hl_min, hl_min_unit)
        hl_max_sec = HLunit_convert(hl_max, hl_max_unit)
        filtered_data = NDfilter.nuclidesFilterHalflife(filtered_data, hl_min_sec, hl_max_sec)

    # 根据衰变模式进行筛选
    if dm_enable_idx > 0 and decay_modes:
        filtered_data = NDfilter.nuclidesFilterDecayModes(filtered_data, dm_enable_idx, decay_modes)

    # 结果处理
    if len(filtered_data) == 0:
        result_text = "没有找到符合条件的核素"
        result_file_path = None
    else:
        result_text = f"找到 {len(filtered_data)} 个符合条件的核素"
        with NamedTemporaryFile(suffix=".json", delete=False, mode='w') as file:
            json.dump(filtered_data, file, indent=2)
            result_file_path = file.name
    return result_text, result_file_path

## 核素查找
def process_search(mode_idx, nom, z, n, a, preview_mode, file_type):
    result = None
    if mode_idx == 0:
        result = NDfilter.nuclidesSearchingNom(nuclides_data, nom.replace(" ", ""))
    elif mode_idx == 1:
        if not (z == None or n == None):
            result = NDfilter.nuclidesSearchingZN(nuclides_data, z, n)
    elif mode_idx == 2:
        if not (z == None or a == None):
            result = NDfilter.nuclidesSearchingZA(nuclides_data, z, a)
    elif mode_idx == 3:
        if not (n == None or a == None):
            result = NDfilter.nuclidesSearchingNA(nuclides_data, n, a)

    ## 结果处理
    if result == None:
        result_text = "没有找到此核素"
        result_dataframe = None
        result_file_path = None
    else:
        ## 文本
        name = result["name"]
        result_text = f"{name}"
        if len(result["levels"]) == 0:
            result_text = result_text + "\n此核素无数据"
        result_text = result_text + "\n\nnndc页面：\n\ngetdataset:\n" + f"https://www.nndc.bnl.gov/nudat3/getdataset.jsp?nucleus={name}&unc=NDS"
        ## temp
        with open ("data/haveDecayPage.json",'r', encoding='utf-8') as file:
            haveDecayPage = json.load(file)
        if haveDecayPage[name]:
            result_text = result_text + "\n\ndecaysearchdirect:\n" + f"https://www.nndc.bnl.gov/nudat3/decaysearchdirect.jsp?nuc={name}&unc=NDS"
        
        ## 预览表格
        if preview_mode == 0:
            result_dataframe = NDfilter.nuclideData_dict2dataframeCompact(result)
        elif preview_mode == 1:
            result_dataframe = NDfilter.nuclideData_dict2dataframe(result)

        ## 文件
        if file_type == 0:
            with NamedTemporaryFile(suffix=".json", delete=False, mode='w') as file:
                json.dump(result, file, indent=2)
                result_file_path = file.name
        elif file_type == 1:
            if preview_mode == 1:
                tmpDataframe = result_dataframe
            else:
                tmpDataframe = NDfilter.nuclideData_dict2dataframe(result)
            with NamedTemporaryFile(suffix=".csv", delete=False, mode='w') as file:
                tmpDataframe.to_csv(file, index=False)
                result_file_path = file.name
        else:
            result_file_path = None

    return result_text, result_dataframe, result_file_path

## 核素图绘制
def process_plot(plot_mode, text_mode, have_legend_idx, file_type3, using_filter, Z_min=0, Z_max=118, N_min=0, N_max=177):
    if have_legend_idx == 0:
        have_legend = True
    else:
        have_legend = False

    area = ((0,0),(118,177))
    if using_filter == 1:
        area = ((Z_min,N_min), (Z_max, N_max))

    result_fig = NDplot.nucildesChartPlotPLT(plot_mode, area, text_mode, have_legend)

    with NamedTemporaryFile(suffix=".svg", delete=False, mode='wb') as file:
        result_fig.savefig(file, format="svg")
        preview_svg_path = file.name

    if file_type3 == "png":
        with NamedTemporaryFile(suffix=".png", delete=False, mode='wb') as file:
            result_fig.savefig(file, format="png")
            result_file_path = file.name
    elif file_type3 == "svg":
        result_file_path = preview_svg_path
    else:
        result_file_path = None

    return result_fig, result_file_path


## 半衰期单位转换
def HLunit_convert(hl, hl_unit):
    if hl_unit == "Stable":
        hl_sec = None
    else:
        hl_sec = hl * HL_UNITS[hl_unit]
    return hl_sec

## 处理查找页面输入组件激活情况
def update_inputs2(mode_idx):
    
    if mode_idx == 0:  # 核素名称模式
        return [gr.Textbox(interactive=True), gr.Number(interactive=False, value=None), gr.Number(interactive=False, value=None), gr.Number(interactive=False, value=None)]
    elif mode_idx == 1:  # Z+N模式
        return [gr.Textbox(interactive=False, value=None), gr.Number(interactive=True), gr.Number(interactive=True), gr.Number(interactive=False, value=None)]
    elif mode_idx == 2:  # Z+A模式
        return [gr.Textbox(interactive=False, value=None), gr.Number(interactive=True), gr.Number(interactive=False, value=None), gr.Number(interactive=True)]
    elif mode_idx == 3:  # N+A模式
        return [gr.Textbox(interactive=False, value=None), gr.Number(interactive=False, value=None), gr.Number(interactive=True), gr.Number(interactive=True)]
    else:
        return [gr.Textbox(interactive=False, value=None), gr.Number(interactive=False, value=None), gr.Number(interactive=False, value=None), gr.Number(interactive=False, value=None)]


def update_inputs3(using_filter):
    if using_filter == 0:
        return [gr.Number(interactive=False, value=None)] * 4
    elif using_filter == 1:
        return [gr.Number(interactive=True)] * 4
    else:
        return [gr.Number(interactive=False, value=None)] * 4

## 导入数据集
nuclides_data_path = "data/nndc_nudat_data_export.json"
with open (nuclides_data_path,'r', encoding='utf-8') as file:
    nuclides_data = json.load(file)

## 半衰期单位转换字典
HL_UNITS = {"fs": 1e-15, "ps": 1e-12, "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400, "y": 31557600, "ky": 31557600e3, "My": 31557600e6, "Gy": 31557600e9}

with gr.Blocks(title="核数据工具") as demo:
    gr.Markdown("""
                ## 核数据工具
                可能是用来处理核数据的相关工具？？

                目前功能有：核素筛选、核素查找、核素图绘制。
                """)
    
    with gr.Tab("核素筛选"):
        gr.Markdown("""
                    ## 核素筛选
                    可以通过质子数(Z)、中子数(N)、质量数(A)以及母核半衰期、衰变模式等进行筛选
                    """)
        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("根据质子数(Z)、中子数(N)、质量数(A)进行筛选")
            with gr.Column(scale=4):
                with gr.Row():
                    with gr.Column(min_width=240):
                        gr.Markdown("质子数(Z)")
                        Z_min = gr.Number(label="最小值", precision=0)
                        Z_max = gr.Number(label="最大值", precision=0)
                        Z_oe = gr.Dropdown(["任意", "奇Z", "偶Z"], label="奇偶", type="index", interactive=True)
                    with gr.Column(min_width=240):
                        gr.Markdown("中子数(N)")
                        N_min = gr.Number(label="最小值", precision=0)
                        N_max = gr.Number(label="最大值", precision=0)
                        N_oe = gr.Dropdown(["任意", "奇N", "偶N"], label="奇偶", type="index", interactive=True)
                    with gr.Column(min_width=240):
                        gr.Markdown("质量数(A)")
                        A_min = gr.Number(label="最小值", precision=0)
                        A_max = gr.Number(label="最大值", precision=0)
                        A_oe = gr.Dropdown(["任意", "奇A", "偶A"], label="奇偶", type="index", interactive=True)

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("根据母核半衰期进行筛选")
            with gr.Column(scale=4):
                with gr.Row():
                    hl_enable = gr.Radio(["不使用", "使用"], value="不使用", type="index", show_label=False)
                    hl_min = gr.Number(label="最小值", minimum=0.)
                    hl_min_unit = gr.Dropdown(["fs", "ps", "ns", "us", "ms", "s", "m", "h", "d", "y", "ky", "My", "Gy", "Stable"], value="fs", interactive=True)
                    hl_max = gr.Number(label="最大值", minimum=0.)
                    hl_max_unit = gr.Dropdown(["fs", "ps", "ns", "us", "ms", "s", "m", "h", "d", "y", "ky", "My", "Gy", "Stable"], value="Stable", interactive=True)

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("根据衰变模式进行筛选")
            with gr.Column(scale=4):
                with gr.Row():
                    dm_enable_idx = gr.Radio(["不使用", "筛选包含所有以下所选衰变模式的核素(and)", "筛选包含任意以下所选衰变模式的核素(or)"], value="不使用", type="index", interactive=True, show_label=False)
                with gr.Row():
                    decayModes = gr.CheckboxGroup(['B-', 'N', '2N', 'B-N', 'P', 'B-A', 'B-2N', 'B-3N', '2P', 'EC', 'A', 'B-4N', 'EC+B+', 'ECA', 'ECP', 'IT', 'EC2P', 'EC3P', 'ECAP', '3P', '2B-', 'ECSF', '14C', 'B-SF', '24NE', 'SF', '20O', '20NE', '25NE', '28MG', 'NE', '22NE', 'SI', 'MG', '34SI'], label="decayModes", interactive=True)

        with gr.Row():
            submit_btn = gr.Button("筛选", variant="primary")
            reset_btn = gr.Button("重置条件", variant="primary")
        
        with gr.Row():
            result_text = gr.Textbox(label="筛选结果", interactive=False, show_copy_button=True)
            result_file = gr.File(label="结果文件", interactive=False)

        inputs = [
            Z_min, Z_max, Z_oe,
            N_min, N_max, N_oe,
            A_min, A_max, A_oe,
            hl_enable, hl_min, hl_min_unit, hl_max, hl_max_unit,
            dm_enable_idx, decayModes
        ]
        
        submit_btn.click(
            fn=process_filters,
            inputs=inputs,
            outputs=[result_text, result_file]
        )
        
        reset_btn.click(
            fn=lambda: [None,None,"任意"]*3 + ["不使用", None, "fs", None, "Stable"] + ["不使用", []],
            outputs=inputs
        )


    with gr.Tab("核素查找"):
        gr.Markdown("## 核素查找")
        with gr.Row():
            with gr.Column(scale=3):
                searchingMode = gr.Radio(["核素名称", "质子数(Z)、中子数(N)", "质子数(Z)、质量数(A)", "中子数(N)、质量数(A)"], value="核素名称", label="查找模式", interactive=True, type="index")
                gr.Markdown("### 根据核素名称查找")
                nuclide_in = gr.Textbox(value=None, info="请输入由质量数及元素名称所组成的核素名称，示例：232Th、232TH、th232、232-Th、th-232等。\n只要不太离谱就能识别……大概？")
                gr.Markdown("### 根据质子数(Z)、中子数(N)、质量数(A)查找")
                with gr.Row():
                    Z_in = gr.Number(value=0, label="质子数(Z)", precision=0, interactive=False)
                    N_in = gr.Number(value=0, label="中子数(N)", precision=0, interactive=False)
                    A_in = gr.Number(value=0, label="质量数(A)", precision=0, interactive=False)
                previewMode = gr.Radio(["紧凑", "常规"], value="紧凑", type="index", label="预览模式")
                outputFileType = gr.Radio(["json", "csv"], value="json", type="index", label="导出文件格式")
                with gr.Row():
                    submit_btn2 = gr.Button("查找", variant="primary")
                    reset_btn2 = gr.Button("重置条件", variant="primary")

            with gr.Column(scale=2):
                result_text2 = gr.Textbox(interactive=False, show_label=False)
                preview_df = gr.Dataframe(label="数据预览", interactive=False)
                result_file2 = gr.File(interactive=False)

        searchingMode.change(
            fn=update_inputs2,
            inputs=searchingMode,
            outputs=[nuclide_in, Z_in, N_in, A_in]
        )

        inputs2 = [
            searchingMode,
            nuclide_in,
            Z_in, N_in, A_in,
            previewMode, outputFileType
        ]

        submit_btn2.click(
            fn=process_search,
            inputs=inputs2,
            outputs=[result_text2, preview_df, result_file2]
        )

        reset_btn2.click(
            fn=lambda: [None]*4,
            outputs=[nuclide_in, Z_in, N_in, A_in]
        )

    with gr.Tab("核素图绘制"):
        gr.Markdown("""
                    ## 核素图绘制

                    可根据半衰期、衰变模式、合成方法等分类模式绘制核素图。

                    部分代码参考了Ming-Hao-Zhang的[Nuclei-Chart-Generator](https://github.com/Ming-Hao-Zhang/Nuclei-Chart-Generator)
                    """)
        with gr.Row():
            with gr.Column(scale=3):
                img_preview = gr.Plot(label="图片预览")
            with gr.Column(scale=2):
                plot_mode = gr.Radio(["寿命", "衰变模式", "合成方法"], value="寿命", label="核素分类模式", type="index")
                text_mode = gr.Radio(["无", "元素名称", "核素名称", "详细信息"], value="无", label="显示信息", type="index")
                have_legend_idx = gr.Radio(["显示图例", "隐藏图例"], value="显示图例", label="图例", type="index")
                file_type3 = gr.Radio(["svg", "png"], value="svg", label="导出格式", info="png为位图格式，svg为矢量图格式。\n显示详细信息时，受分辨率限制，png格式将会失真。如需高清晰度图像，请使用svg。")

                with gr.Accordion(open=False, label="更多选项") as filter3:
                    gr.Markdown("根据质子数(Z)、中子数(N)筛选 (未启用)")
                    using_filter = gr.Radio(["不使用", "使用"], value="不使用", show_label=False, type='index', interactive=False)
                    with gr.Row():
                        with gr.Column(min_width=120):
                            gr.Markdown("质子数(Z)")
                            Z_min = gr.Number(label="最小值", precision=0, interactive=False)
                            Z_max = gr.Number(label="最大值", precision=0, interactive=False)
                        with gr.Column(min_width=120):
                            gr.Markdown("中子数(N)")
                            N_min = gr.Number(label="最小值", precision=0, interactive=False)
                            N_max = gr.Number(label="最大值", precision=0, interactive=False)

                with gr.Row():
                    submit_btn3 = gr.Button("绘制", variant="primary")
                    reset_btn3 = gr.Button("重置条件", variant="primary")

                result_file3 = gr.File(interactive=False)

        using_filter.change(
            fn=update_inputs3,
            inputs=using_filter,
            outputs=[Z_min, Z_max, N_min, N_max]
        )

        filter3.collapse(
            fn=lambda: "不使用",
            outputs= using_filter
        )

        inputs3 = [plot_mode, text_mode, have_legend_idx, file_type3, using_filter, Z_min, Z_max, N_min, N_max]

        submit_btn3.click(
            fn=process_plot,
            inputs=inputs3,
            outputs=[img_preview, result_file3]
        )

        reset_btn3.click(
            fn=lambda: ["寿命", "无", "显示图例", "svg", "不使用"] + [None]*4,
            outputs=inputs3
        )

    demo.launch()