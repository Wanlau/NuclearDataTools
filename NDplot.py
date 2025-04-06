import copy
import json
import math
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
from matplotlib.patches import Patch

## the synthesis methods: Mass Spectroscopy, Radioactive Decay, Light Particles, Fission, Fusion, Spallation, Projectile Fragmentation, and Transfer/Deep Inelastic Scattering
colors_synthesisMethods = {
    "MS" : (0, 0, 0),
    "RD" : (0, 255, 255),
    "LP" : (255, 165, 0),
    "FI" : (255, 255, 0),
    "FU" : (255, 0, 0),
    "SP" : (0, 0, 255),
    "PF" : (0, 127, 0),
    "UN" : (127, 0, 127)
}

names_synthesisMethods = {
    "MS" : "Mass Spectroscopy",
    "RD" : "Radioactive Decay",
    "LP" : "Light Particles",
    "FI" : "Fission",
    "FU" : "Fusion",
    "SP" : "Spallation",
    "PF" : "Projectile Fragmentation",
    "UN" : "Transfer/Deep Inelastic"
}

colors_halflife = {
    "l100ns": (247, 189, 222),
    "100ns" : (255, 198, 165),
    "1us"   : (255, 231, 198),
    "10us"  : (255, 255, 156),
    "100us" : (255, 255, 16),
    "1ms"   : (231, 247, 132),
    "10ms"  : (214, 239, 57),
    "100ms" : (173, 222, 99),
    "1s"    : (82, 181, 82),
    "10s"   : (99, 189, 181),
    "100s"  : (99, 198, 222),
    "1ks"   : (0, 165, 198),
    "10ks"  : (8, 154, 148),
    "100ks" : (0, 132, 165),
    "10Ms"  : (49, 82, 165),
    "1e10s" : (41, 0, 107),
    "1e15s" : (0, 0, 0),
    "ST"    : (0, 0, 0),
    "SU"    : (255, 148, 115),
    "UN"    : (224, 224, 224)
}

legends_halflife = {
    "l100ns": "<100ns",
    "100ns" : "100ns ~ 1us",
    "1us"   : "1us ~ 10us",
    "10us"  : "10us ~ 100us",
    "100us" : "100us ~ 1ms",
    "1ms"   : "1ms ~ 10ms",
    "10ms"  : "10ms ~ 100ms",
    "100ms" : "100ms ~ 1s",
    "1s"    : "1s ~ 10s",
    "10s"   : "10s ~ 100s",
    "100s"  : "100s ~ 1ks",
    "1ks"   : "1ks ~ 10ks",
    "10ks"  : "10ks ~ 100ks",
    "100ks" : "100ks ~ 10Ms",
    "10Ms"  : "10Ms ~ 1e10s",
    "1e10s" : "1e10s ~ 1e15s",
    "1e15s" : ">1e15s",
    "ST"    : "STABLE",
    "SU"    : "SpecialUnit",
    "UN"    : "UNKNOWN"
}

colors_DecayModes = {
    "P"     : (255, 148, 115),
    "N"     : (156, 123, 189),
    "A"     : (255, 255, 66),
    "B-"    : (231, 140, 198),
    "EC+B+" : (99, 198, 222),
    "EC"    : (0, 132, 165),
    "SF"    : (82, 181, 82),
    "STABLE": (0, 0, 0),
    "UNKNOWN":(224, 224, 224)
}

legends_DecayModes = {
    "P"     : "Proton",
    "N"     : "Neutron",
    "A"     : "Alaph",
    "B-"    : "Beta-",
    "EC+B+" : "Beta+ ElectronCapture",
    "EC"    : "ElectronCapture",
    "SF"    : "SpontaneousFission",
    "STABLE": "STABLE",
    "UNKNOWN":"UNKNOWN"
}

nuclides_data_path = "data/nndc_nudat_data_export.json"
data_synthesisMethods_path = "data/Nuclides_synthesisMethods.json"
ElementsList_path = "data/ElementsList.json"
data_NuclidesClassifiedHalflife_path = "data/NuclidesClassifiedHalflife.json"
data_NuclidesClassifiedDecayModes_path = "data/NuclidesClassifiedDecayModes.json"

## 使用matplotlib绘图
## 入参为核素分类模式、所绘制核素区域、显示信息、有无图例
## area部分待完善
def nucildesChartPlotPLT(plot_mode=0, area=((0,0),(118,177)), text_mode=0, have_legend=True):
    ## 确定绘制核素区域
    z_min, n_min = area[0]
    z_max, n_max = area[1]
    area_ysize = z_max - z_min + 1
    area_xsize = n_max - n_min + 1

    ## 计算坐标时使用的单位为像素，而matplotlib的方法使用的单位为英寸，执行前须进行单位转换
    nbwidth = 40
    nbheight = 40
    dpi = 100

    ## 根据显示信息确定绘图大小
    if text_mode == 0:
        resize_ratio = 5
        ticks_step = 10
        fontsize_label = 75
        fontsize_ticks = 50
    elif text_mode == 1:
        resize_ratio = 5
        ticks_step = 10
        fontsize_label = 75
        fontsize_ticks = 50
    elif text_mode == 2:
        resize_ratio = 10
        ticks_step = 10
        fontsize_label = 150
        fontsize_ticks = 100
    elif text_mode == 3:
        resize_ratio = 10
        ticks_step = 10
        fontsize_label = 150
        fontsize_ticks = 100
    else:
        resize_ratio = 1
        ticks_step = 10
        fontsize_label = 15
        fontsize_ticks = 10

    ## 确定绘制范围
    ## 据质子数、中子数取整十
    ## 对于默认的(118,177)，绘制范围为(130, 200)
    x_min = math.floor(n_min / 10) * 10
    x_max = math.ceil(n_max / 10) * 10
    y_min = math.floor(z_min / 10) * 10
    y_max = math.ceil(z_max / 10) * 10

    fig = plt.figure(figsize=(10*resize_ratio*(x_max-x_min)/200, 6*resize_ratio*(y_max-y_min)/130), dpi=dpi)
    ax = plt.gca()


    ## 核素区域绘制，默认白色背景
    bg_color = (255, 255, 255)
    color_data = np.full((area_ysize, area_xsize, 3), bg_color)

    dx = nbwidth  * resize_ratio
    dy = nbheight * resize_ratio
    xposg = np.arange(0, area_xsize*dx + dx, dx)
    yposg = np.arange(0, area_ysize*dy + dy, dy)

    xgrid = np.tile(xposg.reshape(1, -1), (area_ysize + 1, 1))
    ygrid = np.tile(yposg.reshape(-1, 1), (1, area_xsize + 1))

    xgridI = xgrid / dpi
    ygridI = ygrid / dpi

    ## 根据核素分类模式上色
    color_data = nucildesChartPlotPLTColor(copy.deepcopy(color_data), plot_mode, z_min, n_min ,z_max, n_max)

    ## 停用imshow转用pcolormesh以便控制各网格大小以及绘制边框
    #ax.imshow(color_data, origin="lower", extent=[n_min-0.5, n_max+0.5, z_min-0.5, z_max+0.5])
    ax.pcolormesh(xgridI, ygridI, color_data.astype(np.uint8), edgecolors="white", linewidth=4*resize_ratio/dpi)


    ## 绘制坐标轴
    ax.set_xlim(((x_min-5)*dx/dpi, (x_max+20)*dx/dpi))
    ax.set_ylim(((y_min-5)*dy/dpi, (y_max+10)*dy/dpi))
    ax.set_xlabel("Neutron number", fontsize=fontsize_label, font='Times New Roman')
    ax.set_ylabel("Proton number" , fontsize=fontsize_label, font='Times New Roman')

    ## 刻度绘制，同样取整十 (ticks_step = 10)
    xticks = np.arange(x_min, x_max + 20 + ticks_step, ticks_step)
    yticks = np.arange(y_min, y_max + 10 + ticks_step, ticks_step)
    xtickspos = (xticks * dx + np.ones(len(xticks)) * 0.5 * dx) / dpi
    ytickspos = (yticks * dy + np.ones(len(yticks)) * 0.5 * dy) / dpi
    ax.set_xticks(xtickspos, xticks, fontsize=fontsize_ticks)
    ax.set_yticks(ytickspos, yticks, fontsize=fontsize_ticks)

    ## 图例添加
    ## 默认大小匹配全图绘制
    if have_legend:
        legend_handles = legendHandlesGet(plot_mode)
        if len(legend_handles) < 13:
            fontsize_legend = fontsize_ticks
        else:
            fontsize_legend = fontsize_ticks * 0.6
        plt.legend(handles=legend_handles, loc="lower right", fontsize=fontsize_legend)

    ## 添加显示信息
    if not text_mode == 0:
        color_weight = np.array((0.299, 0.587, 0.114))
        text_data = nucildesChartPlotPLTText(plot_mode, z_min, n_min ,z_max, n_max)

        ## 显示元素名称
        if text_mode == 1:
            for tdata in text_data:
                ## 根据方块颜色选择文本颜色(黑或白)
                base_color = color_data[tdata["pos"][1]][tdata["pos"][0]]
                if np.dot(color_weight, base_color) > 128:
                    text_color = "black"
                else:
                    text_color = "white"

                txposi = (tdata["pos"][0] + 0.5) * dx / dpi
                typosi = (tdata["pos"][1] + 0.5) * dy / dpi
                ax.text(txposi, typosi, tdata["text01"], ha="center", va="center", color=text_color, fontsize=6)

        ## 显示核素名称
        elif text_mode == 2:
            for tdata in text_data:
                ## 根据方块颜色选择文本颜色(黑或白)
                base_color = color_data[tdata["pos"][1]][tdata["pos"][0]]
                if np.dot(color_weight, base_color) > 128:
                    text_color = "black"
                else:
                    text_color = "white"

                txposi = (tdata["pos"][0] + 0.5) * dx / dpi
                typosi = (tdata["pos"][1] + 0.5) * dy / dpi
                ax.text(txposi, typosi, tdata["text02"], ha="center", va="center", color=text_color, fontsize=6)

        ## 显示详细信息
        elif text_mode == 3:
            for tdata in text_data:
                ## 根据方块颜色选择文本颜色(黑或白)
                base_color = color_data[tdata["pos"][1]][tdata["pos"][0]]
                if np.dot(color_weight, base_color) > 128:
                    text_color = "black"
                else:
                    text_color = "white"

                txposi = (tdata["pos"][0] + 0.5) * dx / dpi
                typosi = (tdata["pos"][1] + 0.85) * dy / dpi
                text = tdata["text02"] + "\n" + tdata["text03"]
                ax.text(txposi, typosi, text, ha="center", va="top", ma="center", color=text_color, fontsize=2.4)
                
        
    
    #fig.savefig("test.svg", format="svg")
    #fig.savefig("test.png", format="png")
    #plt.show()

    return fig

def nucildesChartPlotPLTColor(color_data, mode, z_min, n_min ,z_max, n_max):
    ## 据半衰期上色
    ## 默认使用基态数据
    ## nndc上的nudat3绘制时若基态无数据则会使用激发态的数据，此处与其不同 比如：137Pm、154Lu、161Ta 等
    if mode == 0:
        ## 自NDfilter.nuclidesClassifyHalflife()
        with open(data_NuclidesClassifiedHalflife_path, "r", encoding="utf8") as file:
            data = json.load(file)

        for row in data:
            if (row["z"] >= z_min and row["z"] <= z_max) and (row["n"] >= n_min and row["n"] <= n_max):
                ypos = row["z"] - z_min
                xpos = row["n"] - n_min
                if row["type"] in colors_halflife:
                    color_data[ypos][xpos] = np.array(colors_halflife[row["type"]])

    ## 据衰变模式上色
    elif mode == 1:
        ## 自NDfilter.nuclidesClassifyDecayMode()
        with open(data_NuclidesClassifiedDecayModes_path, "r", encoding="utf8") as file:
            data = json.load(file)

        for row in data:
            if (row["z"] >= z_min and row["z"] <= z_max) and (row["n"] >= n_min and row["n"] <= n_max):
                ypos = row["z"] - z_min
                xpos = row["n"] - n_min
                if row["type"] in colors_DecayModes:
                    color_data[ypos][xpos] = np.array(colors_DecayModes[row["type"]])
                elif row["type"] in ("2B-", "β⁻"):
                    color_data[ypos][xpos] = np.array(colors_DecayModes["B-"])
                elif row["type"] in ("2P", "3P"):
                    color_data[ypos][xpos] = np.array(colors_DecayModes["P"])
                elif row["type"] in ("2N"):
                    color_data[ypos][xpos] = np.array(colors_DecayModes["N"])
                else:
                    color_data[ypos][xpos] = np.array(colors_DecayModes["UNKNOWN"])

    ## 据合成方法上色
    elif mode == 2:
        with open(data_synthesisMethods_path, "r", encoding="utf8") as file:
            data = json.load(file)

        for row in data:
            if (row["z"] >= z_min and row["z"] <= z_max) and (row["n"] >= n_min and row["n"] <= n_max):
                ypos = row["z"] - z_min
                xpos = row["n"] - n_min
                if row["type"] in colors_synthesisMethods:
                    color_data[ypos][xpos] = np.array(colors_synthesisMethods[row["type"]])

    return color_data

def legendHandlesGet(plot_mode):
    handles = []
    if plot_mode == 0:
        for hl_tag in colors_halflife:
            if hl_tag == "ST":
                pass
            elif hl_tag == "1e15s":
                handles.append(Patch(facecolor=np.array(colors_halflife[hl_tag])/255., label=">1e15s or Stable"))
            else:
                handles.append(Patch(facecolor=np.array(colors_halflife[hl_tag])/255., label=legends_halflife[hl_tag]))

    elif plot_mode == 1:
        for decay_mode in colors_DecayModes:
            handles.append(Patch(facecolor=np.array(colors_DecayModes[decay_mode])/255., label=legends_DecayModes[decay_mode]))

    elif plot_mode == 2:
        for synthesis_method in colors_synthesisMethods:
            handles.append(Patch(facecolor=np.array(colors_synthesisMethods[synthesis_method])/255., label=names_synthesisMethods[synthesis_method]))

    return handles

def nucildesChartPlotPLTText(plot_mode, z_min, n_min ,z_max, n_max):
    with open(ElementsList_path, "r", encoding="utf8") as file:
        elements_list = json.load(file)

    text_data = []
    if plot_mode == 0 or plot_mode == 1:
        with open(nuclides_data_path, "r", encoding="utf8") as file:
            data = json.load(file)

        ## 填充半衰期及衰变模式信息
        ## 对于衰变模式多于三种的，只显示其前三种(仿nndc)
        ## 对于过长的数字，将会对其进行截断
        for nom, ndata in data.items():
            if (ndata["z"] >= z_min and ndata["z"] <= z_max) and (ndata["n"] >= n_min and ndata["n"] <= n_max):
                ypos = ndata["z"] - z_min
                xpos = ndata["n"] - n_min
                text01 = elements_list[str(ndata["z"])]
                text02 = str(ndata["z"]+ndata["n"]) + text01

                hlt = ""
                dmst = []
                if len(ndata["levels"]) == 0:
                    hlt = ""
                    dmst = []
                else:
                    if not "halflife" in ndata["levels"][0]:
                        hlt = ""
                    elif not "value" in ndata["levels"][0]["halflife"]:
                        hlt = ""
                    elif ndata["levels"][0]["halflife"]["value"] == "STABLE":
                        hlt = "STABLE"
                    else:
                        hlv = ndata["levels"][0]["halflife"]["value"]
                        if hlv > 1e4:
                            hlt = f"{hlv:.2e}"
                        else:
                            hlt = str(hlv)

                        if ndata["levels"][0]["halflife"]["unit"] == "m":
                            hlt = hlt + " min"
                        else:
                            hlt = hlt + " " + ndata["levels"][0]["halflife"]["unit"]

                    if not "decayModes" in ndata["levels"][0]:
                        dmst = []
                    else:
                        decay_modes = ndata["levels"][0]["decayModes"]["observed"] + ndata["levels"][0]["decayModes"]["predicted"]
                        if len(decay_modes) == 0:
                            dmst = []
                        else:
                            ## 据分支比排序
                            dms1 = []
                            dms2 = []
                            dmsd = {}
                            dmsd1 = []
                            for decay_mode in decay_modes:  
                                if not "value" in decay_mode:
                                    dms2.append(decay_mode)
                                else:
                                    dmsd[decay_mode["mode"]] = decay_mode
                                    dmsd1.append((decay_mode["value"], decay_mode["mode"]))
                            dmsdf = pd.DataFrame(dmsd1, columns=["value", "mode"])
                            dmsdf.sort_values(by="value", ascending=False, inplace=True)
                            for mode in dmsdf["mode"]:
                                dms1.append(dmsd[mode])
                            dms = dms1 + dms2
                            ###
                            for decay_mode in dms:
                                text = decay_mode["mode"]
                                if not "value" in decay_mode:
                                    text = text + " ?"
                                else:
                                    dmv = decay_mode["value"]
                                    dmt = str(dmv)
                                    if len(dmt) > 7:
                                        if re.search(rf"([eE])", dmt) == None:
                                            if dmv > 1e-3:
                                                dmt = dmt[:7]
                                            else:
                                                match00 = re.fullmatch(rf"(0+)(\.)(0+)([1-9]+)", dmt)
                                                if not match00 == None:
                                                    ne = match00.group(4)
                                                    if len(ne) < 3:
                                                        dmt = f"{dmv:e}"
                                                    else:
                                                        dmt = f"{dmv:.2e}"
                                    if decay_mode["uncertainty"]["type"] == "limit":
                                        if decay_mode["uncertainty"]["limitType"] == "lower":
                                            if decay_mode["uncertainty"]["isInclusive"] == True:
                                                text = text + " ≥ " + dmt + "%"
                                            else:
                                                text = text + " > " + dmt + "%"
                                        elif decay_mode["uncertainty"]["limitType"] == "upper":
                                            if decay_mode["uncertainty"]["isInclusive"] == True:
                                                text = text + " ≤ " + dmt + "%"
                                            else:
                                                text = text + " < " + dmt + "%"
                                    else:
                                        text = text + " = " + dmt + "%"
                                dmst.append(text)
                text03 = hlt + "\n"
                if not len(dmst) > 0:
                    pass
                elif len(dmst) <= 3:
                    for dmt in dmst:
                        text03 = text03 + "\n" + dmt
                else:
                    for idx in range(0, 3):
                        text03 = text03 + "\n" + dmst[idx]
                
                text_data.append({"pos":(xpos, ypos), "text01":text01, "text02":text02, "text03":text03})

    ## 填充合成方法信息
    elif plot_mode == 2:       
        with open(data_synthesisMethods_path, "r", encoding="utf8") as file:
            data = json.load(file)

        for row in data:
            if (row["z"] >= z_min and row["z"] <= z_max) and (row["n"] >= n_min and row["n"] <= n_max):
                ypos = row["z"] - z_min
                xpos = row["n"] - n_min
                text01 = elements_list[str(row["z"])]
                text02 = str(row["z"]+row["n"]) + text01
                text03 = row["type"]
                text_data.append({"pos":(xpos, ypos), "text01":text01, "text02":text02, "text03":text03})
                
    return text_data


## test
