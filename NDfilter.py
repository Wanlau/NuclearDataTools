import json
import os
import re
import copy
import pandas as pd

## 半衰期单位转换字典
HL_UNITS = {"fs": 1e-15, "ps": 1e-12, "ns": 1e-9, "us": 1e-6, "ms": 1e-3, "s": 1, "m": 60, "h": 3600, "d": 86400, "y": 31557600, "ky": 31557600e3, "My": 31557600e6, "Gy": 31557600e9}
nuclides_data_path = "data/nndc_nudat_data_export.json"

def nuclidesFilterZNA(nuclides_data, Z_min=None, Z_max=None, Z_oe_idx=0, N_min=None, N_max=None, N_oe_idx=0, A_min=None, A_max=None, A_oe_idx=0):
    filtered = {}
    for name, data in nuclides_data.items():
        fh = True
        if not (Z_min is None or data['z'] >= Z_min):
            fh = False
        if not (Z_max is None or data['z'] <= Z_max):
            fh = False
        if not (N_min is None or data['n'] >= N_min):
            fh = False
        if not (N_max is None or data['n'] <= N_max):
            fh = False
        if not (A_min is None or data['a'] >= A_min):
            fh = False
        if not (A_max is None or data['a'] <= A_max):
            fh = False
        
        if Z_oe_idx == 1:
            if data['z'] % 2 == 0:
                fh = False
        elif Z_oe_idx == 2:
            if data['z'] % 2 == 1:
                fh = False

        if N_oe_idx == 1:
            if data['n'] % 2 == 0:
                fh = False
        elif N_oe_idx == 2:
            if data['n'] % 2 == 1:
                fh = False

        if A_oe_idx == 1:
            if data['a'] % 2 == 0:
                fh = False
        elif A_oe_idx == 2:
            if data['a'] % 2 == 1:
                fh = False
        
        if fh:
            filtered[name] = data

    return filtered
        
def nuclidesFilterHalflife(nuclides_data, hl_min_sec=0, hl_max_sec=None):
    filtered = {}
    if hl_min_sec == None:
        if hl_max_sec == None:
            for name, data in nuclides_data.items():
                if "levels" in data:
                    for level in data["levels"]:
                        if "halflife" in level:
                            if level["halflife"]["value"] == "STABLE":
                                filtered[name] = data
                                break
    else:
        if hl_max_sec == None:
            for name, data in nuclides_data.items():
                if "levels" in data:
                    for level in data["levels"]:
                        if "halflife" in level:
                            if level["halflife"]["value"] == "STABLE":
                                filtered[name] = data
                                break
                            else:
                                if not level["halflife"]["unit"] in HL_UNITS:
                                    break
                                hl_sec = level["halflife"]["value"] * HL_UNITS[level["halflife"]["unit"]]
                                if hl_sec > hl_min_sec:
                                    filtered[name] = data
                                    break
        else:
            for name, data in nuclides_data.items():
                if "levels" in data:
                    for level in data["levels"]:
                        if "halflife" in level:
                            if level["halflife"]["value"] == "STABLE":
                                pass
                            else:
                                if not level["halflife"]["unit"] in HL_UNITS:
                                    break
                                hl_sec = level["halflife"]["value"] * HL_UNITS[level["halflife"]["unit"]]
                                if hl_sec > hl_min_sec and hl_sec < hl_max_sec:
                                    filtered[name] = data
                                    break

    return filtered

def nuclidesFilterDecayModes(nuclides_data, dm_enable_idx, decay_modes):
    filtered = {}
    nuclides_decayModes = {}
    for name, data in nuclides_data.items():
        decayModes = []
        if "levels" in data:
            for level in data["levels"]:
                if "decayModes" in level:
                    for decayData in level["decayModes"]["observed"]:
                        if not decayData["mode"] in decayModes:
                            decayModes.append(decayData["mode"])
        nuclides_decayModes[name] = copy.deepcopy(decayModes)

    decay_modes_series = pd.Series(decay_modes)

    ## nndc导出的数据里有两处写作"β⁻"了
    replace_dict = {"β⁻": "B-"}

    for name, decayModes in nuclides_decayModes.items():
        nuclide_decayModes_series = pd.Series(decayModes)
        nuclide_decayModes_series = nuclide_decayModes_series.replace(replace_dict)

        if dm_enable_idx == 1:
            if decay_modes_series.isin(nuclide_decayModes_series).all():
                filtered[name] = nuclides_data[name]
        elif dm_enable_idx == 2:
            if decay_modes_series.isin(nuclide_decayModes_series).any():
                filtered[name] = nuclides_data[name]

    return filtered

def nuclidesSearchingNom(nuclides_data, nuclide_nom):
    nuclide = None
    element = None
    nuclide_A = None
    if not re.fullmatch(rf"([A-Za-z]+)([-_|]*)([0-9]+)", nuclide_nom) == None:
        match00 = re.fullmatch(rf"([A-Za-z]+)([-_|]*)([0-9]+)", nuclide_nom)
        element = match00.group(1)
        nuclide_A = match00.group(3)
    elif not re.fullmatch(rf"([0-9]+)([-_|]*)([A-Za-z]+)", nuclide_nom) == None:
        match00 = re.fullmatch(rf"([0-9]+)([-_|]*)([A-Za-z]+)", nuclide_nom)
        element = match00.group(3)
        nuclide_A = match00.group(1)

    if not (element == None or nuclide_A == None):
        if not element in ("n", "N"):
            if len(element) == 1:
                element = element.upper()
            else:
                element = element[:1].upper() + element[1:].lower()
        while nuclide_A.startswith("0"):
            nuclide_A = nuclide_A[1:]
        nom = nuclide_A + element
        if nom in nuclides_data:
            nuclide = nuclides_data[nom]
    
    return nuclide

def nuclidesSearchingZN(nuclides_data, Z_in=0, N_in=0):
    nuclide = None
    for nom, data in nuclides_data.items():
        if data["z"] == Z_in and data["n"] == N_in:
            nuclide = data
            break

    return nuclide

def nuclidesSearchingZA(nuclides_data, Z_in=0, A_in=0):
    nuclide = None
    for nom, data in nuclides_data.items():
        if data["z"] == Z_in and data["a"] == A_in:
            nuclide = data
            break

    return nuclide

def nuclidesSearchingNA(nuclides_data, N_in=0, A_in=0):
    nuclide = None
    for nom, data in nuclides_data.items():
        if data["n"] == N_in and data["a"] == A_in:
            nuclide = data
            break

    return nuclide

def nuclideData_dict2dataframe(nuclideData_dict):
    nom = nuclideData_dict["name"]
    z = nuclideData_dict["z"]
    n = nuclideData_dict["n"]
    a = nuclideData_dict["a"]
    data = []
    for level in nuclideData_dict["levels"]:
        ## 自旋宇称
        spinParity = None
        if "spinParity" in level:
            spinParity = level["spinParity"]

        ## 质量过剩
        massExcess = None
        massExcess_unit = None
        massExcess_unc =  None
        if not "massExcess" in level:
            pass
        else:
            massExcess = level["massExcess"]["unit"]
            massExcess_unit = level["massExcess"]["value"]
            massExcess_unc = level["massExcess"]["uncertainty"]

        ## 半衰期
        hl = None
        hl_unit = None
        hl_unc =  None
        if not "halflife" in level:
            pass
        elif level["halflife"]["value"] == "STABLE":
            pass
        else:
            hl = level["halflife"]["value"]
            hl_unit = level["halflife"]["unit"]
            if level["halflife"]["uncertainty"]["type"] == "asymmetric":
                hl_unc = "+" + str(level["halflife"]["uncertainty"]["upperLimit"]) + " -" + str(level["halflife"]["uncertainty"]["lowerLimit"])
            elif level["halflife"]["uncertainty"]["type"] == "approximation":
                hl_unc = "≈"
            elif level["halflife"]["uncertainty"]["type"] == "limit":
                if level["halflife"]["uncertainty"]["limitType"] == "lower":
                    if level["halflife"]["uncertainty"]["isInclusive"] == True:
                        hl_unc = "≥"
                    else:
                        hl_unc = ">"
                elif level["halflife"]["uncertainty"]["limitType"] == "upper":
                    if level["halflife"]["uncertainty"]["isInclusive"] == True:
                        hl_unc = "≤"
                    else:
                        hl_unc = "<"
            else:
                hl_unc = level["halflife"]["uncertainty"]["value"]

        ## 衰变模式
        if not "decayModes" in level:
            data.append((nom, z, n, a, level["energy"]["value"], level["energy"]["unit"], spinParity, massExcess, massExcess_unit, massExcess_unc, hl, hl_unit, hl_unc, None, None))
        elif len(level["decayModes"]["observed"]) == 0:
            data.append((nom, z, n, a, level["energy"]["value"], level["energy"]["unit"], spinParity, massExcess, massExcess_unit, massExcess_unc, hl, hl_unit, hl_unc, None, None))
        else:
            for decayMode in level["decayModes"]["observed"]:
                data.append((nom, z, n, a, level["energy"]["value"], level["energy"]["unit"], spinParity, massExcess, massExcess_unit, massExcess_unc, hl, hl_unit, hl_unc, decayMode["mode"], decayMode["value"]))

    nuclideDataframe = pd.DataFrame(data, columns=["Nuclide", "Z", "N", "A", "E(level)", "E(level) unit", "Spin Parity", "Mass Excess", "Mass Excess unit", "Mass Excess uncertainty", "Halflife", "Halflife unit", "Halflife uncertainty", "Decay Mode", "Branch Ratio"])

    return nuclideDataframe

def nuclideData_dict2dataframeCompact(nuclideData_dict):
    data = []
    for level in nuclideData_dict["levels"]:
        ## 自旋宇称
        spinParity = None
        if "spinParity" in level:
            spinParity = level["spinParity"]

        ## 质量过剩
        massExcess = None
        if not "massExcess" in level:
            pass
        else:
            massExcess = level["massExcess"]["formats"]["NDS"]+" "+level["massExcess"]["unit"]

        ## 半衰期
        hl = None
        if not "halflife" in level:
            pass
        elif level["halflife"]["value"] == "STABLE":
            hl = "STABLE"
        else:
            hl = level["halflife"]["formats"]["NDS"]+" "+level["halflife"]["unit"]

        ## 衰变模式
        if not "decayModes" in level:
            data.append((str(level["energy"]["value"])+" "+level["energy"]["unit"], spinParity, massExcess, hl, None))
        elif len(level["decayModes"]["observed"]) == 0:
            data.append((str(level["energy"]["value"])+" "+level["energy"]["unit"], spinParity, massExcess, hl, None))
        else:
            decayModes = []
            for decayMode in level["decayModes"]["observed"]:
                decayModes.append((decayMode["mode"], decayMode["value"]))
            decayModeDF = pd.DataFrame(decayModes, columns=["Decay Mode", "Branch Ratio"])
            data.append((str(level["energy"]["value"])+" "+level["energy"]["unit"], spinParity, massExcess, hl, decayModeDF))

    nuclideDataframe = pd.DataFrame(data, columns=["E(level)", "Spin Parity", "Mass Excess", "Halflife", "Decay Modes"])

    return nuclideDataframe

def nuclidesClassifyHalflife(data_path):
    with open (data_path,'r', encoding='utf-8') as file:
        nuclides_data = json.load(file)

    classified = []
    for nom, data in nuclides_data.items():
        z = data["z"]
        n = data["n"]
        hl_type = "UN"

        if len(data["levels"]) == 0:
            hl_type = "UN"
        elif not "halflife" in data["levels"][0]:
            hl_type = "UN"
        elif data["levels"][0]["halflife"]["value"] == "STABLE":
            hl_type = "ST"
        elif not data["levels"][0]["halflife"]["unit"] in HL_UNITS:
            hl_type = "SU"
        else:
            hl_sec = data["levels"][0]["halflife"]["value"] * HL_UNITS[data["levels"][0]["halflife"]["unit"]]
            if hl_sec < 1e-7:
                hl_type = "l100ns"
            elif hl_sec < 1e-6:
                hl_type = "100ns"
            elif hl_sec < 1e-5:
                hl_type = "1us"
            elif hl_sec < 1e-4:
                hl_type = "10us"
            elif hl_sec < 1e-3:
                hl_type = "100us"
            elif hl_sec < 1e-2:
                hl_type = "1ms"
            elif hl_sec < 1e-1:
                hl_type = "10ms"
            elif hl_sec < 1:
                hl_type = "100ms"
            elif hl_sec < 10:
                hl_type = "1s"
            elif hl_sec < 100:
                hl_type = "10s"
            elif hl_sec < 1e3:
                hl_type = "100s"
            elif hl_sec < 1e4:
                hl_type = "1ks"
            elif hl_sec < 1e5:
                hl_type = "10ks"
            elif hl_sec < 1e7:
                hl_type = "100ks"
            elif hl_sec < 1e10:
                hl_type = "10Ms"
            elif hl_sec < 1e15:
                hl_type = "1e10s"
            else:
                hl_type = "1e15s"
        
        classified.append({"z":z, "n":n, "type":hl_type})

    return classified

def nuclidesClassifyDecayMode(data_path):
    with open (data_path,'r', encoding='utf-8') as file:
        nuclides_data = json.load(file)

    classified = []
    for nom, data in nuclides_data.items():
        z = data["z"]
        n = data["n"]
        dm_type = "UNKNOWN"

        if len(data["levels"]) == 0:
            dm_type = "UNKNOWN"
        elif not "decayModes" in data["levels"][0]:
            dm_type = "UNKNOWN"
            if "halflife" in data["levels"][0]:
                if data["levels"][0]["halflife"]["value"] == "STABLE":
                    dm_type = "STABLE"
        else:
            decay_modes = data["levels"][0]["decayModes"]["observed"] + data["levels"][0]["decayModes"]["predicted"]
            if len(decay_modes) == 0:
                dm_type = "UNKNOWN"
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
                dm_type = dms[0]["mode"]
        classified.append({"z":z, "n":n, "type":dm_type})

    return classified

##test
