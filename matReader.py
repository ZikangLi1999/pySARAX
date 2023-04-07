"""
MAT reader for TULIP

Author: LZK
Date: 2023-3-3
"""
import os
import re
import json
import numpy as np

DATA_START_IDX = 16


class JsonNdarrayEncoder(json.JSONEncoder):

    def default(self, o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        return json.JSONEncoder.default(self, o)


class MatReader:

    def __init__(self) -> None:
        self.data = dict()

    def read(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError("MAT file does NOT exist: {}".format(path))
        
        with open(path, 'r', encoding='utf-8') as file:
            # Header information
            for i in range(5):
                line = file.readline()
                lineSplited = line.rstrip('\n').split()
                if not lineSplited:
                    continue
                self.data[lineSplited[0]] = lineSplited[1:]
            
            # Traverse every nuclide
            nuclide = str()
            matrixFlag = 0
            previousTitle = str()
            for line in file.readlines():
                # Split the line into Title & Data
                lineSplited = line.rstrip('\n').split()
                if not lineSplited:
                    continue
                title, data = lineSplited[0], lineSplited[1:]

                if 'name' in title:
                    nuclide, nucId, density, temperature = data
                    self.data[nuclide] = dict()
                    self.data[nuclide]['id'] = int(nucId)
                    self.data[nuclide]['density'] = float(density)
                    self.data[nuclide]['temperature'] = float(temperature)
                
                elif 'matrix' in title:
                    matrixFlag = int(self.data['groupnum'][0])
                    self.data[nuclide][title] = self._split(line).reshape((1, -1))
                    matrixFlag -= 1
                    previousTitle = title
                
                elif matrixFlag > 0:
                    self.data[nuclide][previousTitle] = np.concatenate((
                        self.data[nuclide][previousTitle],
                        self._split(line).reshape((1, -1))
                    ))
                    matrixFlag -= 1
                
                elif 'groupnum' in title:
                    self.data[nuclide][title] = np.array(data, dtype=np.int)
                
                elif 'Temp' in title:
                    self.data[nuclide][title] = np.array(data, dtype=np.int)
                
                elif 'Burn' in title:
                    self.data[nuclide][title] = np.array(data, dtype=float)

                else:
                    self.data[nuclide][title] = self._split(line)
    

    def _split(self, line) -> np.ndarray:
        rawData = line[DATA_START_IDX:].rstrip('\n')
        return np.array(re.findall(r'.{12}', rawData), dtype=np.float)


    def compare(self, other):
        assert type(other) is MatReader

        headerKeys = ('char4', 'depletion_chain', 'groupnum', 'Temp', 'Burn')
        
        # First, compare the nuclide types
        selfKeys = set(self.data.keys())
        otherKeys = set(other.data.keys())
        selfUnique = selfKeys.difference(otherKeys)
        otherUnique = otherKeys.difference(selfKeys)
        if selfUnique:
            print("This MAT has unique nuclides: {}".format(selfUnique))
        if otherUnique:
            print("The other MAT has unique nuclides: {}".format(otherUnique))
        
        result = dict()

        # Then compare the header
        for key in headerKeys:
            pass

        # Last, compare the xsec data
        allNuclides = tuple((selfKeys | otherKeys) - set(headerKeys))
        for nuclide in allNuclides:
            result[nuclide] = self._compareNuclide(this=self.data[nuclide], other=other.data[nuclide])

        ret = MatReader()
        ret.data = result

        return ret
    

    def _compareNuclide(self, this, other):
        allKeys = ('total', 'absorption', 'fission', 'nu_fission', 'chi', 'n2n', 'kappa_fission', 'scatter_matrixP0', 'scatter_matrixP1', 'chid', 'velocity', 'beta_total', 'decay_constants')

        ret = dict()
        for key in allKeys:
            refval = this[key]
            cmpval = other[key]
            relerr = np.abs(np.divide(cmpval - refval, refval, out=np.zeros_like(refval), where=(refval!=0)))
            ret[key] = relerr
        
        return ret
    

    def save(self, path):
        # Check the existence of path
        basePath, fileName = os.path.split(path)
        if not os.path.exists(basePath):
            os.mkdir(basePath)
        
        # Save self.data as a JSON file
        with open(path, 'w', encoding='utf-8') as file:
            json.dump(self.data, file, cls=JsonNdarrayEncoder, indent=4)


if __name__ == '__main__':
    from pprint import pprint
    path = r"C:\SJTUGraduate\Research\Projects\LoongSARAXVerif\jobs\tulip_mat49_20221129153811\output\xsec_decimal\MAT1"
    mat = MatReader()
    mat.read(path=path)
    pprint(mat.data['U235'])
