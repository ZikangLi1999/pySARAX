import numpy as np
import pandas as pd

from .config import _PERIODIC_TABLE, _INDENT, _BULK_TEMP, _GROUP_STRUCTURE


class Material:

    def __init__(self, name) -> None:
        """
        Create a Material object with name and composition
        
        Input
        -----
        name: str, the identical name of this Material
        """
        self.name = name
        self.id = -1
        self.composition = pd.DataFrame(columns=("Nuclide", "Density", "Temperature"))
        self.removeLowDensity = True
    
    def addNuclide(self, nuclide, density, temperature) -> None:
        """
        Add a nuclide into this Material
        nuclide: str, the nuclide symbol, like "Pu239"
        density: float, the atom density of nuclide in UNIT "atom/(cm*barn)"
        temperature: float, the temperature of material in UNIT "K" 
        """
        # Check input
        if type(nuclide) is not str:
            raise TypeError("Input nuclide is not type str. nuclide={}".format(nuclide))
        try:
            density = float(density)
        except ValueError:
            raise TypeError("Cannot convert density to float type. density={}".format(density))
        try:
            temperature = float(temperature)
        except ValueError:
            raise TypeError("Cannot convert temperature to float type. temperature={}".format(temperature))
        
        self.composition = self.composition.append(\
            {"Nuclide": nuclide, "Density": density, "Temperature": temperature}, ignore_index=True)

    def addElement(self, element, density, temperature) -> None:
        """
        Add an element into this Material
        element: str, the element symbol, like "Fe"
        density: float, the atom density of nuclide in UNIT "atom/(cm*barn)"
        temperature: float, the temperature of material in UNIT "K" 
        """
        # Check input
        if type(element) is not str:
            raise TypeError("Input element is not type str. element={}".format(element))
        try:
            density = float(density)
        except ValueError:
            raise TypeError("Cannot convert density to float type. density={}".format(density))
        try:
            temperature = float(temperature)
        except ValueError:
            raise TypeError("Cannot convert temperature to float type. temperature={}".format(temperature))
        
        # Convert element to nuclide based on natural enrichment
        nuclides = _PERIODIC_TABLE.elem2nuc(element=element, density=density)
        for nuc, den in nuclides:
            self.composition = self.composition.append(\
                {"Nuclide": nuc, "Density": den, "Temperature": temperature}, ignore_index=True)

    def fromCSV(self, csvPath, defaultTemp=_BULK_TEMP) -> None:
        """
        Read material composition from a CSV file
        csvPath: pathLike
        defaultTemp: float, the default temperature of material, only used when "Temperature" is missing in CSV
        """
        self.composition = pd.read_csv(csvPath)

        # If the nuclide is given in ZAIDS format, 
        #     traverse every row of the raw data, then convert the ZAIDS to isotope symbol
        if "ZAIDS" in self.composition:
            self.composition["ZAIDS"] = self.composition["ZAIDS"].apply(lambda row: self.parseZAIDS(int(row)))
            self.composition.columns = ("Nuclide", "Density")
        
        # If temperature not given, set it to defaultTemp
        if "Temperature" not in self.composition:
            self.composition["Temperature"] = \
                defaultTemp * np.ones(shape=len(self.composition), dtype=float)

    def fromDataFrame(self, df, defaultTemp=_BULK_TEMP) -> None:
        """
        Read material composition from a pandas.DataFrame
        df: pandas.DataFrame
        defaultTemp: float, the default temperature of material, only used when "Temperature" is missing in DataFrame
        """
        # Check df type
        if type(df) is not pd.DataFrame:
            raise TypeError("Input df is not type pandas.DataFrame. df={}".format(df))
        
        self.composition = df.copy()

        # Convert ZAIDS to nuclide if necessary
        if "ZAIDS" in self.composition:
            self.composition["ZAIDS"] = self.composition["ZAIDS"].apply(lambda row: self.parseZAIDS(int(row)))
        
        # If temperature not given, set it to defaultTemp
        if "Temperature" not in self.composition:
            self.composition["Temperature"] = \
                defaultTemp * np.ones(shape=len(self.composition), dtype=float)

        # Change the column title
        self.composition.columns = ("Nuclide", "Density", "Temperature")

    def check(self) -> bool:
        """
        Check the data of material composition before converting to TULIP

        Checking Items
        --------------
        - Remove the nuclides with density lower than 1E-15, which may lead to NaN when calculating the backgroud cross-section
        """
        if self.removeLowDensity:
            self.composition = self.composition[self.composition['Density'] >= 1E-15]
        
        return True

    def toTULIP(self) -> str:
        """
        Convert the material composition to TULIP format
        """
        # Check Material ID
        if self.id == -1:
            raise RuntimeError("Material ID not specified.")
        
        # Check Material Composition
        if not self.check():
            raise RuntimeError("Material check NOT passed in m{}".format(self.id))
        
        slist = ["! material: {}".format(self.name), "m{}".format(self.id)]
        # Traverse every row(nuclide), convert the data to TULIP format:
        # m(id)
        #     (nuclide)  (Density)  (Temperature)
        for idx, row in self.composition.iterrows():
            slist.append("{}{:<9}{:>12.3E}{:>10.2f}".format(\
                _INDENT, row["Nuclide"], row["Density"], row["Temperature"]))
        
        slist.append('\n')
        return '\n'.join(slist)

    def parseZAIDS(self, zaids):
        """
        Convert the ZAIDS to isotope symbol
        zaids: str or int, like "96242"

        Remark: This method is weak-coupled, which can be moved out of this class
        Update: Use PeriodicTable.zaids2nuc() instead. This method will be kept for compatibility
        """
        if type(zaids) == str:
            protonNumber = int(zaids[:-3])
            massNumber = int(zaids[-3:])
        if type(zaids) == float:
            zaids = int(zaids)
            protonNumber = zaids // 1000
            massNumber = zaids % 1000
        if type(zaids) == int:
            protonNumber = zaids // 1000
            massNumber = zaids % 1000
        # elementName = _PERIODIC_TABLE[_PERIODIC_TABLE['Nr'] == protonNumber]['Symbol'].item()
        elementName = _PERIODIC_TABLE.proton2elem(proton=protonNumber)
        return ''.join((elementName, str(massNumber)))

    def __str__(self) -> str:
        slist = ["Material"]
        slist.append("{:<15}{}".format("name:", self.name))
        slist.append("{:<15}{:d}".format("ID:", self.id))
        slist.append("{:<15}{:d} nuclides".format("composition:", len(self.composition)))
        return '\n'.join(slist)
    
    def __add__(self, other):
        # Find out all of the nuclides in this & other
        allNuclides = np.unique(np.append(self.composition['Nuclide'].values, other.composition['Nuclide'].values))
        # print(self.composition)
        # print(other.composition)

        result = Material(name='{}+{}'.format(self.name, other.name))
        getItem = lambda df: df.item() if not df.empty else 0.0
        for nuc in allNuclides:
            thisNuc = self.composition[self.composition['Nuclide'] == nuc]
            otherNuc = other.composition[other.composition['Nuclide'] == nuc]
            # print(nuc, thisNuc['Density'], otherNuc['Density'])
            result.addNuclide(
                nuc,
                getItem(thisNuc['Density']) + getItem(otherNuc['Density']),
                (getItem(thisNuc['Temperature']) + getItem(otherNuc['Temperature'])) / sum([1 - df.empty for df in (thisNuc['Temperature'], otherNuc['Temperature'])])
            )
        
        return result
    
    def __rmul__(self, other):
        result = Material(name='{:.4E}*({})'.format(other, self.name))
        result.composition = self.composition.copy()
        result.composition['Density'] = other * result.composition['Density']
        return result
