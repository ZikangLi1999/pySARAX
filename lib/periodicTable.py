"""
Periodic Table

Author: LZK
Date: 2023-1-22

Data Source: 
NIST https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some

File Log:
2023-1-22  File created
2023-1-23  PeriodicTable.elem2nuc(), PeriodicTable.proton2elem() completed
"""
import os
import pandas as pd

__author__ = "Zikang LI"
__source__ = "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some"
print(__file__)
__csvPath__ = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "dataset\\nist_natural_abundance.csv")

class PeriodicTable:

    def __init__(self) -> None:
        """
        Periodic Table based on NIST data
        """
        self.table = pd.read_csv(__csvPath__)

        # Amend data types
        self.table['Number'].astype(int)
        self.table['Mass'] = self.table['Mass'].astype(float)
        self.table['Composition'] = self.table['Composition'].astype(float)
    
    def elem2nuc(self, element, density) -> list:
        """
        Convert element density to nuclide density by multiplying natural enrichment
        
        Input
        -----
        element: str, the element symbol
        density: float, the element density

        Return
        ------
        A list containing all nuclides & their densities according to the input
        The format is like: [(nuc1, den1), (nuc2, den2), ...]

        Example
        -------
        ```python
        >>> table = PeriodicTable()
        >>> table.elem2nuc('H', 1.0)
        [('H1', 0.9997702030564883), ('H2', 0.00022979694351165497)]
        ```
        """
        data = self.table[self.table['Element'] == element].fillna(value=0.0)

        # Calculate the natural enrichment of nuclides
        data['Product'] = data['Mass'] * data['Composition']
        sumup = sum(data['Product'].values)
        data['Enrichment'] = data['Product'] / sumup
        data['Density'] = data['Enrichment'] * density

        # Convert to list format
        result = []
        for element, number, density in zip(data['Element'], data['Nuclide'], data['Density']):
            if density == 0.0:
                continue
            nuclide = ''.join((element, str(number)))
            result.append((nuclide, density))

        return result

    def proton2elem(self, proton) -> str:
        """
        Convert proton number to element symbol

        Input
        -----
        proton: int, the proton number

        Example
        -------
        >>> periodicTable().proton2elem(proton=1)
        'H'
        >>> periodicTable().proton2elem(proton=92)
        'U'
        >>> periodicTable().proton2elem(proton=94)
        'Pu'
        """
        assert type(proton) is int
        return self.table[self.table['Number'] == proton]['Element'].unique().item()

    def zaids2nuc(self, zaids) -> str:
        """
        Convert ZAIDS to nuclide symbol

        Input
        -----
        zaids: str or int, the ZAIDS of nuclide, like "92235" for U235

        Example
        -------
        >>> PeriodicTable().zaids2nuc(92235)
        'U235'
        >>> PeriodicTable().zaids2nuc(94239)
        'Pu239'
        """
        if type(zaids) is str:
            proton = int(zaids[:-3])
            mass = int(zaids[-3:])
        elif type(zaids) is int:
            proton = zaids // 1000
            mass = zaids % 1000
        else:
            raise TypeError("Input ZAIDS has wrong type {}".format(type(zaids)))
        
        return ''.join((self.proton2elem(proton=proton), str(mass)))

    @property
    def allElements(self) -> list:
        return list(self.table['Element'].unique())
    
    @property
    def allNuclides(self) -> list:
        return list(self.table[['Element', 'Nuclide']].astype(str).apply(''.join, axis=1))
    
    def __str__(self) -> str:
        return self.table.__str__()
    

# Test during development
if __name__ == '__main__':
    table = PeriodicTable()
    print(table.elem2nuc('H', 0.1))

    print(table.proton2elem(92))

    print(table.zaids2nuc(92235))
    print(table.zaids2nuc("94239"))
