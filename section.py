import copy

from .material import Material
from .config import _INDENT


class Section:

    def __init__(self, name, ring=0, pitch=0., height=None, eqMethod="homo") -> None:
        """
        A Section of a segement in assembly

        Input
        -----
        name: str, the name of Section
        ring: int, the ring number of rod
        pitch: float, the pitch of rod
        bounds: (float, float), the coordinate height of the lower and higher bound
        height: float, the height of this Section (which should have been de-coupled with Section)
        eqMethod: str, the equivalence method in cross-section calculation: ["homo", "1D", "supercell"]
        scSection: Section, the adjacent section (ONLY ONE) of surrounding assembly in supercell
        """
        self.name = name
        self.id = -1
        self.ring = ring
        self.pitch = pitch
        self.height = height
        self.bounds = (0., 0.)
        self.rod = []
        self.region = []
        self.eqMethod = eqMethod
        self.scSection = None

        # Auxiliary variables
        self.r = -1
        self.k = -1

    def append(self, shape, size, material) -> None:
        """
        Append a geometry into this Section

        Input
        -----
        shape: str, "rod"/"circ"/"c" or "region"/"hex"/"h"
        size: float, the characteristic size
        material: class Material
        """
        # Convert the type of size into float
        if type(size) is not float:
            size = float(size)
        
        # Check material
        if type(material) is not Material:
            raise TypeError("Input material is not type Material. material={}".format(material))
        
        # Modify bounds if given
        # if type(bounds) is tuple or type(bounds) is list:
        #     bounds = tuple(bounds)
        #     if len(bounds) != 2:
        #         raise ValueError("Input bounds has NOT 2 elements.")
        #     if type(bounds[0]) is not float or type(bounds[1]) is not float:
        #         raise TypeError("Input bounds are type {} & {}, which should be float.".format(type(bounds[0], type(bounds[1]))))
        #     self.bounds = bounds

        # Check shape & append data
        if type(shape) is not str:
            raise TypeError("Input shape is not type str. shape={}".format(shape))
        shape = shape.lower()
        if shape == "rod" or shape == "circ" or shape == "c":
            shape = "rod"
            self.rod.append((size, material))
        elif shape == "region" or shape == "hex" or shape == "h":
            shape = "region"
            self.region.append((size, material))
        else:
            raise ValueError("Invalid shape: {}".format(shape))
    
    def appendRod(self, size, material) -> None:
        """
        Append a rod geometry into this Section
        size: float, the characteristic size
        material: class Material
        """
        self.append(shape="rod", size=size, material=material)

    def appendRegion(self, size, material) -> None:
        """
        Append a region geometry into this Section
        size: float, the characteristic size
        material: class Material
        """
        self.append(shape="region", size=size, material=material)
    
    def check(self) -> bool:
        """
        Check this Section model
        """
        isPass = True

        if self.ring != len(self.rod):
            isPass = False
            print("Warning: dismatch of ring_num and rod_num in Section [{}]".format(self.name))

        return isPass

    def toTULIP(self, sc=False) -> str:
        """
        Convert the Section data to TULIP format

        Input
        -----
        sc: bool, add "sc_" to the head of every keyword of marginal assembly in supercell
        """
        # Check Section ID & ring
        if self.id == -1:
            raise RuntimeError("Section {} ID not specified.".format(self.name))
        if self.ring == -1 and len(self.rod) > 0:
            raise RuntimeError("Ring number not specified in Section {}.".format(self.name))
        
        # Add "sc_" to supercell keyword
        keywordHeader = "sc_" if sc else ""

        if not sc:
            slist = ["! Section: {}".format(self.name), "mat{:d}".format(self.id)]
        else:
            slist = []
        
        # Rod
        if len(self.rod) > 0:
            slist.append("{}rod".format(keywordHeader))
            slist.append("{}{:<14}{:d}".format(_INDENT, "{}rod_num".format(keywordHeader), self.ring))
            rod_geo = "{}{:<14}".format(_INDENT, "{}rod_geo".format(keywordHeader))
            rod_mat = "{}{:<14}".format(_INDENT, "{}rod_mat".format(keywordHeader))
            for size, material in self.rod:
                rod_geo += "{:<8.4f}".format(size)
                rod_mat += "{:<8d}".format(material.id)
            slist.append(rod_geo)
            slist.append(rod_mat)
            # Keyword "rod_dis" is not needed in supercell
            if not sc:
                slist.append("{}{:<14}{:<8.4f}".format(_INDENT, "rod_dis", self.pitch))
            if not sc and self.eqMethod != "supercell":
                slist.append("{}{:<14}{}".format(_INDENT, "ring_num", \
                    " ".join(list([str(i) for i in range(1, self.ring+1)]))))
        
        # Region
        if len(self.region) > 0:
            slist.append("{}region   hex".format(keywordHeader))
            region_geo = "{}{:<14}".format(_INDENT, "{}region_geo".format(keywordHeader))
            region_mat = "{}{:<14}".format(_INDENT, "{}region_mat".format(keywordHeader))
            for size, material in self.region:
                region_geo += "{:<8.4f}".format(size)
                region_mat += "{:<8d}".format(material.id)
            slist.append(region_geo)
            slist.append(region_mat)
        
        # Remove the keyword "end" in central section, folloed by "sc_rod"
        if self.eqMethod != "supercell":
            slist.append("end")
            slist.append('\n')
        
        return '\n'.join(slist)

    def copy(self, name, deepcopy=False):
        """
        Create a copy of this Section

        Input
        -----
        name: str, a new name for this copy, use "@ORIGIN" to represent the original name
        deepcopy: bool, whether replace copy.copy() with copy.deepcopy()

        Example
        -------
        ```python
        >>> origin = Section(name="fuel")
        >>> copy1 = origin.copy(name="copy")
        >>> copy1.name
        "copy"
        >>> copy2 = origin.copy(name="@ORIGIN_copy")
        >>> copy2.name
        "fuel_copy"
        ```
        """
        # Check input
        try:
            name = str(name)
        except ValueError:
            raise TypeError("Input name could not be converted to string")
        
        newName = name.replace("@ORIGIN", self.name)
        if newName == self.name:
            raise ValueError("Names conflict between Section and its copy.")

        # Create copy
        # Maybe the copy could be simplified
        sectionCopy = Section(name=newName, ring=self.ring, pitch=self.pitch)
        sectionCopy.id = -1
        sectionCopy.bounds = self.bounds
        sectionCopy.height = self.height
        if deepcopy:
            sectionCopy.rod = copy.deepcopy(self.rod)
            sectionCopy.region = copy.deepcopy(self.region)
        else:
            sectionCopy.rod = copy.copy(self.rod)
            sectionCopy.region = copy.copy(self.region)
        return sectionCopy

    def __str__(self) -> str:
        slist = ["Section"]
        slist.append("{:<15}{}".format("name:", self.name))
        slist.append("{:<15}{:d}".format("ID:", self.id))
        slist.append("{:<15}{:d}".format("ring:", self.ring))
        slist.append("{:<15}{:.4f}".format("pitch: ", self.pitch))
        slist.append("{:<15}{:<10.4f}{:<10.4f}".format("bounds:", *self.bounds))
        slist.append("{:<15}{:d}".format("rod:", len(self.rod)))
        slist.append("{:<15}{:d}".format("region:", len(self.region)))

        return '\n'.join(slist)

